"""
pygenesig's main module: Contains the abstract base classes for signature
generation and testing and provides tools for cross-validation.
"""


from abc import abstractmethod, ABCMeta
import numpy as np
import sklearn.model_selection
from dask import delayed
import os.path
from collections import Counter


def filter_samples(exprs, target, n_splits):
    """
    Remove all categories and the respective samples
    from the input data that have not at least *n_splits*
    samples per category.

    This is useful when performing a cross-validation with
    *n* folds, as we need at least *n* samples per category.

    Args:
        exprs:
        target:
        n_splits:

    Returns:
        exprs:
        target:

    """
    target_cnt = Counter(target)
    keep_tissues = [tissue for tissue, count in target_cnt.items() if count >= n_splits]
    target_mask = np.in1d(target, keep_tissues)
    return exprs[:, target_mask], target[target_mask]


def cv_score(expr_file, target_file,
             signature_generator, signature_tester,
             splitter=sklearn.model_selection.StratifiedKFold(n_splits=10),
             sg_kwargs={}, st_kwargs={}):
    """
    Perform a crossvalidation by generating and testing signatures on a given expression matrix.

    Args:
        expr_file (str): path to numpy.array (binary) file containing expression matrix saved with ``numpy.save()``
        target_file (str): path to plaintext file containing target classes, line-by-line
        signature_generator (SignatureGenerator): ``SignatureGenerator`` used to derive the signatures
        signature_tester (SignatureTester): ``SignatureTester`` used to check the quality of signatures
        splitter (sklearn.model_selection._BaseKFold): crossvalidation method from `scikit-learn`_.

    Returns:
        (list, list, list, list): list of signature dictionaries containing the signatures generated
        on each fold, list of confusion matrices containing the confusion matrices calculated on each fold, 
        list of indices used for training for each fold, list of indices used for testing for each fold. 

    .. _scikit-learn:
        http://scikit-learn.org/stable/modules/classes.html#module-sklearn.model_selection

    """
    # we need the number of samples locally for splitting in test and training sets
    target_file = os.path.realpath(target_file)
    expr_file = os.path.realpath(expr_file)
    target_local = np.genfromtxt(target_file, dtype=str, delimiter=",")

    # the rest we can run delayed on the dask cluster
    signatures = []
    scores = []
    train_list = []
    test_list = []
    for train, test in splitter.split(list(enumerate(target_local)), target_local):
        # due to some bug in dask, it's faster to load all files on every worker separately.
        train_list.append(train)
        test_list.append(test)
        expr = delayed(np.load)(expr_file)
        target = delayed(np.genfromtxt)(target_file, dtype=str, delimiter=",")
        sg = delayed(signature_generator)(expr, target, **sg_kwargs)
        st = delayed(signature_tester)(expr, target, **st_kwargs)
        signature = delayed(sg.mk_signatures)(train)
        signatures.append(signature)
        score = delayed(st.score_signatures, nout=2)(signature, test)
        scores.append(score)
    sig_list = delayed(list)(signatures)
    res_list = delayed(list)(scores)
    return sig_list, res_list, train_list, test_list


class SignatureGenerator(metaclass=ABCMeta):
    """
    Abstract base-class for generating gene signatures.
    Child classes build gene signatures from a gene expression matrix and target labels.

    When initializing ``SignatureGenerator`` the input data is checked for consistency.
    You can override the ``__init__`` method in your implementation of ``SignatureGenerator``
    to pass additional parameters, however, make always sure to call this method
    in your ``__init__`` function to take advantage of the consistency check::

        def __init__(self, expr, target, your_param):
            super(YourSignatureGenerator, self).__init__(expr, target)
            # your code goes here


    Args:
        expr (np.ndarray): m x n matrix with m samples and n genes
        target (array-like): m-vector with true tissue for each sample
        
    .. automethod:: _mk_signatures
    """

    def __init__(self, expr, target):
        if not expr.shape[1] == len(target):
            raise ValueError("The length of target must equal the number of samples (columns) in the expr matrix. ")
        if not np.issubdtype(expr.dtype, int) and not np.issubdtype(expr.dtype, float):
            raise TypeError("Expression needs to be numeric. (dtype=int | dtype=float)")

        self.expr = expr
        self.target = target

    @abstractmethod
    def _mk_signatures(self, expr, target):
        """
        Implement this method.
        Generates signatures for a list of tissues given a gene expression matrix.
        This method is called internally by the public method `mk_signatures`, which
        takes care of subsetting the expression data.

        Args:
            expr: m x n matrix with m samples and n genes
            target: vector with true tissue for each sample

        Returns:
            dict: tissue_name -> [list, of, enriched, genes]. The indices correspond to expr.index.

        Note:
            When implementing this method, make sure, that every tissue in `this.target` has an entry in the
            dictionary, even when the signature does not contain any genes.

        """
        pass

    def mk_signatures(self, subset=None):
        """
        Make gene signatures on a subset of the samples in the gene expression matrix.

        Args:
            subset: Sample indices (columns of expression matrix) to use for signature generation.
                Useful for cross-validation. If omitted, all samples will be used.

        Returns:
            dict: signature dictionary. Maps target labels to a list of row-indices of the gene
            expression matrix.

            Example::

                signatures = {
                    'adipose tissue' : [
                        42,
                        124,
                        256,
                        1038,
                        ... ],
                    'skeletal muscle' : [
                        52,
                        679,
                        ... ],
                    ...
                }


        """
        if subset is None:
            subset = np.array(list(range(len(self.target))))

        return self._mk_signatures(self.expr[:, subset], self.target[subset])


class SignatureTester(metaclass=ABCMeta):
    """
    Abstract base-class for testing gene signatures.
    Child classes test if signatures are able to identify their respective samples,
    given an expression matrix and a list of the actual tissues.

    When initializing ``SignatureTester`` the input data is checked for consistency.
    You can override the ``__init__`` method in your implementation of ``SignatureTester``
    to pass additional parameters, however, make always sure to call this method
    in your ``__init__`` function to take advantage of the consistency check::

        def __init__(self, expr, target, your_param):
            super(YourSignatureTester, self).__init__(expr, target)
            # your code goes here

    Args:
        expr (np.ndarray): m x n gene expression matrix with m samples and n genes
        target (array-like): m-vector with true tissue for each sample (`ground truth`)
        
    .. automethod:: _score_signatures
    """

    def __init__(self, expr, target):
        if not expr.shape[1] == len(target):
            raise ValueError("The length of target must equal the number of samples (columns) in the expr matrix. ")
        if not np.issubdtype(expr.dtype, int) and not np.issubdtype(expr.dtype, float):
            raise TypeError("Expression needs to be numeric. (dtype=int | dtype=float)")
        self.expr = expr
        self.target = np.array(target)

    @staticmethod
    def sort_signatures(signatures):
        """
        Retrieve the signatures in a consistent order.

        The signature dictionary is a python default dict, which is unordered.
        However, when displaying results, it is desirable to display the
        signatures in a consistent order.

        You can use this method to iterate over the signature dictionary::

            for tissue in sort_signatures(signatures):
                print(tissue, signatures[tissue])

        Confusion matrices generated with ``SignatureTester.confusion_matrix`` are
        also sorted with ``sort_signatures``.

        Args:
            signatures: signature dictionary.

        Returns:
            list: sorted keys of the signature dictionary.

        """
        return sorted(signatures.keys())

    @abstractmethod
    def _score_signatures(self, expr, signatures):
        """
        Implement this method.

        Generates a score for each sample and signature.

        Args:
            expr (np.ndarray): m x n gene expression matrix with m genes and n samples
            signatures (dict of list): Signatures dictionary with j signatures.
                Dictionary: tissue_name -> [list, of, enriched, genes].

        Returns:
            np.ndarray: j x n score matrix with j signatures and n samples containing the
            scores generated by this method for each sample and signature. The rows are sorted by
            ``SignatureGenerator.sort_signatures``.

        .. important::
            Make sure the rows in the result matrix are sorted by ``SignatureGenerator.sort_signatures``.
            Make sure that a row for each signature in ``signatures`` is included in the matrix, even
            *if the signatures is empty* (number of genes = 0).

        """
        pass

    def score_signatures(self, signatures, subset=None):
        """
        Generates a score for each sample and signature. A high score represents
        a high representation of the signature in the sample.

        Args:
            signatures (dict of list): Signatures dictionary with j signatures.
                Dictionary: tissue_name -> [list, of, enriched, genes].
            subset: sample indices (columns of expression matrix) to use for testing. Useful for crossvalidation.

        Returns:
            np.ndarray: j x n score matrix with j signatures and n samples containing the
            scores generated by this method for each sample and signature. The rows are sorted by
            ``SignatureGenerator.sort_signatures``.

        """
        if all([len(genes) == 0 for genes in signatures.values()]):
            raise SignatureTesterException("Signature Set contains only empty signatures. ")

        if subset is None:
            # take all
            subset = np.array(list(range(len(self.target))))

        return self._score_signatures(self.expr[:, subset], signatures)

    def classify(self, signatures, score_matrix, subset=None):
        """
        Test signatures based on the expression matrix. Predicts the class labels
        using the signatures. A sample is predicted as the label associated with
        the highest scoring signature.

        Args:
            signatures (dict[str, list]): Signatures dictionary returned by ``SignatureGenerator.mk_signature``.
            score_matrix: j x n score matrix with j signatures and n samples generated by ``SignatureTester.score_signatures``
            subset: sample indices (columns of expression matrix) to use for testing. Useful for crossvalidation.

        Returns:
            (list, list): list of actual labels and predicted labels.

        """
        if subset is None:
            # take all
            subset = np.array(list(range(len(self.target))))

        actual = self.target[subset]  # standard of truth
        signature_list = self.sort_signatures(signatures)
        predicted = [signature_list[m] for m in np.nanargmax(score_matrix, axis=0)]

        return actual, predicted

    def confusion_matrix(self, signatures, actual, predicted):
        """
        Make a confusion matrix.

        This is a wrapper for the ``sklearn.metrics.confusion_matrix``. Makes the matrix contain all
        labels available in the signature.

        Note:
            The rows and columns are sorted with ``SignatureTester.sort_signatures``.
            You can therefore use ``sort_signatures`` to label you confusion matrix.

        Args:
            signatures: dict of signatures that was used to predict the labels
            actual: list of actual labels
            predicted: list of predicted labels

        Returns:
            numpy.matrix: confusion matrix

        """
        return sklearn.metrics.confusion_matrix(actual, predicted,
                                                labels=self.sort_signatures(signatures))

    def _test_signatures(self, signatures, subset):
        """
        Wrapper for score_signatures and classify. Used for unit testing.
        """
        scores = self.score_signatures(signatures, subset)
        return self.classify(signatures, scores, subset)


class SignatureTesterException(Exception):
    pass

