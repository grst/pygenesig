from abc import abstractmethod, ABCMeta
import numpy as np
import sklearn.model_selection
from dask import delayed
import os.path


def cross_validate_signatures(expr_file, target_file,
                              signature_generator, signature_tester,
                              splitter=sklearn.model_selection.StratifiedKFold(n_splits=10),
                              sg_kwargs={}, st_kwargs={}):
    """
    Perform a crossvalidation by generating and testing signatures on a given expression matrix.

    Args:
        expr_file (str): np matrix (binary) file containing expression matrix
        target_file (str): file containing target classes, line-by-line
        signature_generator (SignatureGenerator): SignatureGenerator used to derive the signatures
        signature_tester (SignatureTester): SignatureTester used to check the quality of signatures
        splitter (sklearn.model_selection._BaseKFold): crossvalidation method from scikit-learn.
            See [here](http://scikit-learn.org/stable/modules/classes.html#module-sklearn.model_selection).
    """
    # we need the number of samples locally for splitting in test and training sets
    target_file = os.path.realpath(target_file)
    expr_file = os.path.realpath(expr_file)
    target_local = np.genfromtxt(target_file, dtype=str, delimiter=",")

    # the rest we can run delayed on the dask cluster
    signatures = []
    results = []
    for train, test in splitter.split(list(enumerate(target_local)), target_local):
        # due to some bug in dask, it's faster to load all files on every worker separately.
        expr = delayed(np.load)(expr_file)
        target = delayed(np.genfromtxt)(target_file, dtype=str, delimiter=",")
        sg = delayed(signature_generator)(expr, target, **sg_kwargs)
        st = delayed(signature_tester)(expr, target, **st_kwargs)
        signature = delayed(sg.mk_signatures)(train)
        signatures.append(signature)
        results.append(delayed(st.test_signatures)(signature, test))
    sig_list = delayed(list)(signatures)
    res_list = delayed(list)(results)
    return sig_list, res_list


class SignatureGenerator(metaclass=ABCMeta):
    """Abstract base-class for generating gene signatures.
    Child classes build gene signatures for a given gene expression matrix and tissue list. """

    def __init__(self, expr, target):
        """
        Args:
            expr (np.ndarray): m x n matrix with m samples and n genes
            target (array-like): m-vector with true tissue for each sample
        """
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

        Args:
            expr: m x n matrix with m samples and n genes
            target: vector with true tissue for each sample

        Returns:
            dict: tissue_name -> [list, of, enriched, genes]. The indices correspond to expr.index.

        Note:
            When implementing this method, make sure, that every in this.target has an entry in the
            dictionary, even when the signature does not contain any genes.

        """
        pass

    def mk_signatures(self, subset=None):
        """
        Make gene signatures based on the expression matrix.

        Args:
            subset: Sample indices (columns of expression matrix) to use for signature generation.
                Useful for cross-validation.

        Returns:
            dict: tissue_name -> [list, of, enriched, genes]. The indices correspond to expr.index.

        """
        if subset is None:
            subset = np.array(list(range(self.target)))

        return self._mk_signatures(self.expr[:, subset], self.target[subset])


class SignatureTester(metaclass=ABCMeta):
    """Abstract base-class for testing gene signatures.
    Child classes test if a signature is able to identify the respective tissue properly,
    given an expression matrix and a list of the actual tissues. """

    def __init__(self, expr, target):
        """
        Args:
            expr (np.ndarray): m x n matrix with m samples and n genes
            target (array-like): m-vector with true tissue for each sample
        """
        if not expr.shape[1] == len(target):
            raise ValueError("The length of target must equal the number of samples (columns) in the expr matrix. ")
        if not np.issubdtype(expr.dtype, int) and not np.issubdtype(expr.dtype, float):
            raise TypeError("Expression needs to be numeric. (dtype=int | dtype=float)")
        self.expr = expr
        self.target = np.array(target)

    @staticmethod
    def sort_signatures(signatures):
        """
        Order the signatures. Same order as the rows/columns in the output matrix.

        Args:
            signatures: signature dictionary.

        Returns:
            list: sorted keys of the signature dictionary.

        Note: use this method when implementing test_signatures to make sure that
            the confusion matrix is sorted according to this order.

        """
        return sorted(signatures.keys())

    @abstractmethod
    def _predict(self, expr, signatures):
        """
        Implement this method.

        Predicts the tissue for each sample (=column) in the expression matrix
        given a set of signatures.

        Args:
            expr (np.array): gene expression matrix
            signatures (dict of list): Signatures dictionary returned by SignatureGenerator.mk_signature.
                Dictionary: tissue_name -> [list, of, enriched, genes].

        Returns:
            list: list of names of highest scoring signature (=predicted class) for each column in expr.

        """
        pass

    def test_signatures(self, signatures, subset=None):
        """
        Test signatures based on the expression matrix.

        Args:
            signatures (dict[str, list]): Signatures dictionary returned by SignatureGenerator.mk_signature.
                Dictionary: tissue_name -> [list, of, enriched, genes].
            subset: sample indices (columns of expression matrix) to use for testing. Useful for crossvalidation.

        Returns:
            (list, list): list of actual labels and predicted labels.

        """
        if all([len(genes) == 0 for genes in signatures.values()]):
            raise SignatureTesterException("Signature Set contains only empty signatures. ")

        if subset is None:
            # take all
            subset = np.array(list(range(len(self.target))))

        actual = self.target[subset]  # standard of truth
        predicted = self._predict(self.expr[:, subset], signatures)

        return actual, predicted

    def confusion_matrix(self, signatures, actual, predicted):
        """
        Make a confusion matrix.

        This is a wrapper for the sklearn.metrics.confusion matrix. Makes the matrix contain all
        labels available in the signature.

        Args:
            signatures: dict of signatures that was used to predict the labels
            actual: list of actual labels
            predicted: list of predicted labels

        Returns:
            np.matrix: confusion matrix

        """
        return sklearn.metrics.confusion_matrix(actual, predicted,
                                                labels=self.sort_signatures(signatures))


class SignatureTesterException(Exception):
    pass

