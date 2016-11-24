from abc import abstractmethod, ABCMeta
import numpy as np
import sklearn.model_selection
from dask import delayed
import os.path


def cross_validate_signatures(expr_file, target_file,
                              signature_generator, signature_tester,
                              splitter=sklearn.model_selection.StratifiedKFold(n_splits=10)):
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
        sg = delayed(signature_generator)(expr, target)
        st = delayed(signature_tester)(expr, target)
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
            expr (np.matrix): m x n matrix with m samples and n genes
            target (array-like): m-vector with true tissue for each sample
        """
        assert expr.shape[1] == len(target), \
            "The length of target must equal the number of samples (columns) in the expr matrix. "
        self.expr = expr
        self.target = target

    @abstractmethod
    def mk_signatures(self, subset):
        """
        Make gene signatures based on the expression matrix.

        Args:
            subset: Sample indices (columns of expression matrix) to use for signature generation.
                Useful for cross-validation.

        Returns:
            dict: tissue_name -> [list, of, enriched, genes]. The indices correspond to expr.index.

        Note:
            When implementing this method, make sure, that every in this.target has an entry in the
            dictionary, even when the signature does not contain any genes.

        """
        pass


class SignatureTester(metaclass=ABCMeta):
    """Abstract base-class for testing gene signatures.
    Child classes test if a signature is able to identify the respective tissue properly,
    given an expression matrix and a list of the actual tissues. """

    def __init__(self, expr, target):
        """
        Args:
            expr (np.matrix): m x n matrix with m samples and n genes
            target (array-like): m-vector with true tissue for each sample
        """
        assert expr.shape[1] == len(target), \
            "The length of target must equal the number of samples (columns) in the expr matrix. "
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
    def test_signatures(self, signatures, subset):
        """
        Test signatures based on the expression matrix.

        Args:
            signatures (dict of list): Signatures dictionary retured by SignatureGenerator.mk_signature.
                Dictionary: tissue_name -> [list, of, enriched, genes].
            subset: sample indices (columns of expression matrix) to use for testing. Useful for crossvalidation.

        Returns:
            np.matrix: Confusion Matrix

        """
        pass
