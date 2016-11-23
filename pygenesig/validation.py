from abc import abstractmethod, ABCMeta
import numpy as np
import sklearn.model_selection
from dask import delayed


def cross_validate_signatures(expr_file, target_file,
                              signature_generator, signature_tester,
                              splitter=sklearn.model_selection.StratifiedKFold(n_splits=10)):
        """

        Args:
            expr_file (str): csv file containing expression matrix
            target_file (str): csv file containin target classes
            signature_generator (SignatureGenerator): SignatureGenerator used to derive the signatures
            signature_tester (SignatureTester): SignatureTester used to check the quality of signatures
            splitter (sklearn.model_selection._BaseKFold): crossvalidation method from scikit-learn.
                See [here](http://scikit-learn.org/stable/modules/classes.html#module-sklearn.model_selection).
        """
        expr = delayed(np.genfromtext)(expr_file)
        target = delayed(np.genfromtext)(target_file)
        sg = delayed(signature_generator)(expr, target)
        st = delayed(signature_tester)(expr, target)
        signatures = []
        results = []
        for train, test in splitter:
            signatures.append(delayed(sg.mk_signature)(train))
            results.append(delayed(st.test_signature)(test))
        sig_list = delayed(list)(signatures)
        res_list = delayed(list)(results)
        return sig_list, res_list


class SignatureGenerator(metaclass=ABCMeta):
    def __init__(self, expr, target):
        """

        Args:
            expr (np.matrix): m x n matrix with m samples and n genes
            target (array-like): m-vector with true tissue for each sample
        """
        self.expr = expr
        self.target = target

    @abstractmethod
    def mk_signatures(self, subset):
        """

        Returns:
            dict: tissue_name -> [list, of, enriched, genes]. The indices correspond to expr.index.

        """
        pass


class SignatureTester(metaclass=ABCMeta):
    def __init__(self, expr, target):
        self.expr = expr
        self.target = target

    @abstractmethod
    def test_signatures(self, subset):
        """

        Returns:
            np.matrix: Confusion Matrix

        """
        pass
