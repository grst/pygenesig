import readline
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
from pygenesig.validation import SignatureTester
from pygenesig.tools import write_gmt
import sklearn.metrics
import numpy as np
from pygenesig.validation import SignatureGenerator
import tempfile
numpy2ri.activate()

base = importr("base")
stats = importr("stats")
edger = importr("edger")

ro.r('''
    model_matrix = function(target, covariates) {
        return(
            model.matrix(~ target + as.mcovariates
        )
    }
    ''')


def limma(expr, target, covariates):
    genes = np.array([str(i) for i in range(len(expr.shape[0]))])
    dg_list = edger.DGEList(counts=expr, genes=genes)
    design_mat = edger.model_matrix


class LimmaSignatureGenerator(SignatureGenerator):
    def __init__(self, expr, target, covariates, fold_change=100):
        """

        Args:
            expr:
            target:
            covariates:
            fold_change:
        """
        super(LimmaSignatureGenerator, self).__init__(expr, target)
        self.fold_change = fold_change
        self.covariates = covariates

    def mk_signatures(self, subset):
        df_aggr = aggregate_expression(self.expr[:, subset], self.target[subset])
        return None