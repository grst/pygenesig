import readline
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri, pandas2ri, Formula
from pygenesig.validation import SignatureTester
from pygenesig.tools import write_gmt
import sklearn.metrics
import numpy as np
from pygenesig.validation import SignatureGenerator
import tempfile
numpy2ri.activate()
pandas2ri.activate()

base = importr("base")
stats = importr("stats")
edge_r = importr("edgeR")
stats = importr("stats")
limma = importr("limma")


class LimmaSignatureGenerator(SignatureGenerator):
    TARGET_COL = "target"

    def __init__(self, expr, target, covariates, fold_change=100):
        """

        Args:
            expr:
            target:
            covariates (pd.DataFrame):
            fold_change:
        """
        super(LimmaSignatureGenerator, self).__init__(expr, target)
        assert covariates.shape[0] == expr.shape[1], "nrow(covariates) must equal ncol(expr)"
        self.fold_change = fold_change
        self.covariates = covariates

    def _mk_signatures(self, expr, target):
        # genenames = indices
        genes = np.array([str(i) for i in range(expr.shape[0])])
        dg_list = edge_r.DGEList(counts=expr, genes=genes)

        assert all([True if not col.startswith(self.TARGET_COL) else False for col in self.covariates.columns]), \
            "covariates must not contain a column named target. "
        fmla = Formula("~ target + " + " + ".join(self.covariates.columns))
        env = fmla.environment
        env[self.TARGET_COL] = self.target
        for col in self.covariates.columns:
            env[col] = self.covariates[col]
        design_mat = stats.model_matrix(fmla)

        cpm = edge_r.cpm(dg_list)
        cpm_cnt = np.sum(np.array(cpm) > 1, axis=1)
        cpm_inds = cpm_cnt >= 10

        expr_fil = expr[cpm_inds, :]
        genes_fil = genes[cpm_inds]

        # overwrite with filtered dglist
        dg_list = edge_r.DGEList(counts=expr_fil, genes=genes_fil)
        dg_list = edge_r.calcNormFactors(dg_list)

        # run voom
        voom_res = limma.voom(dg_list, design_mat, plot=False)

        # run limma
        lm_fit = limma.lmFit(voom_res, design_mat)

        # filter based on pvalue and logfoldchange
        p_adj = .01 / base.nrow(lm_fit)[0]
        fit_e = limma.treat(lm_fit, lfc=np.log2(self.fold_change))
        rslt = limma.decideTests(fit_e, p_value=p_adj)

        # make signatures
        colnames = np.array(base.colnames(rslt))
        np_rslt = np.array(rslt)
        tissue_cols = np.array([True if col.startswith('target') else False for col in colnames])
        np_rslt = np_rslt[:, tissue_cols]
        tissues = np.array([col[len(self.TARGET_COL):] for col in colnames[tissue_cols]])
        assert len(tissues) == np_rslt.shape[1]
        row_inds = np.where(cpm_inds)[0]
        signatures = {}
        for i, tissue in enumerate(tissues):
            signatures[tissue] = list(row_inds[np_rslt[:, i] == 1])

        return signatures