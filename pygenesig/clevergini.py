"""
The Gini-coefficient (gini index) measures the inequality among values of
a frequency distribution. Originally applied to measure income inequality
we can also use it on gene expression data to find genes that are
over-represented in certain samples.
"""

import numpy as np
import pandas as pd
from pygenesig.validation import SignatureGenerator
from pygenesig.tools import collapse_matrix
from pygenesig.gini import gini


def get_gini_signatures(df_aggr, min_gini=0.7, max_rk=3, min_expr=1):
    """
    Generate gene signatures using gini index.

    Finds over-represented genes for each sample group, given the specified thresholds.

    Args:
        df_aggr (pd.DataFrame): m x n pandas DataFrame with m Genes and n tissues
        min_gini (float): gini cutoff, genes need to have a gini index larger than this value.
        max_rk (int): rank cutoff, include genes if they rank <= max_rank among all tissues.
        min_expr (float): minimal expression

    Returns:
        dict of list: A signature dictionary.

        Example::

            {
                "tissue1" : [list, of, gene, ids],
                "tissue2" : [list, of, other, genes],
                ...
            }

    """
    signatures = {tissue: [] for tissue in df_aggr}
    for i in range(df_aggr.shape[0]):
        row = df_aggr.iloc[i, :].sort(inplace=False, ascending=False)
        if row.iloc[0] >= min_expr:
            while gini(row) >= min_gini:
                signatures[row.index[0]].append(i)
                row.drop(row.index[0], inplace=True)
    signatures = {tissue: set(genes) for tissue, genes in signatures.items()}
    return signatures


class CleverGiniSignatureGenerator(SignatureGenerator):
    """
    Use gini index to generate gene signatures.

    Genes, which are specific for a tissue result in a high gini index,
    whereas genes equally present in all tissues have a gini index close to zero.

    The idea is, that genes with a high gini index will reliably identify
    their tissue of origin.

    Args:
        expr (np.ndarray): m x n matrix with m samples and n genes
        target (array-like): m-vector with true tissue for each sample
        min_gini (float): gini cutoff, genes need to have a gini index larger than this value.
        max_rk (int): rank cutoff, include genes if they rank ``<= max_rank`` among all tissues.
        min_expr (float): genes need to have at least an expression ``>= min_expr`` to be included.
        aggregate_fun (function): function used to aggregate samples of the same tissue.
    """

    def __init__(
        self, expr, target, min_gini=0.7, max_rk=3, min_expr=1, aggregate_fun=np.median
    ):
        super(CleverGiniSignatureGenerator, self).__init__(expr, target)
        self.min_gini = min_gini
        self.max_rk = max_rk
        self.min_expr = min_expr
        self.aggregate_fun = aggregate_fun

    def _mk_signatures(self, expr, target):
        df_aggr = collapse_matrix(
            self.expr, self.target, axis=1, aggregate_fun=self.aggregate_fun
        )
        return get_gini_signatures(
            df_aggr, min_gini=self.min_gini, max_rk=self.max_rk, min_expr=self.min_expr
        )
