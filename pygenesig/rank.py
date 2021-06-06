from pygenesig.validation import SignatureGenerator
from pygenesig.tools import collapse_matrix
import numpy as np
import pandas as pd
from collections import Counter
import itertools


class RankSignatureGenerator(SignatureGenerator):
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
        self,
        expr,
        target,
        include_fraction=0.1,
        not_in_fraction=0.5,
        max_signatures_per_gene=1,
    ):
        super(RankSignatureGenerator, self).__init__(expr, target)
        self.include_fraction = include_fraction
        self.not_in_fraction = not_in_fraction
        self.max_signatures_per_gene = max_signatures_per_gene

    def _mk_signatures(self, expr, target):
        df_aggr = collapse_matrix(expr, target, axis=1)
        include_index = np.array(range(int(self.include_fraction * df_aggr.shape[0])))
        not_in_index = np.array(range(int(self.not_in_fraction * df_aggr.shape[0])))
        tissue_series = {
            tissue: series.sort_values() for tissue, series in df_aggr.iteritems()
        }
        signatures = {
            tissue: series.iloc[include_index].index
            for tissue, series in tissue_series.items()
        }
        gene_count = Counter(
            itertools.chain.from_iterable(
                series.iloc[not_in_index].index for series in tissue_series.values()
            )
        )
        signatures_filtered = {
            tissue: [i for i in genes if gene_count[i] <= self.max_signatures_per_gene]
            for tissue, genes in signatures.items()
        }
        return signatures_filtered
