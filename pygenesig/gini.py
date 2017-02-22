"""
The Gini-coefficient (gini index) measures the inequality among values of
a frequency distribution. Originally applied to measure income inequality
we can also use it on gene expression data to find genes that are
over-represented in certain samples.
"""

import numpy as np
import pandas as pd
from pygenesig.validation import SignatureGenerator


def gini(array):
    """
    Calculate the Gini coefficient of a numpy array.

    Based on: https://github.com/oliviaguest/gini

    Args:
        array (array-like): input array

    Returns:
        float: gini-index of ``array``

    >>> a = np.zeros((10000))
    >>> a[0] = 1.0
    >>> '%.3f' % gini(a)
    '1.000'
    >>> a = np.ones(100)
    >>> '%.3f' % gini(a)
    '0.000'
    >>> a = np.random.uniform(-1,0,1000000)
    >>> '%.2f' % gini(a)
    '0.33'
    """
    # based on bottom eq: http://www.statsdirect.com/help/content/image/stat0206_wmf.gif
    # from: http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    if np.amin(array) < 0:
        array -= np.amin(array)             # values cannot be negative
    array += 1e-12                           # values cannot be 0
    array = np.sort(array)                  # values must be sorted
    index = np.arange(1, array.shape[0]+1)  # index per array element
    n = array.shape[0]                      # number of array elements
    return (np.sum((2 * index - n - 1) * array)) / (n * np.sum(array))  # Gini coefficient


def aggregate_expression(expr, target, aggregate_fun=np.median):
    """
    Aggregate expression by annotation (collapse samples of the same tissue)

    Args:
        expr (np.array): m x n gene expression matrix with m genes and n samples.
        target (list-like): list of length n, containing the target annotation for each sample (e.g. tissue)
        aggregate_fun (function): aggregate to apply, defaults to ``numpy.median``

    Returns:
        pd.DataFrame: m x n pandas dataframe with m = #Genes and n = #target labels

    """
    group_series = pd.Series(target)
    group = group_series.groupby(group_series)
    df_aggr = pd.DataFrame()
    for name, series in group:
        df_aggr[name] = np.apply_along_axis(aggregate_fun, 1, expr[:, series.index])
    return df_aggr


def _apply_gini_to_expression(df_aggr, min_gini, min_expr):
    """
    Helper function to filter dataframe by gini index and expression.

    Args:
        df_aggr (pd.DataFrame): m x n pandas dataframe with m Genes and n tissues
        min_gini: gini cutoff
        min_expr: expression cutoff

    Returns:
        df_aggr: input dataframe filtered by gini and expression
        df_aggr_rank: dataframe of the same shape as df_aggr, containing the ranks of the
            tissues for each gene. Ranks are computed using the 'first' method to obtain
            consecutive ranks for the 'gini-format'
        expr_gini: m-series containing the gini-index for each gene

    """
    df_aggr = df_aggr[df_aggr.apply(np.max, axis=1) >= min_expr]
    expr_gini = df_aggr.apply(gini, axis=1)
    df_aggr = df_aggr[expr_gini >= min_gini]
    # min, because consecutive rank needed for gini method.
    df_aggr_rank = df_aggr.rank(axis=1, ascending=False, method='first')
    expr_gini = expr_gini[df_aggr.index]
    return df_aggr, df_aggr_rank, expr_gini


def get_rogini_format(df_aggr, min_gini=.7, max_rk=3, min_expr=1):
    """
    Imitate the *rogini* output format.

    Args:
        df_aggr (pd.DataFrame): m x n pandas DataFrame with m Genes and n tissues
        min_gini (float): gini cutoff, genes need to have a gini index larger than this value.
        max_rk (int): rank cutoff, include genes if they rank <= max_rank among all tissues.
        min_expr (float): minimal expression

    Returns:
        pd.DataFrame: equal to *rogini* output.

        Example::

            GENEID  CATEGORY        VALUE   RANKING GINI_IDX
            54165   Whole blood (ribopure)  491.359000      1       0.441296
            54165   CD34 cells differentiated to erythrocyte lineage        148.124000      2       0.441296
            54165   Mast cell - stimulated  68.973000       3       0.441296
            101927596       CD4+CD25-CD45RA+ naive conventional T cells     15.505000       1       0.948804
            101927596       CD8+ T Cells (pluriselect)      15.376000       2       0.948804
            101927596       CD4+CD25-CD45RA- memory conventional T cells    10.769000       3       0.948804

    """
    df_aggr, df_aggr_rank, expr_gini = _apply_gini_to_expression(df_aggr, min_gini, min_expr)
    gini_rows = []
    for i in range(df_aggr.shape[0]):
        geneid = df_aggr.index[i]
        row = df_aggr.iloc[i]
        rank_row = df_aggr_rank.iloc[i]
        for rk in range(1, max_rk+1):
            gini_rows.append([
                geneid,
                df_aggr.columns[rank_row == rk][0],
                row[rank_row == rk][0],
                rk,
                expr_gini[geneid]
            ])
    columns = ["GENEID", "CATEGORY", "VALUE", "RANKING", "GINI_IDX"]
    df = pd.DataFrame(gini_rows)
    df.columns = columns
    return df


def get_gini_signatures(df_aggr, min_gini=.7, max_rk=3, min_expr=1):
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
    df_aggr, df_aggr_rank, expr_gini = _apply_gini_to_expression(df_aggr, min_gini, min_expr)
    signatures = {}
    for tissue in df_aggr:
        signatures[tissue] = set(df_aggr.loc[df_aggr_rank[tissue] <= max_rk].index)
    return signatures


class GiniSignatureGenerator(SignatureGenerator):
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
    def __init__(self, expr, target, min_gini=.7, max_rk=3, min_expr=1, aggregate_fun=np.median):
        super(GiniSignatureGenerator, self).__init__(expr, target)
        self.min_gini = min_gini
        self.max_rk = max_rk
        self.min_expr = min_expr
        self.aggregate_fun = aggregate_fun

    def _mk_signatures(self, expr, target):
        df_aggr = aggregate_expression(expr, target, aggregate_fun=self.aggregate_fun)
        return get_gini_signatures(df_aggr, min_gini=self.min_gini, max_rk=self.max_rk, min_expr=self.min_expr)

    def get_rogini_format(self):
        df_aggr = aggregate_expression(self.expr, self.target, aggregate_fun=self.aggregate_fun)
        return get_rogini_format(df_aggr, min_gini=self.min_gini, max_rk=self.max_rk, min_expr=self.min_expr)