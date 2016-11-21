import numpy as np
import pandas as pd

def gini(array):
    """
    Calculate the Gini coefficient of a numpy array.

    based on: https://github.com/oliviaguest/gini

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


def aggregate_dataframe(expr_df, group_series, aggregate_fun = np.median):
    """
    Aggregate dataframe by annotation (collapse samples of the same tissue)

    Args:
        expr_df (pd.DataFrame): data frame with expression data
        group_series (pd.Series): series containing the tissue annotation for each index
        aggregate_fun (function): aggregate to apply, e.g. mean, median

    Returns:
        pd.DataFrame: data frame with rows = rows(expr_df) = genes and cols = tissues

    """
    group = group_series.groupby(group_series)
    df_aggr = pd.DataFrame()
    for name, series in group:
        df_aggr[name] = expr_df.iloc[:,series.index].apply(aggregate_fun, axis=1)
    return df_aggr


def jaccard_ind(set1, set2):
    """

    Args:
        set1 (set):
        set2 (set):

    Returns:

    """
    n = len(set1.intersection(set2))
    return np.divide(n, (len(set1) + len(set2) - n))


def get_gini_signatures(df_aggr, min_gini=.7, max_rk=3, min_expr=1):
    """
    Retreive gene signatures from an aggregated data frame

    Args:
        df_aggr (pd.DataFrame): Dataframe rows = genenames, cols = tissue aggregates
        min_gini (float): gini cutoff, genes need to have a gini index larger than this value.
        max_rk (int): rank cutoff
        min_expr (float): minimal expression

    Returns:
        dict of list: tissue -> [list, of, gene, ids]

    """
    df_aggr = df_aggr[df_aggr.apply(np.max, axis=1) >= min_expr]
    expr_gini = df_aggr.apply(gini, axis=1)
    df_aggr = df_aggr[expr_gini >= min_gini]
    df_aggr_rank = df_aggr.rank(axis=1)
    signatures = {}
    for tissue in df_aggr:
        signatures[tissue] = set(df_aggr.loc[df_aggr_rank[tissue] <= max_rk].index)
    return signatures
