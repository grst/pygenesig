import numpy as np
import pandas as pd
import dask
import dask.threaded
import dask.dataframe as dd


def jaccard_ind(set1, set2):
    """

    Args:
        set1 (set):
        set2 (set):

    Returns:

    """
    n = len(set1.intersection(set2))
    return np.divide(n, (len(set1) + len(set2) - n))


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


def aggregate_dataframe(expr_df, group_series, aggregate_fun=np.median, get=dask.threaded.get):
    """
    Aggregate dataframe by annotation (collapse samples of the same tissue)

    Args:
        expr_df (dask.dataframe): data frame with expression data
        group_series (pd.Series): series containing the tissue annotation for each index
        aggregate_fun (function): aggregate to apply, e.g. mean, median

    Returns:
        pd.DataFrame: data frame with rows = rows(expr_df) = genes and cols = tissues

    """
    group = group_series.groupby(group_series)
    df_aggr = pd.DataFrame()
    for name, series in group:
        mask = [True if i in list(series.index) else False for i in range(len(expr_df.columns))]
        df_aggr[name] = expr_df.loc[:, mask].apply(aggregate_fun, axis=1).compute(get=get)
    return df_aggr


def get_gini_signatures(df_aggr, min_gini=.7, max_rk=3, min_expr=1, get=dask.threaded.get):
    """
    Retreive gene signatures from an aggregated data frame

    Args:
        df_aggr (dask.dataframe): Dataframe rows = genenames, cols = tissue aggregates
        min_gini (float): gini cutoff, genes need to have a gini index larger than this value.
        max_rk (int): rank cutoff
        min_expr (float): minimal expression

    Returns:
        dict of list: tissue -> [list, of, gene, ids]

    """
    df_aggr = df_aggr[df_aggr.apply(np.max, axis=1).compute(get=get) >= min_expr]
    expr_gini = df_aggr.apply(gini, axis=1).compute(get=get)
    df_aggr = df_aggr[expr_gini >= min_gini]
    df_aggr_rank = df_aggr.rank(axis=1).compute(get=get)
    signatures = {}
    for tissue in df_aggr:
        signatures[tissue] = set(df_aggr.loc[df_aggr_rank[tissue] <= max_rk].index)
    return signatures


def mk_gini_signature(gct, col_vars, subset=None, tissue_col="tissue",
                      min_gini=.7, max_rk=3, min_expr=1, get=dask.threaded.get):
    """

    Args:
        gct (pd.DataFrame): gct expression matrix (col = sample)
        col_vars (pd.DataFrame): annotation of the gct matrix (row = sample).
            the samples in col_vars must be a subset of those in gct.
        subset: list of indices of col_vars to be used to generate the signatures.
        tissue_col: name of the column in col_vars that contains the tissue annotation.

    Returns:
        dict of list

    """
    assert all(x in list(gct.columns) for x in list(col_vars.index)), "The samples in col_vars are not a subset of gct"
    if subset is not None:
        col_vars = col_vars.iloc[subset, :]
    gct = gct[list(col_vars.index)]
    df_aggr = aggregate_dataframe(gct, col_vars[tissue_col], get=get)
    df_aggr = dd.from_pandas(df_aggr, npartitions=20)
    return get_gini_signatures(df_aggr, min_gini=min_gini, max_rk=max_rk, min_expr=min_expr, get=get)


def mk_de_signatures(gct, col_vars, covariates=()):
    """
    see mk_gini_signatures.
    Maybe it's time to try rpy2?
    Args:
        gct:
        col_vars:

    Returns:

    """
    return None

