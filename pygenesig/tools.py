import numpy as np
import itertools
import pandas as pd
from pprint import pprint
from pygenesig.file_formats import load_gmt


def translate_signatures(signatures, rosetta, ignore_missing=False):
    """
    Translate gene identifiers in a signature dictionary.

    Args:
        signatures (dict of list): signature dictionary
        rosetta (dict): translation table mapping one gene identifier to another
        ignore_missing (boolean): If true, no error will be raised if an identifier is not in
            the translation dictionary. Respective entries will be skipped.

    Returns:
        dict of list: translated signature dictionary

    Raises:
        KeyError: if a gene is not in the rosetta dictionary unless ignore_missing is specified
    """
    if ignore_missing:
        # remove genes from signature which is not in rosetta.
        signatures = {
            tissue: [gene for gene in genes if gene in rosetta]
            for tissue, genes in signatures.items()
        }
    return {
        tissue: [rosetta[gene] for gene in genes]
        for tissue, genes in signatures.items()
    }


def combine_signatures(*args, function=set.intersection):
    """Combine signatures (e.g. by taking the intersection)
    Args:
        *args: list of signature dictonaries
        function: set operation to combine the signatures

    Returns:
        dict of set: combined signature dictionary.

    >>> s1 = {"A": {1, 2, 3}, "B": {2, 3, 4}}
    >>> s2 = {"A": {1, 3, 4, 5}, "B": {42}}
    >>> pprint(combine_signatures(s1, s2))
    {'A': {1, 3}, 'B': set()}
    """
    assert len(args) > 0, "No signatures provided"

    keys = args[0].keys()
    for sig in args:
        assert sig.keys() == keys, "All signature dictonaries must have identical keys."

    return {k: function(*[set(sig[k]) for sig in args]) for k in keys}


def jaccard_ind(set1, set2, *args):
    """
    Computes the Jaccard-Index of two or more sets.

    Args:
        set1 (list-like):
        set2 (list-like):
        *args: arbitrary number of more sets.

    Returns:
        float: jaccard index of all sets

    """
    set1 = set(set1)
    set2 = set(set2)
    i = len(set.intersection(set1, set2, *args))
    u = len(set.union(set1, set2, *args))
    return np.divide(i, u)


def pairwise_jaccard_ind(list_of_signatures):
    """
    Compute the pairwise jaccard index for a list
    of signature sets. Useful for calculating the pariwise overlap
    between the different crossvalidation folds.

    Args:
        list_of_signatures: list of signature dicts.

    Returns:
        dict: signature_name -> [list, of, jaccard, indices]

    Note:
        takes the signature names from the first dict in list_of_signatures
        to build the output dictionary.
    """
    assert len(list_of_signatures) > 0, "no signatures provided."
    pairwise_jacc = {}
    for signame in list_of_signatures[0]:
        pairwise_jacc[signame] = []
        for sigset1, sigset2 in itertools.product(
            list_of_signatures, list_of_signatures
        ):
            pairwise_jacc[signame].append(
                jaccard_ind(sigset1[signame], sigset2[signame])
            )
    return pairwise_jacc


def performance_per_tissue(list_of_confusion_matrices, sig_labels, perf_fun):
    """
    Compute per-tissue performance measures from all-against-all confusion matrices.

    Args:
        list_of_confusion_matrices (list of np.array):  list of confusion matrices
        sig_labels (array-like): list of signatures in the same order as in the confusion matrices.
        perf_fun (function): ``(TP, FN, TP, TN)`` computing a performance measure from the binary confusion matrix.
            See ``perfmeasures`` module.

    Returns:
        dict: signature_name -> list of performance meausures for each confusion matrix provided.

    """
    assert len(list_of_confusion_matrices) > 0, "no matrices provided."
    res = {}
    for i, sig in enumerate(sig_labels):
        res[sig] = []
        for confmat in list_of_confusion_matrices:
            TP = confmat[i, i]
            FN = np.sum(confmat[i, :]) - TP
            FP = np.sum(confmat[:, i]) - TP
            TN = np.sum(confmat) - TP - FN - FP
            res[sig].append(perf_fun(TP, FN, FP, TN))
    return res


def jaccard_mat(sigs1, sigs2, colname1="set_1", colname2="set_2", as_matrix=False):
    """
    Compute a matrix of jaccard indices to compute the overlap of two signature sets.

    Args:
        sigs1: signature dictionary
        sigs2: signature dictionary
        colname1: Name of the column for sigs1 in the dataframe
        colname2: Name of the column for sigs2 in the dataframe
        as_matrix: if False, a long-form dataframe will be returned, if True, a 2d matrix will be returned instead.

    Returns:
        pd.DataFrame: Matrix of Jaccard indices in long format

    Plot the overlap of signatures:
    >>> import seaborn as sns
    >>> signatures = load_gmt("tests/bioqc/test_bioqc_log_pvalue.gmt")
    >>> df = jaccard_mat(signatures, signatures)
    >>> sns.heatmap(df.pivot(*df.columns))  # doctest: +ELLIPSIS
    <AxesSubplot:...>
    """
    jaccard_list = []
    for name1, genes1 in sigs1.items():
        for name2, genes2 in sigs2.items():
            jaccard_list.append((name1, name2, jaccard_ind(set(genes1), set(genes2))))

    df = pd.DataFrame(jaccard_list, columns=(colname1, colname2, "jaccard index"))
    if as_matrix:
        return df.pivot(index=colname1, columns=colname2)
    return df


def collapse_matrix(mat, group_by, axis=0, aggregate_fun=np.median):
    """
    Aggregate expression by annotation (collapse samples of the same tissue)

    Args:
        mat (np.array): m x n gene expression matrix with m genes and n samples.
        group_by (list-like): list of length m (if axis=0) or list of length n (if axis=1)
        axis (int): 0 for rows, 1 for columns
        aggregate_fun (function): aggregate to apply, defaults to ``numpy.median``

    Returns:
        pd.DataFrame: collapsed matrix with annotation from `group_by`.

    """
    if axis == 0:
        assert mat.shape[0] == len(group_by)
    elif axis == 1:
        assert mat.shape[1] == len(group_by)

    mat_df = pd.DataFrame(mat)
    group_by = list(group_by)  # strip index from series
    return mat_df.groupby(group_by, axis=axis).aggregate(aggregate_fun)


def normalize(array):
    """normalize a vector to values between 0 and 1"""
    amax = float(np.nanmax(array))
    amin = float(np.nanmin(array))
    if amax - amin == 0:
        return [0] * len(array)
    array = [(x - amin) / (amax - amin) for x in array]
    return array


def normalize_sum(array):
    """normalize a vector to sum to 1"""
    s = sum(array)
    return [x / s for x in array]
