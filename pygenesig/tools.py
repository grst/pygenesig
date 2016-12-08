import numpy as np
import itertools


def write_gmt(signatures, file, description="na"):
    """
    Writes signatures to a `GMT file`_.

    Args:
        signatures (dict of list): signature dictionary
        file: path to output file
        description: text to fill in the gmt description field.

    .. _GMT file:
        http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29

    """
    with open(file, 'w') as f:
        for sig, genes in sorted(signatures.items()):
            genes = [str(g) for g in signatures[sig]]
            f.write("\t".join(itertools.chain([sig, description], genes)) + "\n")


def jaccard_ind(set1, set2, *args):
    """
    Computes the Jaccard-Index of two or more sets.

    Args:
        set1 (set):
        set2 (set):
        *args: arbitrary number of more sets.

    Returns:
        float: jaccard index of all sets

    """
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
        for sigset1, sigset2 in itertools.product(list_of_signatures, list_of_signatures):
            pairwise_jacc[signame].append(jaccard_ind(sigset1[signame], sigset2[signame]))
    return pairwise_jacc


def performance_per_tissue(list_of_confusion_matrices, sig_labels, perf_fun):
    """
    Compute per-tissue performance measures from all-against-all confusion matrices.

    Args:
        list_of_confusion_matrices (list of np.array):  list of
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



