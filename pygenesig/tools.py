import numpy as np
import itertools


def which(cond):
    """cond is boolean mask"""
    return [i for i, b in enumerate(cond) if b]


def write_gmt(signatures, file, description="na", order=None):
    """
    Writes signatures to a GMT file.

    Args:
        signatures (dict of iterable): dictionary 'signature name' -> ['list', 'of', 'gene', 'names']
        file: path to output file
        description: text to fill in the gmt description field.

    Raises:
        KeyError: if there is an item in order which is not in signatures.

    Note:
        File format specification: http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29

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

    Note: takes the signature names from the first dict in list_of_signatures
        to build the output dictionary.
    """
    assert len(list_of_signatures) > 0, "no signatures provided."
    pairwise_jacc = {}
    for signame in list_of_signatures[0]:
        pairwise_jacc[signame] = []
        for sigset1, sigset2 in itertools.product(list_of_signatures, list_of_signatures):
            pairwise_jacc[signame].append(jaccard_ind(sigset1[signame], sigset2[signame]))
    return pairwise_jacc

