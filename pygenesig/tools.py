import numpy as np
import itertools


def which(cond):
    """cond is boolean mask"""
    return [i for i, b in enumerate(cond) if b]


def write_gmt(signatures, path, description="na"):
    """
    Writes signatures to a GMT file.

    Args:
        signatures (dict of iterable): dictionary 'signature name' -> ['list', 'of', 'gene', 'names']
        path: output file

    Note:
        File format specification: http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29

    """
    with open(path, 'w') as f:
        for sig, genes in sorted(signatures.items()):
            f.write("\t".join(itertools.chain([sig, description], genes)) + "\n")


def jaccard_ind(set1, set2):
    """
    Computes the Jaccard-Index of two sets.

    Args:
        set1 (set):
        set2 (set):

    Returns:
        float: jaccard index of set1 and set2

    """
    n = len(set1.intersection(set2))
    return np.divide(n, (len(set1) + len(set2) - n))
