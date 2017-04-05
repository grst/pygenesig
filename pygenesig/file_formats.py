import pandas as pd
import numpy as np
import itertools


def write_expr(expression_matrix, file):
    """Store a m x n gene expression matrix as numpy object. """
    np.save(file, expression_matrix)


def read_expr(expr_file):
    """Read a m x n gene expression matrix from a numpy object. """
    return np.load(expr_file)


def write_target(target_array, file):
    """Given a m x n gene expression matrix with m genes and n samples. Write an n-array with
    one target annotation for each sample. """
    np.savetxt(file, target_array, delimiter=",", fmt="%s")


def read_target(target_file):
    """Given a m x n gene expression matrix with m genes and n samples. Read an n-array with
    one target annotation for each sample. """
    return np.genfromtxt(target_file, dtype=str, delimiter=",")


def write_rosetta(rosetta_array, rosetta_file):
    """Given a m x n gene expression matrix with m genes and n samples. Write a m-array
    with one identifier for each gene.

    This can be used to map the index-based signature back to gene symbols. """
    np.savetxt(rosetta_file, rosetta_array, delimiter=",", fmt="%s")


def read_rosetta(rosetta_file, inverse=False):
    """Given a m x n gene expression matrix with m genes and n samples. Read an m-array
    with one identifier for each gene.

    This will be converted into a dictionary mapping the gene-index
    to the gene identifier. This can be used to map the index-based signatures
    back to gene symbols.

    Args:
        rosetta_file (str):
        inverse (boolean): If true, map gene-symbol to index.

    Important::
        Be carefule when using `inverse = True` when the list in
        rosetta_file is not unique. In that case only the last
        entry makes it into the list!

    """
    fdata = np.genfromtxt(rosetta_file, dtype=str, delimiter=",")
    if inverse:
        return {
            # '-' is an artifact from ribiosAnnotation
            gene_symbol: i for i, gene_symbol in enumerate(fdata) if gene_symbol != '-'
        }
    else:
        return dict(enumerate(fdata))


def read_gct(file):
    """
    Read a `GCT file`_ to a gene expression matrix.

    Args:
        file: path to GCT file

    Returns:
        np.array: gene expression matrix as 2d numpy array. (rows = genes, cols = samples)

    .. _GCT file:
        http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#gct

    """
    gct = pd.read_csv(file, sep="\t", skiprows=2, index_col=0)
    return gct.iloc[:, 1:].as_matrix()  # get rid of description column


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
            genes = sorted([str(g) for g in signatures[sig]])
            f.write("\t".join(itertools.chain([sig, description], genes)) + "\n")


def load_gmt(file):
    """
    Read a `GMT file`_ into a signature dictionary.

    Args:
        file: path to GMT file

    Returns:
        dict of list: signature dictionary

        Example::

            {
                "tissue1" : [list, of, gene, ids],
                "tissue2" : [list, of, other, genes],
                ...
            }

    .. _GMT file:
        http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29


    """
    signatures = {}
    with open(file) as f:
        for line in f.readlines():
            cols = line.strip().split("\t")
            name = cols[0]
            genes = cols[2:]
            signatures[name] = genes
    return signatures


