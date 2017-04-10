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


def read_rosetta(rosetta_file, as_dict=True, inverse=False):
    """Given a m x n gene expression matrix with m genes and n samples. Read an m-array
    with one identifier for each gene.

    If `as_dict` is True, this will be converted into a dictionary mapping the gene-index
    to the gene identifier. This can be used to map the index-based signatures
    back to gene symbols.

    If `inverse` is True, the mapping is inverse, and the gene-symbol will be mapped to the
    gene-index.

    Args:
        rosetta_file (str): path to file
        as_dict (boolean): If True, a dictionary will be returned. Else it will be a flat numpy.array.
        inverse (boolean): If true, map gene-symbol to index.

    Important::
        Be carefule when using `inverse = True` when the list in
        rosetta_file is not unique. In that case only the last
        entry makes it into the list!

    Returns:
        dict or np.array: mapping or numpy array.

    """
    fdata = np.genfromtxt(rosetta_file, dtype=str, delimiter=",")
    if as_dict:
        return make_dict(fdata, inverse)
    else:
        return fdata


def make_dict(array, inverse=False):
    """
    convert an array to a dictonary mapping the index to the corresponding array entry.
    Use `inverse` to reverse the mapping, i.e. the array entry to the index.

    Args:
        array (array-like):
        inverse (boolean):

    Returns:
        dict: the map

    """
    if inverse:
        return {
            # '-' is an artifact from ribiosAnnotation
            gene_symbol: i for i, gene_symbol in enumerate(array) if gene_symbol != '-'
        }
    else:
        return dict(enumerate(array))


def read_gct(file, include_pfdata=False, rosetta_as_dict=True, inverse=False):
    """
    Read a `GCT file`_ to a gene expression matrix.

    Args:
        file (str): path to GCT file
        include_pfdata (boolean): If true, target will be read from the column names and gene symbols from the description column.
        rosetta_as_dict (boolean): If True, fdata/rosetta will be a dict, otherwise an np.array.
            See `pygenesig.file_formats.read_rosetta`.
        inverse (boolean): reverse the rosetta dict. Not applicable if rosetta_as_dict=False. See
            `pygenesig.file_formats.read_rosetta`.

    Returns:
        np.array or (np.array, np.array, dict or np.array): gene expression matrix as 2d numpy array. (rows = genes, cols = samples)

    .. _GCT file:
        http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#gct

    """
    gct = pd.read_csv(file, sep="\t", skiprows=2, index_col=0)
    exprs = gct.iloc[:, 1:].as_matrix()  # get rid of description column
    if include_pfdata:
        fdata = gct.Description.values
        # pandas suffixes .NNN to duplicate column names which we will remove here.
        split_name = np.vectorize(lambda x: ".".join(x.split(".")[:-1]) if "." in x else x)
        target = split_name(gct.columns[1:].values)
        if rosetta_as_dict:
            fdata = make_dict(fdata, inverse=inverse)
        return exprs, target, fdata
    else:
        return exprs


def write_gct(file, exprs, target=None, description=None, name=None):
    """
    Write a gct file.

    Args:
        file (str): path to output file
        exprs (np.array): m x n matrix with m genes and n samples
        target: n-array with the labels for the n samples
        description: m array with a description for each gene (e.g. gene symbol)
        name: m array with the name for each gene (e.g. gene index)

    """
    if description is None:
        description = np.repeate("na", exprs.shape[0])
    if name is None:
        name = np.arange(0, exprs.shape[0])
    if target is None:
        target = np.arange(0, exprs.shape[1])
    assert exprs.shape[0] == description.size == name.size
    assert exprs.shape[1] == target.size

    gct = pd.DataFrame(exprs)
    gct.columns = target
    fdata = pd.DataFrame({'NAME': name, 'Description': description})
    gct = pd.concat((fdata, gct), axis=1)
    gct.set_index('NAME', inplace=True)

    with open(file, 'w') as f:
        f.write("#1.2\n")
        f.write("{} {}\n".format(*exprs.shape))

    gct.to_csv(file, mode='a', sep="\t")


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
    Deprecated.
    Alias for read_gmt.
    """
    return read_gmt(file)


def read_gmt(file):
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


