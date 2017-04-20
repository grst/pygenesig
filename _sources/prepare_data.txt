# Prepare data

## For creating signatures
For creating signatures, we essentially need two pieces of information: 
* The *gene expression* data and 
* the target annotation (e.g. tissue). 

Depending on the method, we may need additional covariates, but let's focus on the simple case here: 

### Gene expression matrix
A `$m \times n$` gene expression matrix with `$m$` genes and `$n$` samples. 
Depending on the method you may want to use either RPKM or raw counts. The `expr` object
needs to be a 2d `numpy.array`. 

*pygenesig* provides functions to read the standardized [gct](http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#_Creating_Input_Files_GCT)
format into a gene expression matrix:

```python
from pygenesig.file_formats import read_gct
expr = read_gct('expression.gct')
```

You can also come up with your own methods for reading the expression data, which can be useful if you use different 
formats. Here, we achive the same using `pandas`:

```python
import pandas as pd
gct = pd.read_csv("expression.gct", sep="\t", skiprows=2, index_col=0)
# the first column contains the description, which we do not want to include in the matrix. 
expr = gct.iloc[:, 1:].values
```

You may want to save the matrix for later use:

```python
from pygenesig.file_formats import write_expr, read_expr
write_expr(expr, "exprssion_matrix.npy")

# later, read it with
expr = read_expr("expression_matrix.npy")
```

### Target annotation
A `numpy.array` of length `$n$` containing the target classes for each sample in the gene expression matrix. 
The annotation needs to be in the same order as the columns in the matrix. 

Say, you have a tab-separated file containin meta information for each sample:
```
id  tissue  gender 
1   adipose male
2   adipose female
3   colon   female
...
```

You can again use `pandas` to extract the target annotation:

```python
import pandas as pd
from pygenesig.file_formats import write_target
pheno_data = pd.read_csv("meta.tsv", sep="\t", index_col=0)
target = pheno_data.tissue.values
write_target(target, "target.txt")
```

```eval_rst
.. Important::
    Make sure that

    - the number of columns in the gene expression matrix equals the length of the target annotation array. You can check that with ``expr.shape[1] == target.shape[0]``
    - the elements in ``target`` are in the same order as the columns in ``exprs``
```



## For testing signatures
For testing signatures, we again need the [target annotation](#target-annotation) as *ground truth*. 
Additionally we need to provide the signatures we want to test. 

### Signature Dictionary
Signatures are represented as a simple dictionary, mapping the signature name to a list of genes. A gene is represented by it's row-index in the gene expression matrix:

```python
signatures = {
    'adipose tissue' : [
        42,
        124,
        256,
        1038,
        ... ],
    'skeletal muscle' : [
        52,
        679,
        ... ],
    ...
}
```

If you use *pygenesig* for generating the signatures, you will get such a dictionary by default. 

However, the indices are not meaningful outside the context of *pygenesig*, therefore you probably want 
to change them to something meaningful (e.g. gene symbols) instead before you export them as 
[gmt](http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GMT)  file
for further use in other software. 

This is accomplished using a "rosetta file", containing one gene identifier per line in the same order 
as the rows in the gene expression matrix: 

*rosetta.txt:* 
```
HAS2
TSTA3
ADH1
...
```

```python
from pygenesig.file_formats import read_rosetta, write_gmt
from pygenesig.tools import translate_signatures
rosetta = read_rosetta('rosetta.txt')
signatures_symbol = translate_signatures(signatures, rosetta)

# export as gmt
write_gmt(signatures_symbol, 'my_signatures.gmt')
```

### Reading `.gmt` files. 
This also works in the other direction: 
You can read gene signatures in [gmt](http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GMT) 
format into a signature dictionary. Now the signature is gene-symbol based and we need to translate
it into the index-based representation for use within *pygenesig*

```python
from pygenesig.file_formats import read_gmt
signatures = read_gmt("my_signatures.gmt")
signatures
{
    'adipose tissue' : [
        'HAS2',
        'TSTA3',
        ... ],
    'skeletal muscle' : [
        'CXXC6',
        'TNKS2',
        ... ],
    ...
}
```

```python
from pygenesig.file_formats import read_rosetta
from pygenesig.tools import translate_signatures

rosetta_inv = read_rosetta('rosetta.txt', inverse=True)
signatures_ind = translate_signatures(signatures, rosetta_inv)
```