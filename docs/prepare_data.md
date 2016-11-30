# Prepare data

## Creating signatures
For creating signatures, we essentially need two pieces of information: The gene expression data and the target annotation (e.g. tissue). Depending on the method, we may need additional covariates, but let's focus on the simple case here: 

### Gene expression matrix
A `$m \times n$` gene expression matrix with `$m$` genes and `$n$` samples. 
Depending on the method you may want to use either RPKM or raw counts. 

The object needs to be a 2d `numpy.array`. You can easily derive that matrix from a gct file using `pandas`

```python
import pandas as pd
gct = pd.read_csv("expression.gct", sep="\t", skiprows=2, index_col=0)
# the first column contains the description, which we do not want to include in the matrix. 
expr = gct.iloc[:, 1:].values
```

You may want to save the matrix for later use:

```python
import numpy as np
np.save("exprssion_matrix.npy", expr)
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
meta = pd.read_csv("meta.tsv", sep="\t", index_col=0)
target = np.array(meta.tissue)
```

```eval_rst
.. Important::
    Make sure that

    - the number of columns in the gene expression matrix equals the length of the target annotation array. You can check that with ```expr.shape[1] == target.shape[0]```
    - the elements in `target` are in the same order as the columns in `exprs`
```



## Testing signatures
For testing signatures, we again need the [target annotation](#target-annotation) as *standard of truth*. Additionally we need to provide the signatures we want to test. 

### Signature Dictionary
Signatures are represented as a simple dictionary, mapping the signature name to a list of gene indices, for example:

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

