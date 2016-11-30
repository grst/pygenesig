# Using `pygenesig` to validate signatures. 

## Concept 
The workflow for validating signatures is illustrated in the following figure: 
<flowchart> 

An instance of the `SignatureGenerator` is used to generate signatures. 
An instance of the `SignatureTester` is used to predict the annotation using the signatures and compares them to the standard of truth. 

The classes can easily be combined to perform a cross-validation. 

## Data objects
### Gene expression matrix
A $m \times n$ gene expression matrix with $m$ genes and $n$ samples. 
Depending on the method you may want to use either RPKM or raw counts. 

### Target annotation
A numpy array of length $n$ containing the target classes for each sample in the gene expression matrix. 

The annotation needs to be in the same order as the columns in the matrix. 

### Signature Dictionary
Signatures are stored in a simple dictionary, mapping the signature name to a list of gene indices, for example:

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

<!-- case studies -->
## Perform crossvalidation on a dataset. 
< include output of jupyter notebook > 


## Validate signatures on a different dataset. 
< include output of jupyter notebook >
