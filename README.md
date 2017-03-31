# *pygenesig*, a framework to generate and validate tissue-specific gene signatures. 

Gene signatures are sets of genes derived from gene expression data, which identify a certain tissue, cell type, pathway, etc. This package provides a framework to create and validate such signatures. The package is easily extensible to add new methods for signature creation and testing.

Documentaion is available [here](http://grst.github.io/gene-set-study). 

### Signature Checker Pipeline
A wrapper around the *pygenesig* library for
 * generating signatures
 * cross-validatin signatures
 * applying signatures to real-world data
 * exploring the results in a browser

can be found in a [dedicated github repository](https://github.com/grst/pygenesig-pipeline). 

### Case studies
can be found in [pygenesig-example](https://github.com/grst/pygenesig-example)


## Installation
```
git clone git@github.com:grst/pygenesig.git
cd pygenesig
pip install --user -e pygenesig
```
