# Getting started with pygenesig
*Pygenesig* is a python package to create and validate gene signatures. 

Gene signatures are sets of genes derived from gene expression data, which identify a certain tissue, cell type, pathway, etc. This package provides a framework to create and validate such signatures. The package is easily extensible to add new methods for signature creation and testing. 


```eval_rst
.. toctree::
    :maxdepth: 2 
    :caption: User's Guide
    :hidden:

    self
    prepare_data.md
    use_pygenesig.md

.. toctree::
    :maxdepth: 2
    :caption: Developer's Guide
    :hidden:

    developers_guide.md
    apidoc.md
```

## Dependencies
The following python packages are required (install with `pip install`)
```
rpy2
dask
numpy
pandas
scikit-learn
```

**rpy2**

Signature testing is based on [BioQC](https://accio.github.io/BioQC) which is an `R` package. Pygenesig relies on [rpy2](http://rpy2.bitbucket.org/) to use the R package. If you don't have a working R installation on your system, you can get one through [Anaconda](https://www.continuum.io/conda-for-r). 

For pygenesig to work, you need to install the BioQC package in R. You can get the package from [bioconductor](https://bioconductor.org/packages/release/bioc/html/BioQC.html) or the development version from [github](https://github.com/Accio/BioQC):

```r
R> source("http://bioconductor.org/biocLite.R")
R> biocLite("BioQC")
```


**dask**

Crossvalidation can be comutationally intensive. Pygenesig uses [dask](http://dask.readthedocs.io/en/latest/) to support parallel cross-validation. Depending on your setup, you can run code with multi threading on a single PC or on multiple nodes on a high performance cluster. 


## Installing pygenesig

First, clone the repository from github:
```
git clone git@github.com:grst/gene-set-study.git
```

Second, install the package:
```
python3 setup.py develop
```
The `develop` keywork links the python files from the repository to the packages-folder instead of copying them. Like this, any change you make in the source code is immediately available to the python interpreter. 

Now you can use the package, for example from a `jupyter notebook`. 

