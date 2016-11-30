# pygenesig
A python package to create and validate gene signatures. 

Gene signatures are sets of genes derived from gene expression data, which identify a certain tissue, cell type, pathway, etc. This package provides a framework to create and validate such signatures. The package is easily extensible to add new methods for signature creation and testing. 

## Contents:

```eval_rst
.. toctree::
    :maxdepth: 2 

    prepare_data.md
    pygenesig.rst
```

## Getting started:

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

