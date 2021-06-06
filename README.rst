Pygenesig, a framework to generate and validate tissue-specific gene signatures.
===================================================================================

|tests| |docs| |pypi| |black|

.. |tests| image:: https://github.com/grst/pygenesig/actions/workflows/python-package.yml/badge.svg
    :target: https://github.com/grst/pygenesig/actions/workflows/python-package.yml
    :alt: Build Status

.. |docs| image:: https://readthedocs.org/projects/pygenesig/badge/?version=latest
    :target: https://pygenesig.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. |pypi| image:: https://img.shields.io/pypi/v/pygenesig?logo=PyPI
    :target: https://pypi.org/project/pygenesig/
    :alt: PyPI

.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: The uncompromising python formatter


Gene signatures are sets of genes derived from gene expression data, which identify
a certain tissue, cell type, pathway, etc. *Pygenesig* provides a framework to create
and validate such signatures. The package is easily extensible to add new methods
for signature creation and testing.

Getting started
^^^^^^^^^^^^^^^
Please refer to the `documentation <https://pygenesig.readthedocs.io>`_, in particular
the sections

- `Preparing data <https://pygenesig.readthedocs.io/en/latest/prepare_data.html>`_
- `Creating signatures <https://pygenesig.readthedocs.io/en/latest/use_pygenesig.html>`_, and
- `Testing signatures <https://pygenesig.readthedocs.io/en/latest/use_pygenesig.html#testing-signatures>`_.

If you want to implement own methods for creating and testing signatures, please take a
look at the

- `Developer's guide <https://pygenesig.readthedocs.io/en/latest/developers_guide.html>`_ and the
- `API documentation <https://pygenesig.readthedocs.io/en/latest/apidoc.html>`_.


Installation
^^^^^^^^^^^^

You need to have Python 3.7 or newer installed on your system. Some methods for creating
or testing signatures additionally require R. If you don't have
Python installed, we recommend installing `Miniforge <https://github.com/conda-forge/miniforge/releases>`_.

There are several alternative options to install pygenesig:

1) Install pygenesig in a self-contained conda environment:

   This is the most reliable option to make both R and Python work. Make sure you
   have the `conda-forge` and the `bioconda` channels set-up with the correct priorities
   as `described in the Bioconda documentation <https://bioconda.github.io/user/install.html#set-up-channels>`_.

   .. code-block::

      # use `mamba` instead of `conda` for more speed
      mamba create -n pygenesig python=3.8 pip bioconductor-edger bioconductor-bioqc
      conda activate pygenesig
      pip install pygenesig

2) Install pygenesig via pip and R packages manually

   .. code-block::

     pip install pygenesig

   Then, in an R:

   .. code-block:: r

     install.packages("BiocManager")
     BiocManager::install(c("edgeR", "BioQC"))

   Usually, if `R` is in your `PATH`, `rpy2 <https://rpy2.github.io/>`_ automatically
   detects your R installation. If you get an error message while importing `pygensig`,
   try setting the `R_HOME` environment variable before importing pygenesig:

   .. code-block:: python

    import os
    os.environ["R_HOME"] = "/usr/lib/R"
    import pygenesig

3) Install the latest development version from GitHub:

   .. code-block::

     pip install git+https://github.com/grst/pygenesig.git@master

   You'll need to separately install R packages as described above.


Release notes
^^^^^^^^^^^^^
See the `release section <https://github.com/grst/pygenesig/releases>`_.

Contact
^^^^^^^
Please use the `issue tracker <https://github.com/grst/pygenesig/issues>`_.
