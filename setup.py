from distutils.core import setup
from setuptools import find_packages

# read the contents of your README file
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.rst"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="pygenesig",
    version="0.2.1",
    packages=find_packages(),
    description="Create and Validate Gene Signatures",
    author="Gregor Sturm",
    author_email="mail@gregor-sturm.de",
    url="https://github.com/grst/pygenesig",  # use the URL to the github repo
    keywords=["bioinformatics", "gene expression", "signatures"],
    license="MIT",
    install_requires=[
        "pyyaml",
        "seaborn",
        "numpy",
        "pandas",
        "rpy2",
        "scikit-learn",
        "dask",
        "distributed",
    ],
    extras_require=dict(
        tests=["nose", "black"],
        docs=["recommonmark", "sphinx", "sphinx_rtd_theme"],
    ),
    classifiers=[],
    include_package_data=True,
    long_description_content_type="text/x-rst",
    long_description=long_description,
)
