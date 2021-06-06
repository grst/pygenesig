from distutils.core import setup
from setuptools import find_packages

setup(
    name="pygenesig",
    version="0.2.0",
    packages=find_packages(),
    description="Create and Validate Gene Signatures",
    author="Gregor Sturm",
    author_email="gregor.sturm@cs.tum.edu",
    url="https://github.com/grst/chunksub",  # use the URL to the github repo
    keywords=["bioinformatics", "gene expression", "signatures"],
    license="MIT",
    install_requires=[
        "jinja2",
        "pyyaml",
        "docopt",
        "seaborn",
        "numpy",
        "pandas",
        "rpy2",
        "scikit-learn",
        "dask",
        "distributed",
    ],
    extras_require=dict(
        tests=["nose"], docs=["recommonmark", "sphinx", "sphinx_rtd_theme"]
    ),
    classifiers=[],
    include_package_data=True,
)
