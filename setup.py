from distutils.core import setup
from setuptools import find_packages

setup(
    name="pygenesig",
    version="0.2.0",
    packages=find_packages(),
    description="Validate Gene Signatures",
    author="Gregor Sturm",
    author_email="gregor.sturm@cs.tum.edu",
    url="https://github.com/grst/chunksub",  # use the URL to the github repo
    keywords=[
        "bioinformatics",
        "sge",
        "torque",
        "hpc",
        "slurm",
        "parallel",
    ],  # arbitrary keywords
    license="GPLv3",
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
    classifiers=[],
    include_package_data=True,
)
