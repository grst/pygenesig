"""
`BioQC`_ is a R/Bioconductor package to detect tissue heterogeneity in gene expression data. It ships with a set
of 150 tissue signatures and implements a computationally efficient Wilcoxon-Mann-Whitney (WMW) test.

The WMW-test can be applied to check for an over-representation of an arbitrary gene set
in a gene expression sample.

This module implements a python bridge to the R package using ``rpy2`` and implements a
``SignatureTester`` based on BioQC.

.. _BioQC:
    https://accio.github.io/BioQC
"""


import readline  # not used, but needed as a workaround for an import problem.
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
from pygenesig.validation import SignatureTester
from pygenesig.tools import write_gmt
import sklearn.metrics
import numpy as np
import tempfile
numpy2ri.activate()

base = importr("base")
bioqc = importr("BioQC")
biobase = importr("Biobase")


wmw_test = bioqc.wmwTest
"""
Alias to R function ``BioQC::wmwTest`` using rpy2.

Check the `BioQC Documentation`_
for more details.

.. _BioQC Documentation:
    https://bioconductor.org/packages/release/bioc/manuals/BioQC/man/BioQC.pdf
"""

read_gmt = bioqc.readGmt
"""
Alias to R function ``BioQC::readGmt`` using rpy2.

Check the `BioQC Documentation`_
for more details.

.. _BioQC Documentation:
    https://bioconductor.org/packages/release/bioc/manuals/BioQC/man/BioQC.pdf
"""


ro.r('''
    mk_eset = function(exprs, gene_symbols) {
        return(new("ExpressionSet",
            exprs=exprs,
            featureData=new("AnnotatedDataFrame",
                            data.frame(GeneSymbol=gene_symbols))))
     }
     ''')
mk_eset = ro.r['mk_eset']
"""
Generate an R ExpressionSet with ``exprs()`` and ``fData()`` annotation.
fData will contain one single column named GeneSymbol.

Args:
    exprs: ``numpy.array`` m x n expression matrix. ``np.matrix`` does NOT work!
    gene_symbols: ``numpy.array`` m vector containing gene symbols.
"""


class BioQCSignatureTester(SignatureTester):
    """
    Use BioQC to check the validity of gene signatures.

    We use the WMW-test implemented by BioQC to generate a p-value for
    each signature in each sample. If a signature is enriched in a specific sample
    the WMW-test will result in a small p-value.

    We define a sample as `correctly classified` by the signature
    if the signature corresponding to the `true tissue of origin`
    has the lowest p-value of all tested signatures.

    Args:
        expr (np.ndarray): m x n matrix with m samples and n genes
        target (array-like): m-vector with true tissue for each sample
    """

    @staticmethod
    def _reorder(template, target):
        """
        Reorder an array according to a template-array.

        Args:
            template: List of names how the order should be
            target: List of names how the order is

        Returns:
            np.ndarray: Ordered indices to bring target in the order or template.

        .. Note::
            ``template == target[_reorder(template, target)]`` is always ``True``.

        """
        target_array = np.array(target)
        return np.array([np.where(target_array == name)[0][0] for name in template])

    @staticmethod
    def signatures2gmt(signatures):
        """
        Convert signature dictionary into an R GMTList object.

        Args:
            signatures (dict): signature dictionary

        Returns:
            GMTList: R GMTList object containing the signatures.
        """
        with tempfile.NamedTemporaryFile() as gmt_file:
            write_gmt(signatures, gmt_file.name)
            gmt = read_gmt(gmt_file.name)
            return gmt

    @staticmethod
    def run_bioqc(expr, gene_symbols, gmt):
        """
        Apply BioQC to a gene expression matrix.

        Args:
            expr: m x n gene expression matrix with m genes and n samples.
            gene_symbols: m vector containing gene symbols
            gmt: GMTList object containing k signatures. Generate with ``signatures2gmt()``.

        Returns:
            np.array: k x n p-value matrix.
        """
        # although BioQC does support testing on a matrix and
        # list of signatures, I decided to go with the
        # gmt/eset solution. This is easier to debug and it's less easy to confuse things.
        eset = mk_eset(np.array(expr), np.array(gene_symbols))
        bioqc_res = wmw_test(eset, gmt, valType="p.greater")
        return np.array(bioqc_res)

    def _score_signatures(self, expr, signatures):
        gene_symbols = [str(x) for x in range(expr.shape[0])]
        signatures_not_empty = {
            name: genes for name, genes in signatures.items() if len(genes) > 0
        }
        gmt = BioQCSignatureTester.signatures2gmt(signatures_not_empty)
        bioqc_res = self.run_bioqc(expr, gene_symbols, gmt)
        if len(signatures_not_empty) == 1:
            # bug in R which returns a list instead of a matrix for one signature only.
            # make it a matrix again.
            bioqc_res = np.array([bioqc_res])
        bioqc_res_log = -np.log10(bioqc_res)

        reorder = self._reorder(list(base.names(gmt)), self.sort_signatures(signatures))
        # init empty result matrix
        result = np.empty((len(signatures), expr.shape[1]))
        result[:] = np.NAN
        # reassign row by row to result matrix to include empty signatures.
        for i_bioqc, i_result in enumerate(reorder):
            result[i_result, :] = bioqc_res_log[i_bioqc, :]

        return result




