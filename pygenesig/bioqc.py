import readline
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

"""
See R vignette.
"""
wmw_test = bioqc.wmwTest

"""
See R vignette.
"""
read_gmt = bioqc.readGmt


"""
    Generate an R Expression set with exprs() and fData() annotation.
    fData will contain one single column named GeneSymbol.

    Args:
        exprs: numpy.array m x n expression matrix. np.matrix does NOT work!
        gene_symbols: numpy.array m vector containing gene symbols.
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


class BioQCSignatureTester(SignatureTester):
    @staticmethod
    def signatures2gmt(signatures):
        """convert signature dictionary into an R GMTList object. """
        with tempfile.NamedTemporaryFile() as gmt_file:
            write_gmt(signatures, gmt_file.name)
            gmt = read_gmt(gmt_file.name)
            return gmt

    @staticmethod
    def run_bioqc(expr, gene_symbols, gmt):
        """

        Args:
            expr: m x n gene expression matrix with m genes and n samples.
            gene_symbols: m vector containing gene symbols
            gmt: GMTList object containing k signatures. Generate with signatures2gmt().

        Returns:
            np.array: k x n p-value matrix.

        """
        # although BioQC does support testing on a matrix and
        # list of signatures, I decided to go with the
        # gmt/eset solution. This is easier to debug and it's less easy to confuse things.
        eset = mk_eset(np.array(expr), np.array(gene_symbols))
        bioqc_res = wmw_test(eset, gmt, valType="p.greater")
        return np.array(bioqc_res)

    def test_signatures(self, signatures, subset):
        """
        Test Signatures with BioQC. A sample is considered to be
        'true-positive' if the highest scoring signature is the one
        generated from the actual tissue.

        Args:
            signatures:
            subset:

        Returns:
            np.matrix: confusion matrix (generated with sklearn.metrics.confusion_matrix())

        """
        # gene symbols = indices
        gene_symbols = [str(x) for x in range(self.expr.shape[0])]
        gmt = BioQCSignatureTester.signatures2gmt(signatures)
        bioqc_res = self.run_bioqc(self.expr[:, subset], gene_symbols, gmt)

        bioqc_res_log = -np.log10(bioqc_res)
        gmt_signature_names = list(base.names(gmt))

        # predicted signature names for each sample
        predicted = [gmt_signature_names[i] for i in np.argmax(bioqc_res_log, axis=0)]
        actual = self.target[subset]
        return sklearn.metrics.confusion_matrix(actual, predicted)



