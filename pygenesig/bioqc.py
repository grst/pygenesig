import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri
from validation import SignatureTester
numpy2ri.activate()

base = importr("base")
bioqc = importr("BioQC")
biobase = importr("Biobase")

wmw_test = bioqc.wmwTest

read_gmt = bioqc.readGmt



"""
mk_eset

Args:
    exprs:
    gene_symbols:
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
    pass