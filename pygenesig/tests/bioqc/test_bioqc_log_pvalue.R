################
# Create testcase for BioQC "test_log_pvalue" in python. 
################

library(Biobase)
library(BioQC)
expr = matrix(c(5, 5, 0, 0, 0, 0, 5, 5, 3, 2, 2, 2, 0, 0, 0, 1), ncol=4, byrow = T)
eset = new("ExpressionSet", 
           exprs=expr, 
           featureData=new("AnnotatedDataFrame",
                           data.frame(GeneSymbol=as.character(0:3))))

gmt = readGmt("./test_bioqc_log_pvalue.gmt")

bioqc_res = wmwTest(eset, gmt, valType="p.greater", col="GeneSymbol")
bioqc_res = absLog10p(bioqc_res)
write.csv(bioqc_res, "test_bioqc_log_pvalue.csv")
