library(ribiosUtils)
library(ribiosIO)
library(ribiosExpression)
gtex <- read_gct_matrix("/DATA/bi10/biosignature_data/GTex/RocheAnnotated_GTEx_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm__Pilot_2013_01_31.gct")
fd <- read.table("/DATA/bi10/biosignature_data/GTex/RocheAnnotated_GTEx_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm__Pilot_2013_01_31-featureAnnotation.txt",sep="\t", comment.char="", quote="", header=TRUE)
pd <- read.table("/DATA/bi10/biosignature_data/GTex/RocheAnnotated_GTEx_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm__Pilot_2013_01_31-sampleAnnotation.txt", sep="\t", comment.char="", quote="", header=TRUE)
rownames(gtex) <- rownames(fd); colnames(gtex) <- rownames(pd)
eset <- new("ExpressionSet",
            exprs=gtex,
            featureData=new("AnnotatedDataFrame", fd),
            phenoData=new("AnnotatedDataFrame", pd))

visGene <- function(geneSymbol) {
  ind <- which(fData(eset)$GeneSymbol %in% geneSymbol)
  stopifnot(length(ind)==1)
  indExprs <- exprs(eset)[ind,]
  pheno <- pData(eset)
  indFac <- names(sort(tapply(indExprs, pheno$SMTSD, median), decreasing=TRUE))
  cols <- rainbow(length(indFac))
  par(mar=c(18,4,1,1))
  boxplot(log10(indExprs)~factor(pheno$SMTSD, levels=indFac), las=3L,
          ylab="log10(RPKM)", col=cols, main=geneSymbol)
  abline(h=log10(c(1,15,50)), lty=2, lwd=c(1,1.5,2))
  boxplot(log10(indExprs)~factor(pheno$SMTSD, levels=indFac), las=3L,
          ylab="log10(RPKM)", col=cols, add=TRUE)
}
