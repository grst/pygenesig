library(ribiosUtils)
library(ribiosExpression)
library(BioQC)
library(ribiosAnnotation)
library(ribiosPlot)
library(limma)
library(ribiosIO)
library(edgeR)
library(data.table)
library(readr)
source("funcs.R")

genesetFile <- "data/RocheAnnotated_GTEx_Analysis_RNA-seq_RNA-SeQCv1.1.8_gene_reads__Pilot_2013_01_31.RData"

if(!loadFile(genesetFile)) {
  graw <- read_gct_matrix("data/GTEx_Analysis_RNA-seq_RNA-SeQCv1.1.8_gene_reads__Pilot_V3_patch1.gct")
  ## convert count
  pt.raw <- read.table("data/RocheAnnotated_GTEx_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm__Pilot_2013_01_31-sampleAnnotation.txt",
                   sep="\t", header=TRUE, comment.char="", quote="")
  tissues <- read.table("udis/TissueCVTranslate.txt",
                        sep="\t", header=FALSE, comment.char="", quote="",
                        col.names=c("Group", "SMTSD", "UDISCV"))
  pt <- merge(pt.raw, tissues[,c("SMTSD", "UDISCV")], by="SMTSD")
  pt <- matchColumn(colnames(graw), pt, "SAMPID")
  
  ## annotate feat data by annotateEnsemble
  grawID <- gsub("\\.[0-9]*$", "", rownames(graw))
  ftFile <- "data/RocheAnnotated_GTEx_Analysis_RNA-seq_RNA-SeQCv1.1.8_gene_reads__Pilot_2013_01_31.txt"
  ft <- annotateEnsemble(grawID)
  write.table(ft,ftFile,
              sep="\t", row.names=FALSE, quote=FALSE)

  rownames(graw) <- rownames(ft)
  
  ## annotate pheno data
  stopifnot(identical(colnames(graw),as.character(pt$SAMPID)))
  rownames(pt) <- colnames(graw)
  
  eset <- new("ExpressionSet",
              exprs=graw,
              featureData=new("AnnotatedDataFrame", ft),
              phenoData=new("AnnotatedDataFrame", pt))

  gset <- eset[!is.na(fData(eset)$GeneID),]
  gset_dt = data.table(exprs(gset), keep.rownames = FALSE)
  gset_dt[,GeneID:=as.character(fData(gset)$GeneID)]
  gset_dt = gset_dt[,lapply(.SD, function(x) sum(x, na.rm=TRUE)), by=GeneID]
  
  ft = annotateGeneIDs(gset_dt$GeneID)
  
  geneset = new("ExpressionSet",
                exprs=as.matrix(as.data.frame(gset_dt)[,colnames(exprs(gset))]),
                featureData=new("AnnotatedDataFrame", ft),
                phenoData=new("AnnotatedDataFrame", pData(gset)))

  save(geneset, file=genesetFile)
  write_gct(exprs(geneset), file="data/roche_annotated_read_counts.gct")
  write_tsv(fData(geneset), "data/roche_annotated_fdata.tsv")
  write_tsv(pData(geneset), "data/roche_annotated_pdata.tsv")
}


##----------------------------------------##
## Task 1: average expression by tissue
##----------------------------------------##

fitFile <- "data/RocheAnnotated_GTEx_Analysis_RNA-seq_RNA-SeQCv1.1.8_gene_reads__Pilot_2013_01_31-limma.RData"

if(!loadFile(fitFile)) {
  ageComb <- geneset$AGE
  ageComb[ageComb=="70-79 years"] <- "60-69 years"
  ageComb <- droplevels(ageComb)
  levels(ageComb)[5] <- "Age 60-79"
  mm <- model.matrix(~0+geneset$UDISCV+factor(geneset$GENDER)+geneset$SMRIN+droplevels(ageComb))
  colnames(mm) <- c(levels(geneset$UDISCV),
                    "Gender", "RIN", "Age30-39", "Age40-49", "Age50-59", "Age60-79")
  ## system.time(y <- estimateGLMCommonDisp(y, mm, verbose=TRUE))
  
  ## limma
  nf <- calcNormFactors(geneset)
  expy <- voom(geneset, mm, plot=FALSE, lib.size=colSums(exprs(geneset))*nf)
  fit <- lmFit(expy, mm)
  efit <- eBayes(fit)
  head(test <- topTable(efit, coef="Gender", n=20000))
  ufInd <- seq(from=nlevels(geneset$UDISCV)+1, to=ncol(mm))
  be <- efit$coefficients[,ufInd] %*% t(mm[,ufInd])
  
  log2cpmExp <- expy$E-be
  cpmExp <- 2^log2cpmExp

  log2cpmAvg <- summarizeColumns(log2cpmExp, geneset$UDISCV, median)
  cpmAvg <- 2^log2cpmAvg

  write_gct(cpmExp, file="data/roche_annotated_cpm.gct")
  
  save(mm, efit, log2cpmExp, cpmExp, log2cpmAvg, cpmAvg,
       file=fitFile)
}

cpmFile <- "data/RocheAnnotated_GTEx_Analysis_RNA-seq_RNA-SeQCv1.1.8_gene_reads__Pilot_2013_01_31-cpm.RData"
cpmGct <- "data/RocheAnnotated_GTEx_Analysis_RNA-seq_RNA-SeQCv1.1.8_gene_CPM__Pilot_2013_01_31.gct"

if(!loadFile(cpmFile)) {
  cpmEset <- new("ExpressionSet",
                 exprs=cpmAvg,
                 featureData=featureData(geneset),
                 phenoData=new("AnnotatedDataFrame",
                   data.frame(Tissue=colnames(cpmAvg),
                              row.names=colnames(cpmAvg))))

  options(scipen=3, digits=4)
  ## export GCT
  writeGct(cpmEset, file=cpmGct,
           feat.name="GeneID", feat.desc="GeneSymbol")
  
  save(cpmEset, file=cpmFile)
}
  
##----------------------------------------##
## Task 3: Gini index
##----------------------------------------##
cpmExpGini <- apply(cpmExp, 1L, gini)
cpmAvgGini <- apply(cpmAvg, 1L, gini)
isExp <- rowSums(cpm(exprs(geneset))>=1)>=10L
hist(cpmAvgGini)
 ## Gini index does not differ much when calculated with summarized statistics or individuals
smoothScatter(cpmExpGini, cpmAvgGini,
              xlab="Gini index using individual expression values",
              ylab="Gini index using median expression of tissues")
smoothScatter(cpmExpGini[isExp], cpmAvgGini[isExp],
              xlab="Gini index using individual expression values",
              ylab="Gini index using median expression of tissues")

ipdf("Gini-scatter-indi-avg.pdf") 

cpmAvgRank <- apply(cpmAvg, 1L, function(x) rank(x))
topRnk <- 3
thrRnk <- ncol(cpmAvg)-topRnk
thrGini <- 0.8
cpmIsSig <- t(cpmAvgRank)>thrRnk
cpmThr <- 5

## Gini index larger than the threshold, is generally expressed, and is expressed in the specific tissue (>=1 CPM)
cpmIsSig[cpmAvgGini<thrGini | !isExp,] <- FALSE
cpmIsSig[cpmAvg<cpmThr] <- FALSE


cpmSigList <- apply(cpmIsSig, 2L, which)
cpmSigGeneID <- sapply(cpmSigList, function(x) fData(geneset)$GeneID[x])
cpmSigGeneSymbol <- sapply(cpmSigList, function(x) sort(setdiff(unique(fData(geneset)$GeneSymbol[x]),"-")))
mybar <- function(...) {
  panel.grid(v=-1, h=0)
  panel.barchart(...)
}
summary(cpmSigLen <- sapply(cpmSigList, length))
# barchart(cpmSigLen, panel=mybar,
#         xlab="Signature gene counts")
ipdf("GTEx-Gini_0.8_top3_countBarchart.pdf")

cpmSigGeneSymbol[["Pancreas"]]
cpmSigGeneSymbol[["Adipose - Subcutaneous"]]
cpmSigGeneSymbol[["Adipose - Visceral (Omentum)"]]
cpmSigGeneSymbol[["Artery - Coronary"]]
cpmSigGeneSymbol[["Artery - Aorta"]]
cpmSigGeneSymbol[["Fibroblast"]]

## export GMT
list2gmt <- function(list, description="") {
  res <- lapply(seq(along=list),
                function(x)
                list(name=names(list)[x],
                     description=description,
                     genes=as.character(list[[x]])))
  names(res) <- names(list)
  return(res)
}
cnames <- names(cpmSigGeneSymbol)
transformCnames <- function(x) {
  x <- gsub(" - ", " ", x)
  x <- gsub("\\(|\\)", "", x)
  x <- gsub(" ", "_", x)
  paste(x, "_NGS_GTEx_0.8_3",sep="")
}
ncnames <- transformCnames(cnames)
names(cpmSigGeneSymbol) <- names(cpmSigGeneID) <- ncnames
cpmGmt <- list2gmt(cpmSigGeneSymbol, description="Roche")
cpmGIGmt <- list2gmt(cpmSigGeneID, description="Roche")
write_gmt(cpmGmt, "data/exp.tissuemark.gtex.roche.symbols.gmt")
write_gmt(cpmGIGmt, "data/exp.tissuemark.gtex.roche.GeneID.progmt")


