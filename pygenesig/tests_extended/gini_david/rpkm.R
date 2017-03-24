library(ribiosUtils)
library(ribiosIO)
library(ribiosAnnotation)
library(ribiosPlot)

geneRpkms <- read_gct_matrix("data/GTEx_Analysis_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm__Pilot_V3_patch1.gct")
sp <- read.table("data/GTEx_Analysis_Annotations_Sample_DS__Pilot_V3.txt",sep="\t", head=TRUE, quote="", comment.char="")

## pheno Data
pheno <- matchColumn(colnames(geneRpkms), sp, "SAMPID")
pheno$SMTS <- factor(pheno$SMTS)
pat <- read.table("data/GTEx_Analysis_Annotations_Subject_DS__Pilot_V3.txt",sep="\t", head=TRUE, quote="")
sampleSub <- sapply(strsplit(as.character(pheno$SAMPID), "-"), function(x) paste(x[1], x[2], sep="-"))
sampleSub[sampleSub=="NA12878-SM"] <- "NA12878_C"
sampleSub[sampleSub=="NA12878_C-SM"] <- "NA12878_C"
stopifnot(all(sampleSub %in% pat[,1]))
pheno$SUBJID <- sampleSub
pheno <- merge(pheno, pat, by="SUBJID")

## make a similar plot as on the GTEx website
ind <- grep("ENSG00000101986", rownames(geneRpkms))
indExprs <- geneRpkms[ind,]
indFac <- names(sort(tapply(indExprs, pheno$SMTSD, median), decreasing=TRUE))
cols <- rainbow(length(indFac))
par(mar=c(18,4,1,1))
boxplot(log10(indExprs)~factor(pheno$SMTSD, levels=indFac), las=3L,
        ylab="log10(RPKM)", col=cols, main="ABCD1")
abline(h=log10(1), lty=2, lwd=2)
boxplot(log10(indExprs)~factor(pheno$SMTSD, levels=indFac), las=3L,
        ylab="log10(RPKM)", col=cols, main="ABCD1", add=TRUE)
ipdf("ABCD1_GTex_Roche.pdf", height=12, width=10)

## translate GENE names
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl",
                   useMart("ensembl"))
ensRoot <- gsub("\\.[0-9]*$", "", rownames(geneRpkms))
geneRpkmSyms <- attr(geneRpkms, "desc")
## by symbol 
human <- gtiTaxAnnotation(9606)
feat.raw <- matchColumn(geneRpkmSyms,
                    human, "GeneSymbol")

test <- geneRpkms[with(feat.raw, is.na(GeneID)),]
badHigh <- apply(test, 1L, median)
badGenes <- subset(feat.raw[is.na(feat.raw$GeneID),][which(badHigh>1),1:3],
                  !grepl("\\.|-", GeneSymbol))$GeneSymbol
writeLines(badGenes,
           "problemCases.txt")

## Load into GTI, map alias to GeneSymbols
fix <- read.table("problemCases-GTI.xls",sep="\t", header=TRUE, quote="", comment.char="")
stopifnot(!any(duplicated(fix[,1])))
fixId <- match(fix[,1], geneRpkmSyms)
geneRpkmSyms[fixId] <- fix$Gene.Symbol
feat <- matchColumn(geneRpkmSyms,
                        human, "GeneSymbol")
feat$GTex_EnsEMBL_GeneID <- ensRoot
feat <- feat[,c("GeneID", "GeneSymbol", "GeneName", "GeneType", "GTex_EnsEMBL_GeneID")]

##----------------------------------------##
## Export
##----------------------------------------##
write_gct(geneRpkms,
          "data/RocheAnnotated_GTEx_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm__Pilot_2013_01_31.gct",
          feat.name=feat$GeneSymbol,
          feat.desc=feat$GeneID)
write.table(feat,
            "data/RocheAnnotated_GTEx_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm__Pilot_2013_01_31-featureAnnotation.txt",
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(pheno,
            "data/RocheAnnotated_GTEx_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm__Pilot_2013_01_31-sampleAnnotation.txt",
            sep="\t", row.names=FALSE, quote=FALSE)

stop()

##----------------------------------------##
## Data for UDIS: only GeneID
##----------------------------------------##
library(ribiosExpression)
rpkms <- read_gct_matrix("data/RocheAnnotated_GTEx_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm__Pilot_2013_01_31.gct")
feat <- read.table("data/RocheAnnotated_GTEx_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm__Pilot_2013_01_31-featureAnnotation.txt",
                   header=TRUE, sep="\t", comment.char="", quote="")

isEntrez <- !is.na(feat$GeneID)
entRpkms <- rpkms[isEntrez,]
entFeat <- feat[isEntrez,]
rownames(entRpkms) <- rownames(entFeat) <- entFeat$GTex_EnsEMBL_GeneID

entEset <- new("ExpressionSet",
               exprs=entRpkms, featureData=new("AnnotatedDataFrame", entFeat))
entEsetSumm <- summarizeProbesets(entEset, index.name="GeneID", fun=sum, keep.nonindex=FALSE, keep.featureName=FALSE)

writeGct(entEsetSumm,
          "udis/RocheAnnotated_GTEx_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm__Pilot_2013_01_31-GeneID.gct",
          feat.name="GeneID", feat.desc="GeneSymbol")

##humanGenes <- gtiTaxAnnotation(9606)
##hOut <- cbind(ID=humanGenes$GeneID, EG_ID=humanGenes$GeneID)
##hgFile <- "udis/HumanGeneID.annot"
##writeLines("!platform_table_begin", hgFile)
##write.table(hOut, hgFile, sep="\t", append=TRUE, quote=FALSE, row.names=FALSE)
##cat()
## by EnsEMBL/BioMart
##featFile <- "BioMart-Ensembl-EntrezID.txt"
##if(file.exists("BioMart-Ensembl-EntrezID.txt")) {
##  feat <- read.table(featFile)
##} else {
##  feat.raw <- getBM(filters="ensembl_gene_id",
##                attributes=c("ensembl_gene_id", "entrezgene"),
##                values=ensRoot,
##                mart=mart)
##  feat <- matchColumn(ensRoot, feat.raw, "ensembl_gene_id")
##  colnames(feat) <- c("ensembl_gene_id", "GeneID")
##  write.table(feat, file=featFile)
##}
#### compare
##isNA <- is.na(feat$GeneID) | is.na(idBySym$GeneID)
##diffID <- which(feat$GeneID != idBySym$GeneID & !isNA)
##diffLines <- data.frame(Ens=feat[diffID, "ensembl_gene_id"],
##                        BioMart=feat[diffID, "GeneID"],
##                        GTI_ID=idBySym$GeneID[diffID],
##                        Symbol=geneRpkmSyms[diffID],
##                        GTI=idBySym$GeneSymbol[diffID])
