id2rownames <- ribiosAnnotation:::id2rownames
ORACLE.BIN.USER <- ORACLE.BIN.PWD <- "genome"

annotateEnsemble <- function (ids, orthologue = FALSE, multiOrth = FALSE) 
{
  comm <- paste("SELECT e.ENSEMBL_ID, c.RO_GENE_ID,c.GENE_SYMBOL, c.DESCRIPTION, c.TAX_ID", 
                " FROM GTI_GENES c INNER JOIN ENSEMBL_GENE e ON c.RO_GENE_ID=e.GENE_ID ", sep = "")
  ann <- querydbTmpTbl(comm, "e.ENSEMBL_ID", ids, "bin",
                       ORACLE.BIN.USER, 
                       ORACLE.BIN.PWD)
  cnames <- c("ENSEMBL_ID", "GeneID", "GeneSymbol", "GeneName", "TaxID")
  conames <- c("OrigGeneID", "OrigGeneSymbol", "OrigGeneName", 
               "OrigTaxID")
  if (!orthologue) {
    colnames(ann) <- cnames
    cn <- "ENSEMBL_ID"
    res <- ann
  } else {
    colnames(ann) <- conames
    cn <- "OrigGeneID"
    ort <- annotateHumanOrthologsNoOrigTax(ann$OrigGeneID, 
                                           multiOrth = multiOrth)
    if (multiOrth) {
      res <- merge(ann, ort, by = "OrigGeneID", all.x = TRUE)
    }
    else {
      ort.re <- matchColumn(ann$OrigGeneID, ort, "OrigGeneID", 
                            multi = FALSE)
      res <- cbind(ann, ort.re[, -1L])
    }
    res <- putColsFirst(res, c("GeneID", "GeneSymbol", "TaxID", 
                               "OrigTaxID", "OrigGeneID", "OrigGeneSymbol", "OrigGeneName"))
  }
  res <- matchColumn(ids, res, cn, multi = orthologue && multiOrth)
  rownames(res) <- id2rownames(res[, cn])
  return(res)
}
