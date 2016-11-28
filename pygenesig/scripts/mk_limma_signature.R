expr = npyLoad("../data/gtex/exprs_counts.npy", dotranspose = F)
colData = read_csv("../data/gtex/covariates.csv")
colData$SMTS = as.factor(colData$SMTS)
colData$Gender = as.factor(colData$Gender)

# RIN = RNA quality index
# SMTS = tissue
dgList = DGEList(counts=expr, genes=as.character(1:nrow(expr)))
designMat = model.matrix(~colData$SMTS + colData$RIN + colData$Gender)

# filter for cpm > 1 for at least 10 samples per gene
CPM = cpm(dgList)
cpmCnt = apply(CPM > 1, 1, sum)
cpmInds = cpmCnt >= 10
dgList.f = dgList[cpmInds,]

dgList.f = calcNormFactors(dgList.f)

v = voom(dgList.f, designMat, plot=FALSE)

fit = lmFit(v, designMat)
save(fit, file="limma_model.RData")
}

# apply cutoffs. 
fit.e = treat(fit, lfc=log2(100))
rslt = decideTests(fit.e, p.value=(0.01 / nrow(fit)))


```{r annotation}
ensembl = useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="www.ensembl.org")  
ensemble2hgnc = data.table(getBM(mart=ensembl, attributes=(c("ensembl_gene_id", "hgnc_symbol"))))
ensemble2hgnc = ensemble2hgnc[hgnc_symbol!='',]

# remove the .digit from the gene name. 
salitizeEnid = function(x) { 
  return(str_split(x, "\\.")[[1]][1])
}

# append sanitized EnsemblID to data table 
rslt.tab = data.table(rslt, keep.rownames=TRUE)
names = make.names(colnames(rslt.tab))
setnames(rslt.tab, names)
rn.2 = lapply(rslt.tab[,rn], salitizeEnid)
rslt.tab = rslt.tab[,rn2:=rn]

# append hgnc to data table 
rslt.matched = data.table(matchColumn(ensemble2hgnc$ensembl_gene_id, rslt.tab, 'rn2'))
rslt.matched = rslt.matched[,hgnc:=ensemble2hgnc[,hgnc_symbol]]
```

### Create signatures
Finally, we create lists of signatures for each tissue by collating the list of gene symbols that
are marked as overexpressed by *limma*. 
```{r create_signatures}
signatures = list()
sigCols =  names[grepl("^colData", names)]
sigCols = sigCols[1:(length(sigCols)-2)] # remove RIN and Gender
for(tCol in sigCols) {
  signatures = append(signatures, list(rslt.matched[get(tCol)==1,hgnc]))
}

names(signatures) = str_match(sigCols, "colData.SMTS(.*)")[,2]
```

### Export to gmt
```{r write_gmt, include=FALSE}
writeGmt = function(signatures, filename, description) { 
  if(file.exists(filename)) {
    file.remove(filename)
  }
  for(signame in names(signatures)) {
    sig = signatures[[signame]]
    genes = if(length(sig) > 0) str_c(sig, collapse="\t") else ""
    line = sprintf("GTEx_%s\t%s\t%s", signame, description, genes)
    write(line,file=filename,append=TRUE)
  }
}
writeGmt(signatures, "gmt/gtex_signatures.gmt", "GTEx_limma")