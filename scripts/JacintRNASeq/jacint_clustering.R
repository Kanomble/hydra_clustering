## standard in R 4, setting for R<4
options(stringsAsFactors=FALSE)

## input paths
datp <-  file.path("results/jacint_results/reg_rnaseq/orthofinder_ogs_deseq/")

## create output path for results
outp <- file.path("results/jacint_results/reg_rnaseq/")

outfile <- file.path(outp,"hydra_all_deseq2_lg2.tsv")
outfilep <- file.path(outp,"hydra_all_deseq2_pvalues.tsv")

## get deseq2 result files
files <- list.files(pattern="MergedOGs.csv$", path=datp)

datl <- datv <- nam <- NULL
for ( file in files ) {
  fdat <- read.csv(file.path(datp,file))
  nam[[file]] <- sub("_MergedOGs\\.csv$", "", file)
  
  datl[[file]] <- cbind.data.frame(ID=fdat$X,FC=fdat$log2FoldChange)
  datv[[file]] <- cbind.data.frame(ID=fdat$X,P=fdat$padj)
}

## collect IDs
ids <- unique(unlist(lapply(datl, function(x) x$ID)))

## full matrix
mdat <- matrix(NA, nrow=length(ids), ncol=length(nam))
colnames(mdat) <- nam
rownames(mdat) <- ids

pdat <- mdat

for ( j in 1:ncol(mdat) ) {
  mdx <- match( datl[[j]]$ID, ids)
  mdat[mdx, j] <- datl[[j]]$FC
  pdat[mdx, j] <- datv[[j]]$P
}

write.table(cbind.data.frame(ID=rownames(mdat),mdat), row.names=FALSE,
            quote=FALSE, file=outfile, sep="\t")
write.table(cbind.data.frame(ID=rownames(pdat),pdat), row.names=FALSE,
            quote=FALSE, file=outfilep, sep="\t")