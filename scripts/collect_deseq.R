## standard in R 4, setting for R<4
options(stringsAsFactors=FALSE)

## options
## use deduplicated data?
use.dedup <- TRUE

### PATHS: adapt for run on different machine! ###

## input paths
datp <-  file.path("results/deseq2_rsem_final/tables/")

## create output path for results
outp <- file.path("results/deseq2_rsem_final/")
dir.create(outp, showWarnings=FALSE)
outfile <- file.path(outp,"hydra_all_deseq2_lg2.tsv")
outfilep <- file.path(outp,"hydra_all_deseq2_pvalues.tsv")
outfiler <- file.path(outp,"hydra_all_counts.tsv")

## get deseq2 result files
## TODO: use only those specified in samples_controls.csv!?
files <- list.files(pattern="results.csv$", path=datp)

datl <- datv <- nam <- NULL
for ( file in files ) {
  fdat <- read.csv(file.path(datp,file))
  nam[[file]] <- sub("hydra_","",sub("_results.csv","",file))
  
  ## TODO: filter here!
  
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

## collect raw data
rawp <- file.path("data/rsem_counts_final/")
files <- list.files(pattern="genes.results$", recursive=TRUE, path=rawp)
# files <- files[grep("kallistoCountsDedup",files, invert=!use.dedup)]

datr <- NULL
for ( file in files ) {
  id <- strsplit(file,"/")[[1]][2]
  id <- strsplit(id,'.genes')[[1]][[1]]
  tmp <- read.delim(file.path(rawp, file))
  colnames(tmp)[4:7] <- paste0(id,"_",colnames(tmp)[4:7])
  if ( is.null(datr) ) 
    datr <- tmp
  else {
    if ( any(datr[,1]!=tmp[,1]) ) {
      stop("wrong ids")
    }
    datr <- cbind(datr,tmp[4:7])
  }
}
write.table(datr,row.names=FALSE,
            quote=FALSE, file=outfiler, sep="\t")
