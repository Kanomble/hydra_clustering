## parse and process clustering from suppl.table S9 from Siebert et al. 2019
## scRNAseq of stem cell differentiation trajectories in Hydra magnipapillata

## TODO: analyze overlaps between clusters
## TODO: generate TRUE/FALSE table of unique clusterings
## TODO: use swissprot/PFAM/NR annotations to get GO or
##       other systematic annotation
## TODO: add subclustered interstitial and neuronal cells on sheets 2 and 3

library(readxl)
library(segmenTools)

options(stringsAsFactors=FALSE)

# INPUT
# cls.file <- file.path("./data/hvaep_cell_type_to_gene_cluster_table_just_gos.tsf")
#cls.file <- file.path("./data/hvaep_cell_type_to_gene_cluster_table.tsf")
#hvaep_cell_type_to_gene_cluster_table.tsf
cls.file <- file.path("./results/processedData/normalized_mean_final/hvaep_cell_type_to_gene_cluster_table.tsf")
gof <- file.path("./results/processedData/normalized_mean_final/hvaep_uniprot_kegg_go.tsf")

## OUTPUT
out.path <- file.path("./results/processedData/normalized_mean_final/")
fig.path <- file.path("./results/figures/normalized_mean_final/")

### PARAMETERS
## add super classes
add.super <- TRUE # FALSE #  

## read supplemental file
clsf <- read.csv(cls.file, sep="\t")
## investigate
table(clsf[,1])

clsf <- data.frame(clsf)

## get transcript IDs
tids <- sub("\\|.*","",clsf[,2])
sum(duplicated(tids)) # NOTE: many duplicates - no unique clustering

## export
res <- cbind.data.frame(ID=tids, CLS=clsf[,1])
out.file <- file.path(out.path, "hvaep_clustering_normalized.tsv")
write.table(res, file=out.file, sep="\t",row.names=FALSE, quote=FALSE)


## analyze overlaps
tids <- unique(res$ID)
cids <- unique(res$CLS)

## generate T/F table
ctab <- matrix(FALSE, nrow=length(tids), ncol=length(cids))
rownames(ctab) <- tids
colnames(ctab) <- cids
for ( j in 1:ncol(ctab) ) 
  ctab[,j] <- tids%in%res$ID[res$CLS==cids[j]]

## add super-categories
if ( add.super ) {
  super <- list()

  super[["C_neuron"]] <- grep("N$",cids, value=TRUE) # 8 
  super[["C_germline"]] <- grep("GC",cids, value=TRUE) # 2
  super[["C_nematoblast"]] <- grep("NB$",cids, value=TRUE) # 3
  super[["C_nematocytes"]] <- c(grep("NC$", cids, value=TRUE), grep("I_EarlyNem",cids, value=TRUE)) # --> this is not in the superclass: I_EarlyNem --> it has to be included later
  super[["C_gland_cell"]] <- c(grep("GL$",cids, value=TRUE),
                              grep("I_SpumMucGl",cids, value=TRUE)) # I_GranGl I_ZymoGl I_SpumMucGl | excluded for now: I_GlProgen 
  
  super[["C_progenitor"]] <- grep("I_GlProgen",cids, value=TRUE)
  super[["C_neuroprogenitor"]] <- grep("I_Neuro", cids, value=TRUE)
  super[["C_EC_foot"]] <- c(grep("Ec_BasalDisk", cids, value=TRUE), grep("Ec_Peduncle", cids, value=TRUE))
  super[["C_EN_foot"]] <- grep("En_Foot", cids, value=TRUE)
  super[["C_tentacle"]] <- grep("_Tentacle$", cids, value=TRUE)
  super[["L_foot"]] <- c(grep("Ec_BasalDisk", cids, value=TRUE), 
                         grep("Ec_Peduncle", cids, value=TRUE), grep("En_Foot", cids, value=TRUE))
  super[["C_head"]] <- grep("Head$", cids, value=TRUE)
  
  ## stem cell lineages # first level of hierachy
  super[["L_ENDODERM"]] <- grep("En_BodyCol/SC", cids, value=TRUE)
  super[["L_ECTODERM"]] <- grep("Ec_BodyCol/SC", cids, value=TRUE)
  super[["L_INTERSTITIAL"]] <- grep("I_ISC",cids, value=TRUE)
  
  # added drop = FALSE if there is just one column in x
  supm <- do.call(cbind,lapply(super, function(x) apply(ctab[, x, drop = FALSE], 1, any)))
  ctab <- cbind(ctab,supm)
}

## store T/F table (slim format: 1 for TRUE, empty for FALSE)
out.file <- file.path(out.path,"hvaep_clustering_table.tsv")
clt <- apply(ctab,2,as.numeric)
clt[clt==0] <- NA
write.table(cbind.data.frame(ID=rownames(ctab),clt),
            file=out.file, sep="\t",row.names=FALSE, quote=FALSE, na="")

## analyze internal overlaps in a T/F clustering table
ovl <- annotationOverlap(ctab)

## overlap count
plotdev(file.path(fig.path,"hvaep_clustering"),type="png",
        width=11, height=11)
par(mai=c(2.5,2.5,.6,.6))
plotOverlaps(ovl, p.min=1e-250, p.txt=1e-125, n=100,
             show.total=TRUE, text.cex=.5, xlab=NA, ylab=NA)
dev.off()

## UNIQUE GENES per cell type and lineage categories
cltu <- clt

## super clusters present?
scidx <- min(grep("^L_", colnames(clt)))
pridx <- min(grep("^C_", colnames(clt)))

## column indexes for separate categories
upos <- sort(unique(c(1,scidx,pridx,1+ncol(clt))))
urng <- list()
for ( j in seq_along(upos[-1]) ) 
  urng[[j]] <- upos[j]:(upos[j+1]-1)

for ( j in seq_along(urng) ) {
  ## NOT unique to cell type
  rmc <- apply(clt[,urng[[j]]],1,sum,na.rm=TRUE) > 1#1 #4
  cltu[rmc,urng[[j]]] <- NA
}

out.file <- file.path(out.path,"hvaep_clustering_table_unique.tsv")
write.table(cbind.data.frame(ID=rownames(ctab),cltu),
            file=out.file, sep="\t",row.names=FALSE, quote=FALSE, na="")


### ANALYZE ANNOTATION ###

## parse all GO terms
got <- read.delim(gof, stringsAsFactors=FALSE)

## filter available?
filter.go <- TRUE
if ( filter.go )
  got <- got[match(rownames(ctab), got$ID),]

## get go term mapping
#gtrms <- lapply(strsplit(got$GO.terms,";"),trimws)
gtrms <- lapply(strsplit(got$GO,";"),trimws)
gtrms <- unique(unlist(gtrms))
gtrms <- lapply(strsplit(gtrms, " \\[GO:"),
                function(x) {x[2] <- paste0("GO:",x[2]);x})
gtrms <- do.call(rbind, gtrms)
gtrms[,2] <- sub("\\]$","",gtrms[,2])
terms <- gtrms[,1]
names(terms) <- gtrms[,2]
gtrms <- terms

# Dont overwrite the first entry with NA but extend dataframe with NA
gtrms <- c("NA", gtrms)
names(gtrms)[1] <- gtrms[1] <-"NA"

## get GO terms for genes in list
gol <- got$GO[match(rownames(ctab), got$ID)]
gol <- lapply(strsplit(gol,";"),trimws)
gol <- lapply(gol, function(x) gsub("\\]$","",gsub(".*\\[GO:","GO:", x)))
gol <- lapply(gol, function(x) { x[is.na(x)] <- "NA"; x})

## construct T/F table
gotab <- matrix(FALSE, nrow=nrow(ctab), ncol=length(gtrms))
colnames(gotab) <- names(gtrms)
for ( i in 1:length(gol) ){
  gotab[i,gol[[i]]] <- TRUE  # test gol[[1]] %in% gtrms[1]
}


## prepare T/F table - TODO: it was a T/F table above, avoid back and forth
clt.table <- cltu
clt.table[is.na(clt.table)] <- 0
clt.table <- apply(clt.table, 2, as.logical)

## use only cell types with at least 1 unique gene
clt.table <- clt.table[,apply(clt.table,2,sum) > 0]


## super clusters present?
scidx <- min(grep("^L_", colnames(clt.table)))
pridx <- min(grep("^C_", colnames(clt.table)))

## column indexes for separate categories
upos <- sort(unique(c(1,scidx,pridx,1+ncol(clt.table))))
urng <- list()
for ( j in seq_along(upos[-1]) ) 
  urng[[j]] <- upos[j]:(upos[j+1]-1)

## collapse to clusterings
ctc <- list()
for ( j in seq_along(urng) ) {
  ctd <- sapply(1:nrow(clt.table), function(i) {
    nm <- colnames(clt.table[,urng[[j]]])[clt.table[i,urng[[j]]]]
    if ( length(nm)==0) nm <- "non-unique"
    if ( length(nm)>1) nm <- "ERROR" # should be unique!
    nm
  })
  ctc[[j]] <- ctd
}
ctc <- do.call(cbind, ctc)

ovg <- clusterAnnotation(cls=ctc[,3], 
                         data=gotab, terms=gtrms, replace.terms=TRUE,
                         verbose=FALSE)
ovgc <- sortOverlaps(ovg, p.min=5e-7, cut=TRUE)

png(file.path(fig.path,"hvaep_GO_lineage_I_superclasses.png"),
    res=100, units="in",
    width=ncol(ovgc$p.value)/4 + 4,
    height=nrow(ovgc$p.value)/4 + 2)
par(mai=c(1.5,3.5,.5,.5), mgp=c(1,.4,0), tcl=-.25, las=2)
plotOverlaps(ovgc, p.min=1e-10, p.txt=1e-5,
             xlab="", ylab=NA,show.total=TRUE,text.cex=.75)
dev.off()

ovg <- clusterAnnotation(cls=ctc[,2], 
                         data=gotab, terms=gtrms, replace.terms=TRUE,
                         verbose=FALSE)
ovgc <- sortOverlaps(ovg, p.min=5e-7, cut=TRUE)
png(file.path(fig.path,"hvaep_GO_lineage_C_superclasses.png"),
    res=100, units="in",
    width=ncol(ovgc$p.value)/4 + 4,
    height=nrow(ovgc$p.value)/4 + 2)
par(mai=c(1.5,3.5,.5,.5), mgp=c(1,.4,0), tcl=-.25, las=2)
plotOverlaps(ovgc, p.min=1e-10, p.txt=1e-5,
             xlab="", ylab=NA,show.total=TRUE,text.cex=.75)
dev.off()


ovg <- clusterAnnotation(cls=ctc[,1], 
                         data=gotab, terms=gtrms, replace.terms=TRUE,
                         verbose=FALSE)
ovgc <- sortOverlaps(ovg, p.min=1e-4, cut=TRUE)

png(file.path(fig.path,"hvaep_GO_celltype.png"),
    res=100, units="in",
    width=ncol(ovgc$p.value)/4 +3.5+.5,
    height=nrow(ovgc$p.value)/4 +1.5+.5)
par(mai=c(1.5,3.5,.5,.5), mgp=c(1,.4,0), tcl=-.25, las=2)
plotOverlaps(ovgc, p.min=1e-20, p.txt=1e-10,
             xlab="", ylab=NA,show.total=TRUE,text.cex=.75)
dev.off()

## REVERSE: check "structural constituent of the ribosome" vs. cell types

ridx <- names(which(gtrms=="structural constituent of ribosome"))
rcls <- rep("other",nrow(ctc))
rcls[gotab[,ridx]] <- "ribosome"

clt.table <- clt
clt.table[is.na(clt.table)] <- 0
clt.table <- apply(clt.table, 2, as.logical)

ovr <- clusterAnnotation(cls=rcls, data=clt.table, verbose=TRUE)
png(file.path(fig.path,"hvaep_ribosomes.png"),
    res=100, units="in",width=5,height=nrow(ovr$p.value)/4+1)
par(mai=c(1.5,3.5,.5,.5), mgp=c(1,.4,0), tcl=-.25, las=2)
plotOverlaps(ovr, p.min=1e-20, p.txt=1e-10,
             xlab="", ylab=NA,show.total=TRUE,text.cex=.75)
dev.off()


## GO ANALYSIS for clustering of all cell types

## cluster genes by presence in cell types

K=6
clk <- kmeans(ctab, center=K)
clk.srt <- 1:K

ovg <- clusterAnnotation(cls=clk$cluster, cls.srt=clk.srt,
                         data=gotab, terms=gtrms, replace.terms =TRUE,
                         verbose=FALSE)
ovgc <- sortOverlaps(ovg, p.min=0.0001, cut=TRUE)

png(file.path(fig.path,"hvaep_GO.png"),
    res=100, units="in",width=6,height=nrow(ovgc$p.value)/4)
par(mai=c(.5,3.5,.5,.5), mgp=c(1,.4,0), tcl=-.25, las=2)
plotOverlaps(ovgc, p.min=1e-20, p.txt=1e-10,
             xlab="kmeans clusters", ylab=NA,show.total=TRUE,text.cex=.5)
dev.off()