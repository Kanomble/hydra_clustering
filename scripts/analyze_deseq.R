### ANALYZE DESEQ2 RESULTS ###
## analyze transcriptomes from four experiments (Bathia/Fraune) by
## clustering from Siebert et al. 2019, mapped to transcript
## models by blasts

## TODO:
## * optimize clustering?

library(segmenTools)
library(viridis)
library(umap)
options(stringsAsFactors=FALSE)

### PATHS ###

## input paths
resp <- file.path("results/deseq2_rsem_final/") ## deseq2 results
figp <- file.path("results/figures/normalized_mean_final/")
datp <- file.path("results/processedData/normalized_mean_final/")

expf <- file.path("data/rsem_counts_final/sample_descriptions.csv")
rawf <- file.path("results/deseq2_rsem_final/hydra_all_counts_t_to_g.tsv") #result/deseq2_rsem/ # _t_to_g
datf <- file.path("results/deseq2_rsem_final/hydra_all_deseq2_lg2_t_to_g.tsv") # 
pvlf <- file.path("results/deseq2_rsem_final/hydra_all_deseq2_pvalues_t_to_g.tsv") # _t_to_g

## annotation
cltf <- file.path("./results/processedData/normalized_mean_final/hvaep_clustering_table.tsv") #cltf <- file.path(datp,"siebert19_clustering_table.tsv") # cell types
gof <- file.path("./results/processedData/normalized_mean_final/hvaep_uniprot_kegg_go.tsf") #gof <- file.path(datp,"aepLRv2_kegg_go.tsv") # KEGG/GO

#anf <- file.path("./data/new_mapping/") #anf <- file.path(dath, "originalData", "aepLRv2_pfam_sp_nr.tsv") # Siebert et al

### ANALYSIS PARAMETERS ###

## data set filter: remove some redundant comparisons
## data order: 
selected.data <- c("Wild_vs_GF","Cvbct_vs_GF","Conventionalized_vs_GF",
                   "3OHC12_vs_control","3OC12_vs_control","22.C_vs_18.C","12.C_vs_18.C","08.C_vs_18.C","Eco1KD_B8_vs_control_B8") # "22.C_vs_18.C","12.C_vs_18.C","08.C_vs_18.C","Eco1KD_B8_vs_control_B8","3OC12_vs_control"
#"08.C_vs_18.C", ,"Cvbct_vs_GF","Conventionalized_vs_GF"
#"12.C_vs_18.C", "12.C_vs_18.C","08.C_vs_18.C"
#"22.C_vs_18.C"
#"Wild_vs_GF","Cvbct_vs_GF"                   
#
## SET THIS to FALSE if you use different combinations of controls!
## TODO: use this also for other plots!
short.class.names <- TRUE 

## data filter: cells with cell-type annotations, etc?
filter <- "gfzero" #"gfzero"#"gfzero" # "celltype" #  "outlier" # 
filt.nap <- TRUE# FALSE ## filter 1.8k genes with p_adj=NA in all experiments

## celltypes: only analyze unique gene-cell type combis
unique.celltype <- TRUE #FALSE #FALSE #TRUE # FALSE #

## filter for up/down classes
lg2.cutoff <- 1
p.cutoff <- 0.05

## volcano plot cutoffs
maxx <- 4 # 8
maxy <- 10


## parse data
dat <- read.delim(datf)
pvl <- read.delim(pvlf)
clt <- read.delim(cltf) # sieber et al. cell types
got <- read.delim(gof)
#ann <- read.delim(anf)

# set all values lg2fold changes of dat < -15 and > 15 to NA
#identifier <- dat$ID
#dat[, -1][abs(dat[, -1]) >= 15] <- NA
#dat$ID <- identifier
#pvl[is.na(dat)] <- NA



## clean names: remove experiment names from column names
colnames(dat) <- sub(".*?_","",colnames(dat))
#colnames(pvl) <- sub(".*?_","",colnames(pvl))

## clean IDs in Siebert et al. supplement
#ann[,1] <- sub("\\|.*","", ann[,1])
#colnames(ann)[1] <- "ID"

## check identity
if ( any(dat$ID!=pvl$ID) ) stop("wrong IDs in p-value file")

## expand Annotation data to RNAseq data file
clte <- clt[match(dat$ID, clt$ID),]
gote <- got[match(dat$ID, got$ID),]
#anne <- ann[match(dat$ID, ann$ID),### FILTER DATA ###
K <- 6

## custom filter 
RM <- logical()
if ( filter=="ct_or_go" ) {
  ## rm: not present in any annotation file
  RM <- is.na(clte$ID) & is.na(gote$ID)
} else if ( filter=="celltype" ) {
  ## rm: not present in cell type annotation
  RM <- is.na(clte$ID)
  K <- 6 #6 # 8 # cluster number
} else if ( filter=="outlier" ) {
  ## rm: present in cell type annotation
  RM <- !is.na(clte$ID) 
  maxx <- 8
} else if ( filter=="gfzero" ) {
  
  ## rm: 0 counts in germ-free conditions
  
  ## get experiment sample IDs
  experiments <- read.csv(expf)
  gfid <- experiments$ID[experiments$treatment=="GF"]
  ## get raw data
  raw <- read.delim(rawf)
  raw <- raw[match(dat$ID, raw$gene_id),]
  ## transcript lengths
  taep <-raw[,1:3]
  ## germ free: GF experiments
  raw <- raw[,paste0(sub("-S1",".S1",gfid),"_expected_count")]
  RM <- apply(raw,1,sum)==0
}


## SELECTED DATA, incl. order 
datf <- dat[,selected.data]
colnames(pvl) <- colnames(dat)
pvlf <- pvl[,selected.data]
ids <- dat$ID

if ( any(RM) ) { #18638
  
  cat(paste("removing", sum(RM) , "features with filter:\t", filter,"\n",
            "KEEPING", sum(!RM), "features for cluster analysis\n\n"))
  
  ## filter all data
  ids <-   ids[!RM]
  datf <- datf[!RM,]
  pvlf <- pvlf[!RM,]
  clt <-  clte[!RM,]
  got <-  gote[!RM,]
  #ann <-  anne[!RM,]
}


if ( short.class.names ) {
  warning("using short class names; this can cause errors ",
          "when combining different controls")
  colnames(datf) <- sub("_vs_.*","", colnames(datf))
  colnames(pvlf) <- sub("_vs_.*","", colnames(pvlf))
}

## FURTHER FILTERS

## remove all with p.value
if ( filt.nap ) {
  
  rmp <- which(apply(is.na(pvlf),1,sum)>=1) #ncol(pvlf)-1 --> removing all columns with p_value = NA
  
  cat(paste("removing", length(rmp) , "features with NA p-values\n",
            "KEEPING", nrow(pvlf)-length(rmp),
            "features for cluster analysis\n\n\n\n"))
  
  ids <-   ids[-rmp]
  datf <- datf[-rmp,]
  pvlf <- pvlf[-rmp,]
  clt <-   clt[-rmp,]
  got <-   got[-rmp,]
  #ann <-   ann[-rmp,]
}

### COMPARISON TO CELL TYPES ###

## p values cutoffs for plots
cp.min <- 1e-20
cp.txt <- 1e-10
cp.srt <- 5*cp.txt

## FILTER CELL TYPE ANNOTATION:

## super categories present?
## TODO: do this smarter, e.g. store super clustering in separate table
## or use prefix super:?

## tag cells w/o annotation - TODO:after removing non-unique??
cna <- apply(clt[,2:ncol(clt)],1,function(x) !any(!is.na(x)))

## remove annotations that are NOT unique to a cell type, per hierachy level
cltu <- clt

## super clusters present?
scidx <- min(grep("^L_", colnames(clt)))
pridx <- min(grep("^C_", colnames(clt)))

## column indexes for separate categories
upos <- sort(unique(c(2,scidx,pridx,1+ncol(clt))))
urng <- list()
for ( j in seq_along(upos[-1]) ) 
  urng[[j]] <- upos[j]:(upos[j+1]-1)

for ( j in seq_along(urng) ) {
  ## NOT unique to cell type
  ## print(colnames(clt[, urng[[j]]]))
  ## print(1)
  rmc <- apply(clt[,urng[[j]]],1,sum,na.rm=TRUE) > 2
  cltu[rmc,urng[[j]]] <- NA
}

if ( unique.celltype ) clt <- cltu

## test: 412
sum(cltu$C_neuron,na.rm=TRUE)


outfile_uc <- paste0(datp,"unique_celltypes.table")
write.table(x=cltu,file=outfile_uc,quote = FALSE, sep = ';')

## prepare T/F table from Siebert et al. clustering
clt.table <- clt[,2:(ncol(clt))]
clt.table[is.na(clt.table)] <- 0
clt.table <- apply(clt.table, 2, as.logical)

## use only cell types with at least 1 unique gene
## sum over columns and check if summation value is greater than 0
clt.table <- clt.table[,apply(clt.table,2,sum) > 0]

## add "NA"
## TODO: tag non-annotated cells AFTER removing non-unique?
## or add NU class
clt.table <- cbind(clt.table, "NA"=cna)

## t-test profile
ovt <- clusterProfile(x=datf, cls=clt.table)
ov.stat <- ovt$statistic # store for clustering
ovt$statistic[] <- round(ov.stat) # round for plot

## sort
ovs <- sortOverlaps(ovt, p.min=cp.srt, cut=FALSE) # 1e-25
## sort & cut
ovsc <- sortOverlaps(ovt, p.min=cp.srt, cut=TRUE) # 1e-25

## split by up/down - hack: sort separately and fuse order
ovsup <- ovsdo <- ovs
ovsup$p.value[ovs$sign==-1] <- 1
ovsup <- sortOverlaps(ovsup, p.min=cp.srt, cut=TRUE) # 1e-25
ovsdo$p.value[ovs$sign==1] <- 1
ovsdo <- sortOverlaps(ovsdo, p.min=cp.srt, cut=TRUE) # 1e-25

nsrt <- unique(c(names(ovsup$num.query[,1]),names(ovsdo$num.query[,1])))

ovsup <- sortOverlaps(ovs, srt=names(ovsup$num.query[,1]))
ovsdo <- sortOverlaps(ovs, srt=names(ovsdo$num.query[,1]))
ovsall <- sortOverlaps(ovs, srt=nsrt)

statistic_file <- paste0(datp, "ttest_all_statistic.table")
pval_file <- paste0(datp, "ttest_all_pvalue.table")
target_number_file <- paste0(datp, "ttest_all_target_number.table")
query_number_file <- paste0(datp, "ttest_all_query_number.table")

write.table(x=ovt$statistic,file=statistic_file)
write.table(x=ovt$p.value,file=pval_file)
write.table(x=ovt$num.target,file=target_number_file)
write.table(x=ovt$num.query,file=query_number_file)

## order by clustering

## hclust - cell types
dst <- dist(ov.stat)
clh <- hclust(dst)
new.srt <- clh$order

## hclust - experiments
dst <- dist(t(ov.stat))
clh <- hclust(dst)
new.col <- clh$order

## k-means - cell types
##clk <- kmeans(ov.stat, center=3)
###new.srt <- order(clk$cluster)

## sort t-test profiles by hclust order
ovl <- ovll <- sortOverlaps(ovt, srt=new.srt)

## sort columns in ovll by hclust order
for (i in 1:length(ovll))
  if (class(ovll[[i]]) == "matrix") 
    if (ncol(ovll[[i]]) == ncol(ovt$p.value)) 
      ovll[[i]] <- ovll[[i]][, new.col, drop = FALSE]


statistic_file_sc <- paste0(datp, "ttest_sorted_cut_statistic.table")
pval_file_sc <- paste0(datp, "ttest_sorted_cut_pvalue.table")
target_number_file_sc <- paste0(datp, "ttest_sorted_cut_target_number.table")
query_number_file_sc <- paste0(datp, "ttest_sorted_cut_query_number.table")

## plot cell type profiles
write.table(x=ovsc$statistic,file=statistic_file_sc)
write.table(x=ovsc$p.value,file=pval_file_sc)
write.table(x=ovsc$num.target,file=target_number_file_sc)
write.table(x=ovsc$num.query,file=query_number_file_sc)

## bare w/o legend
png(file.path(figp,"celltype_ttest_sorted_cut.png"),res=100, units="in",width=ncol(ovsc$p.value)/3+2, height=nrow(ovsc$p.value)/4+1.5)
par(mai=c(1.5,1.9,0.1,.6), tcl=-.25, mgp=c(1.3,.4,0))
plotOverlaps(ovsc,p.min=cp.min, p.txt=cp.txt, show.total=TRUE, text.cex=.8,
             xlab=NA, ylab=NA)
dev.off()

png(file.path(figp,"celltype_ttest_updown_cut.png"), res=100, units="in",
    width=ncol(ovsc$p.value)/3+2, height=nrow(ovsc$p.value)/4+1.5)
par(mai=c(1.5,1.9,0.1,.6), tcl=-.25, mgp=c(1.3,.4,0))
plotOverlaps(ovsall,p.min=cp.min, p.txt=cp.txt, show.total=TRUE, text.cex=.8,
             xlab=NA, ylab=NA)
dev.off()

## full, unsorted
png(file.path(figp,"celltype_ttest.png"), res=100, units="in",
    width=ncol(ovt$p.value)/3+5, height=nrow(ovt$p.value)/4+1.5)
layout(matrix(1:4,ncol=2),
       widths=c(2, ncol(ovt$p.value)/3),
       heights=c(nrow(ovt$p.value)/4,1.5))
par(mgp=c(.5,.4,0), tcl=-.25)
plot.new()
par(mai=c(.25,.25,0,.01))
plotOverlapsLegend(p.min=cp.min, p.txt=cp.txt, type=2, dir=2, las=1)
par(mai=c(.2,.15,0.01,.5))
plotOverlaps(ovt,p.min=cp.min, p.txt=cp.txt, show.total=TRUE, text.cex=.8,
             xlab=NA, ylab=NA)
dev.off()

png(file.path(figp,"celltype_ttest_sorted.png"), res=100, units="in",
    width=ncol(ovs$p.value)/3+2, height=nrow(ovs$p.value)/4+1.5)
layout(matrix(1:4,ncol=2),
       widths=c(2, ncol(ovs$p.value)/3),
       heights=c(nrow(ovs$p.value)/4,1.5))
par(mgp=c(.5,.4,0), tcl=-.25)
plot.new()
par(mai=c(.25,.25,0,.01))
plotOverlapsLegend(p.min=cp.min, p.txt=cp.txt, type=2, dir=2, las=1)
par(mai=c(.2,.15,0.01,.5))
plotOverlaps(ovs,p.min=cp.min, p.txt=cp.txt, show.total=TRUE, text.cex=.8,
             xlab=NA, ylab=NA)
dev.off()

png(file.path(figp,"celltype_ttest_clustered.png"), res=100, units="in",
    width=ncol(ovl$p.value)/3+2, height=nrow(ovl$p.value)/4+1.5)
layout(matrix(1:4,ncol=2),
       widths=c(2, ncol(ovl$p.value)/3),
       heights=c(nrow(ovl$p.value)/4,1.5))
par(mgp=c(.5,.4,0), tcl=-.25)
plot.new()
par(mai=c(.25,.25,0,.01))
plotOverlapsLegend(p.min=cp.min, p.txt=cp.txt, type=2, dir=2, las=1)
par(mai=c(.2,.15,0.01,.5))
plotOverlaps(ovl,p.min=cp.min, p.txt=cp.txt, show.total=TRUE, text.cex=.8,
             xlab=NA, ylab=NA)
dev.off()



png(file.path(figp,"celltype_ttest_clustered_both.png"), res=100, units="in",
    width=ncol(ovll$p.value)/3+2, height=nrow(ovll$p.value)/4+1.5)
layout(matrix(1:4,ncol=2),
       widths=c(2, ncol(ovll$p.value)/3),
       heights=c(nrow(ovll$p.value)/4,1.5))
par(mgp=c(.5,.4,0), tcl=-.25)
plot.new()
par(mai=c(.25,.25,0,.01))
plotOverlapsLegend(p.min=cp.min, p.txt=cp.txt, type=2, dir=2, las=1)
par(mai=c(.2,.15,0.01,.5))
plotOverlaps(ovll,p.min=cp.min, p.txt=cp.txt, show.total=TRUE, text.cex=.8,
             xlab=NA, ylab=NA)
dev.off()

## WILCOX TESTS - SLOW, confirms t-test profiles
if ( TRUE ) {
  
  w.test <- function(x,y) {
    res <- wilcox.test(x,y)
    ## normalized U-statistic
    tt <- res$statistic/(sum(!is.na(x))*sum(!is.na(y))) -0.5
    rt <- list()
    rt$statistic <- unlist(tt)
    rt$p.value <- unlist(res$p.value)
    rt
  }
  
  ovw <- clusterProfile(x=datf, cls=clt.table, test=w.test)
  ov.stat <- ovw$statistic
  ovw$statistic[] <- round(ov.stat,1)
  
  ## hclust 
  dst <- dist(ov.stat)
  clh <- hclust(dst)#, center=3)
  #new.srt <- clh$order
  
  ## k-means
  clk <- kmeans(ov.stat, center=3)
  ##new.srt <- order(clk$cluster)
  
  ## sort
  ovws <- sortOverlaps(ovw, srt=new.srt)
  
  WIDTHS <- c(.3,.7)
  
  png(file.path(figp,"celltype_wilcox.png"), res=100, units="in", width=7, height=7)
  layout(matrix(1:4,ncol=2),  widths=WIDTHS, heights=c(.8,.25))
  par(mgp=c(.5,.4,0), tcl=-.25)
  plot.new()
  par(mai=c(.25,.25,0,.01))
  plotOverlapsLegend(p.min=1e-50, p.txt=1e-25, type=2, dir=2, las=1)
  par(mai=c(.2,.15,0.01,.5))
  plotOverlaps(ovw,p.min=1e-50, p.txt=1e-25, show.total=TRUE, text.cex=.8,
               xlab=NA, ylab=NA)
  dev.off()
  
  png(file.path(figp,"celltype_wilcox_clustered.png"), res=100, units="in",
      width=7, height=7)
  layout(matrix(1:4,ncol=2), widths=WIDTHS, heights=c(.8,.25))
  par(mgp=c(.5,.4,0), tcl=-.25)
  plot.new()
  par(mai=c(.25,.25,0,.01))
  plotOverlapsLegend(p.min=1e-50, p.txt=1e-25, type=2, dir=2, las=1)
  par(mai=c(.2,.15,0.01,.5))
  plotOverlaps(ovws,p.min=1e-50, p.txt=1e-25, show.total=TRUE, text.cex=.8,
               xlab=NA, ylab=NA)
  dev.off()
}

### CLUSTERING et al.: cross-correlation, clustering, etc. ###

## set NA to LFC=0 and p=1 for algorithms to work
datz <- datf
datz[is.na(datz)] <- 0
pvlz <- pvlf
pvlz[is.na(pvlz)] <- 1

cat(paste("setting", sum(is.na(datf)), "NA to 0 for algorithms\n"))

## CROSS-CORRELATE EXPERIMENTS

plotdev(file.path(figp,"cormatrix_experiments"), type="png", res=100, width=5, height=5)
par(mai=c(2.2,2.2,.1,.1), mgp=c(1,.4,0), tcl=-.25)
image_matrix(cor(datz),col=viridis::viridis(100), axis=1:2, xlab=NA, ylab=NA)
dev.off()

if ( FALSE ) {
  ## by genes - large and not useful without sorting
  image_matrix(cor(t(datz)),col=viridis::viridis(100))
}


## CLUSTER TRANSCRIPTOME

## by cutoffs
cll <- sapply(1:ncol(datf),function(i) {
  
  nm <- colnames(datz)[i]
  x <- datz[,i]
  p <- pvlz[,i]
  cl <- rep("non",length(x))
  cl[x>  lg2.cutoff] <- "up"
  cl[x< -lg2.cutoff] <- "down"
  cl[which(p > p.cutoff)] <- "non"
  cl[is.na(p)] <- "non"
  cl
  ##paste0(nm,"_",cl)
})
colnames(cll) <- paste0(colnames(datf))#,"_class")


## by kmeans
set.seed(1)
clk <- kmeans(datz, center=K)
cls.srt <- 1:K 
cls <- clk$cluster # clustering!
cls.col <- cls.srt
names(cls.col) <- cls.srt

write.table(cls,paste0(datp,"kmeans_on_datz.table"))
write.table(datz,paste0(datp,"kmeans_on_datz.table"))
write.table(clt.table,paste0(datp,"clt.table"))
## test clustering alternatives
if ( FALSE ) {
  
  ## test flow clust - NOTE: all BIC are NaN
  clf <- flowClust::flowClust(datz, K=3:10)
  ## inspect BIC 
  lapply(clf, function(x) x@BIC)
  
  ## TODO: test mclust results?
  ## test mclust - NOTE: best BIC always for max K
  library(mclust)
  clm <- mclust::Mclust(datz, G=3:7)
  ## inspect BIC
  clm$BIC
  plot(clm, what="BIC") # TODO: use plateau-method to select?
  
  ## M3C - TAKES TOO LONG, run at hpc?
  pdf("m3c_output.pdf")
  clc <- M3C::M3C(data.frame(t(datz)), clusteralg="pam", maxK=9)
  dev.off()
}



## T-TEST PROFILES
ovt <- clusterProfile(x=datf, cls=factor(cls,levels=cls.srt))
ov.stat <- ovt$statistic # store for clustering
ovt$statistic[] <- round(ov.stat) # round for plot

## transpose
ovtt <- lapply(ovt, function(x) if ( class(x)=="matrix" ) t(x) else x)
ovtt$num.query <- ovt$num.target
ovtt$num.target <- ovt$num.query

png(file.path(figp,"cluster_ttest.png"), res=100, units="in", width=5, 
    height=nrow(ovtt$p.value)/4+1)
par(mai=c(0.5,2.2,.5,.5), mgp=c(1,.4,0), tcl=-.25)
plotOverlaps(ovtt,p.min=1e-50, p.txt=1e-10, show.total=TRUE, text.cex=1,
             ylab=NA, xlab="kmeans clusters",axis1.col=cls.col)
dev.off()

## OVERLAP WITH CELL TYPES CLUSTERS
ovc <- clusterAnnotation(cls=cls, cls.srt=cls.srt,
                         data=clt.table)
ovc <- sortOverlaps(ovc, p.min=0.001)

png(file.path(figp,"cluster_celltype.png"), res=100, units="in", width=5, 
    height=nrow(ovc$p.value)/4)
par(mai=c(.5,2.5,.5,.5), mgp=c(1,.4,0), tcl=-.25, las=2)
plotOverlaps(ovc, p.min=1e-20, p.txt=1e-10,ylab=NA,show.total=TRUE,text.cex=.7,
             xlab="kmeans clusters", axis1.col=cls.col)
dev.off()

ovcc <- sortOverlaps(ovc, p.min=5e-5, cut=TRUE)

png(file.path(figp,"cluster_celltype_cut.png"), res=100, units="in", width=5, 
    height=nrow(ovcc$p.value)/5+1)
par(mai=c(.5,2.5,.5,.5), mgp=c(1,.4,0), tcl=-.25, las=2)
plotOverlaps(ovcc, p.min=1e-20, p.txt=1e-10,ylab=NA,show.total=TRUE,text.cex=.7,
             xlab="kmeans clusters", axis1.col=cls.col)
dev.off()

## CLUSTERS vs. CUTOFF CLASSES

tsrt <- paste0(sub("_vs_.*","",selected.data), "_",
               rep(c("up","non","down"),each=length(selected.data)))

## generate TRUE/FALSE table from cutoff classes
enms <- colnames(cll)
clln <- sapply(1:ncol(cll), function(j) paste0(enms[j],"_",cll[,j]))
cllt <- cbind(ids,apply(clln,1, function(x) paste(x,collapse=";")))
cllt <- parseAnnotationList(cllt)


## vs themselves
tovl <- annotationOverlap(cllt)

if ( short.class.names ) { # manual sorting
  tovl <- sortOverlaps(tovl, srt=tsrt, axis=2)
  tovl <- sortOverlaps(tovl, srt=tsrt, axis=1)
} else { # p-value sorting
  tovl <- sortOverlaps(tovl, p.min=1e-20, cut=FALSE)
  tovl <- sortOverlaps(tovl, srt=rownames(tovl$num.query), axis=1)
}
plotdev(file.path(figp,"classes"),type="png", width=ncol(tovl$p.value)/3+1,
        height=nrow(tovl$p.value)/3+1)
par(mai=c(2,2,.6,.6))
plotOverlaps(tovl, p.min=1e-20, p.txt=1e-10,
             show.total=TRUE, text.cex=.5, xlab=NA, ylab=NA)
dev.off()

## vs kmeans
lgc <- clusterAnnotation(cls=cls, data=cllt)
lgcs <- sortOverlaps(lgc, p.min=0.001, cut=FALSE)
plotdev(paste0(file.path(figp,"cluster_classes")), type="png",
        width=ncol(lgcs$p.value)/3+2,
        height=nrow(lgcs$p.value)/4+1)
par(mai=c(.5,2,.5,.5), mgp=c(1,.4,0), tcl=-.25, las=2)
plotOverlaps(lgcs, p.min=1e-10, p.txt=1e-5,
             xlab="kmeans clusters", ylab=NA,show.total=TRUE,text.cex=.7,
             axis1.col=cls.col)
dev.off()
lgcs <- sortOverlaps(lgc, srt=tsrt)
plotdev(paste0(file.path(figp,"cluster_classes_sorted")), type="png",
        width=ncol(lgcs$p.value)/3+2,
        height=nrow(lgcs$p.value)/4+1)
par(mai=c(.5,2,.5,.5), mgp=c(1,.4,0), tcl=-.25, las=2)
plotOverlaps(lgcs, p.min=1e-10, p.txt=1e-5,
             xlab="kmeans clusters", ylab=NA,show.total=TRUE,text.cex=.7,
             axis1.col=cls.col)
dev.off()

## TODO: loop over all experiments (ncol(cll)) for more focussed plots

if ( FALSE ) {
  eid <- "Cvbct"
  
  ecl <- cll[,eid]
  lgc <- clusterAnnotation(cls=ecl, data=cllt)
  lgcs <- sortOverlaps(lgc, p.min=0.001, cut=FALSE)
  plotdev(paste0(file.path(figp,eid,"_overlap")), type="png", width=5, height=8)
  par(mai=c(.7,3.5,.5,.5), mgp=c(1.3,.4,0), tcl=-.25, las=2)
  plotOverlaps(lgcs, p.min=1e-10, p.txt=1e-5,
               xlab="", ylab=NA,show.total=TRUE,text.cex=.7)
  mtext(sub("_vs_.*","",eid), side=1, line=2.5, las=1, cex=1.2)
  dev.off()
  
  ## profile of max response clusters in all experiments
  ovt <- clusterProfile(x=datf, cls=ecl)
  ovtt <- lapply(ovt, function(x) if ( class(x)=="matrix" ) t(x) else x)
  ovtt$num.query <- ovt$num.target
  ovtt$num.target <- ovt$num.query
  
  plotdev(paste0(file.path(figp,paste0(eid,"_goiprofile"))), type="png", width=5, height=5)
  par(mai=c(.7,3.5,.5,.5), mgp=c(1.3,.4,0), tcl=-.25, las=2)
  plotOverlaps(ovtt, p.min=1e-10, p.txt=1e-5, txt.col=NA,
               xlab="", ylab=NA,show.total=TRUE,text.cex=.7)
  mtext(sub("_vs_.*","",eid), side=1, line=2.5, las=1, cex=1.2)
  dev.off()
}




## QC : VOLCANO PLOTS colored by clustering
for ( j in 1:ncol(datf) ) {
  
  exp <- colnames(datz)[j]
  
  ## cut axes for more informative plots
  ptmp <- -log10(pvlz[,j])
  ptmp[ptmp>maxy] <- maxy
  ltmp <- datz[,j]
  idx <- which(abs(ltmp)>maxx)
  ltmp[idx] <- sign(ltmp[idx])*maxx
  exp <- gsub("\\.","-", exp)
  plotdev(file.path(paste0(figp,"volcano_",exp)), type="png", width=3.5, height=3.5)
  par(mai=c(.5,.5,0.25,.1), mgp=c(1.3,.4,0), tcl=-.25)
  plot(0,col=NA, xlab="log2 fold-change",
       ylab=expression(-log[10](p)), main=exp, xlim=c(-maxx,maxx),
       ylim=c(0,maxy))
  abline(h=c(0,maxy), col=8, lwd=.7)
  abline(v=c(-maxx,0,maxx), col=8, lwd=.7)
  rect(xleft=-maxx, ybottom=-log10(p.cutoff), xright=-lg2.cutoff, ytop=10,
       col="#AAAAAA",border=NA)
  rect(xleft=lg2.cutoff, ybottom=-log10(p.cutoff), xright=maxx, ytop=10,
       col="#AAAAAA",border=NA)
  points(ltmp, ptmp, col=cls, pch=19, cex=.25)
  dev.off()
}


### ANNOTATION ENRICHMENT ANALYSES ###

gp.min <- 1e-10
gp.txt <- 1e-5

## get go term mapping
gtrms <- lapply(strsplit(got$GO,";"),trimws)
gtrms <- unique(unlist(gtrms))
##gtrms <- gtrms[!is.na(gtrms)]
gtrms <- lapply(strsplit(gtrms, " \\[GO:"),
                function(x) {x[2] <- paste0("GO:",x[2]);x})
gtrms <- do.call(rbind, gtrms)
gtrms[,2] <- sub("\\]$","",gtrms[,2])
terms <- gtrms[,1]
names(terms) <- gtrms[,2]
gtrms <- terms
gtrms[is.na(gtrms)] <- names(gtrms)[is.na(gtrms)] <- "NA"


## get GO terms for genes in list
gol <- got$GO # [match(clt$ID, got$ID)]
gol <- lapply(strsplit(gol,";"),trimws)
gol <- lapply(gol, function(x) gsub("\\]$","",gsub(".*\\[GO:","GO:", x)))
gol <- lapply(gol, function(x) { x[is.na(x)] <- "NA"; x})

## construct T/F table
gotab <- matrix(FALSE, nrow=nrow(got), ncol=length(gtrms))
colnames(gotab) <- names(gtrms)
for ( i in 1:length(gol) ) 
  gotab[i,gol[[i]]] <- TRUE

## GO ANALYSIS
ovg <- clusterAnnotation(cls=cls, cls.srt=cls.srt,
                         data=gotab, terms=gtrms, replace.terms =TRUE,
                         verbose=FALSE)

for ( ex in -3:-8 ) {
  
  ovgc <- sortOverlaps(ovg, p.min=10^ex, cut=TRUE)
  png(paste0(file.path(figp,"cluster_GO_1e"),ex,".png"), res=100, units="in", width=6,
      height=nrow(ovgc$p.value)/5 +1)
  par(mai=c(.5,3.5,.5,.5), mgp=c(1,.4,0), tcl=-.25, las=2)
  plotOverlaps(ovgc, p.min=gp.min, p.txt=gp.txt,
               xlab="kmeans clusters", ylab=NA,show.total=TRUE,text.cex=.7,
               axis1.col=cls.col)
  dev.off()
}



## CLUSTERING vs. UMAP, PCA?

## use umap?
ump <- umap::umap(datz)
plotdev(file.path(figp,"cluster_umap"), type="png", width=3.5, height=3.5)
par(mai=c(.5,.5,0.1,.1), mgp=c(1.3,.4,0), tcl=-.25)
plot(ump$layout, col=cls,pch=19, cex=.5)
dev.off()

## use PCA?
pca <- prcomp(t(datz))
plotdev(file.path(figp,"cluster_pca"), type="png", width=3.5, height=3.5)
par(mai=c(.5,.5,0.1,.1), mgp=c(1.3,.4,0), tcl=-.25)
plot(ash(pca$rotation[,1:2]), col=cls,pch=19, cex=.5)
dev.off()

## PCA for experiments
pca <- prcomp(datz)
plotdev(file.path(figp,"experiments_pca"), type="png", width=3.5, height=3.5)
par(mai=c(.5,.5,0.1,.1), mgp=c(1.3,.4,0), tcl=-.25)
plot(pca$rotation[,1:2], pch=19, cex=.5)
text(x=pca$rotation[,1], y=pca$rotation[,2], labels=colnames(datz))
dev.off()



### WRITE RESULTS FILE - lg2 fold changes, clusterings, annotations ##

## expression data
tmp <- datf
colnames(tmp) <- paste0(colnames(tmp),"_log2ratio")
tpp <- pvlf
colnames(tpp) <- paste0(colnames(tpp),"_p_adjusted")

## cell type data
#ctd <- sapply(1:nrow(clt.table), function(i) {
#  paste0(colnames(clt.table)[clt.table[i,]],collapse=";")
#})
#ctd[ctd=="NA"] <- ""
#
#results <- data.frame(ID=ids, tmp, tpp, cluster=cls,
#                      cll,
#                      celltype=ctd,
#                      got[,c("UniProtKB","KEGG.ID")],
#                      ann[,c("Pfam_annot","SP_annot","nr_annot")])
#resfile <- file.path(datp, paste0("clusterings_",filter,".tsv"))
#write.table(results, file=resfile,
#            sep="\t", quote=FALSE, row.names=FALSE, na="")
