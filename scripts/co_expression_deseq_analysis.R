library(tximport)
library(DESeq2)
library(DEGreport)
library(limma)
library(pheatmap)
library(igraph)
library(RColorBrewer)

# define input file paths
sample_file <- file.path("data/rsem_counts_after_sortmerna/sample_descriptions.csv") # "data/rsem_counts/sample_descriptions_without_temp_oxo.csv"
controls_file <- file.path("data/rsem_counts_after_sortmerna/sample_controls.csv") # "data/rsem_counts/sample_controls_without_temp_oxo.csv"

input_dir <- "data/rsem_counts_after_sortmerna/"
output_dir <- "results/deseq2_rsem_sortmerna/"

# settings
count_threshold <- 20

# read samples file
samples <- read.csv(sample_file)
controls <- read.csv(controls_file)

# extract unique experiments for looping
experiments <- unique(samples$EXP)

expr_list <- list() # empty list


for(exp in experiments){
  # get experimental RSEM count file paths
  expSamples <- samples[samples$EXP == exp,]
  files <- file.path(input_dir, exp, paste0(expSamples$ID,".genes.results"))
  names(files) <- expSamples$ID
  # check if all files exist
  all(file.exists(files))
  # recommended tximport vignette cmd
  txi.rsem <- tximport(files, type="rsem", txIn = FALSE, txOut = FALSE)
  # create sampleTable for deseq2
  sampleTable <- data.frame(condition = factor(expSamples$treatment), batch = factor(expSamples$batch))
  rownames(sampleTable) <- colnames(txi.rsem$counts)
  # deseq2 workflow
  if(exp=="HydraRecolonization"){
    dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~batch+condition)
  } else {
    dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition)
  }
  
  keep <- rowSums(counts(dds)) >= count_threshold
  dds <- dds[keep,]
  
  if(exp == "HydraAHL"){
    cat("[*] Releveling To Real Reference\n")
    dds$condition <- relevel(dds$condition, ref = "control")
  } 
  
  if(exp == "EcoKD1"){
    cat("[*] Releveling To Real Reference\n")
    dds$condition <- relevel(dds$condition, ref = "control_B8")
  }
  
  keep <- rowSums(counts(dds)) >= count_threshold
  dds <- dds[keep,]
  dds <- DESeq(dds)
  
  # variance stabilizing transformation
  vsd <- vst(dds, blind=FALSE)
  expr_matrix <- assay(vsd)
  expr_list <- append(expr_list, list(expr_matrix)) 
}

expr_list <- lapply(expr_list, function(mat) {
  df <- as.data.frame(mat)        # ensure it’s a data frame
  df$Gene <- rownames(df)         # add rownames as a column
  df
})

combined_matrix <- Reduce(function(x, y) {
  merge(x, y, by = "Gene", all = TRUE)
}, expr_list)

rownames(combined_matrix) <- combined_matrix$Gene
combined_matrix$Gene <- NULL
#combined_matrix[is.na(combined_matrix)] <- 0

# remove batch effects
counter <- 1
all_datasets <- c()
experiment_traits <- c()
for (exp in experiments) {
  dataset <- paste0("dataset", counter)
  
  for (col in samples[samples$EXP == exp, ]$ID) {
    all_datasets <- c(all_datasets, dataset)
    experiment_traits <- c(experiment_traits,exp)
  }
  
  counter <- counter + 1
}

combined_matrix_corrected <- removeBatchEffect(combined_matrix, batch = all_datasets)
datExpr <- t(combined_matrix_corrected)
# Define traits as dataset labels
traits <- data.frame(Dataset = factor(experiment_traits))
rownames(traits) <- rownames(datExpr)  # sample names must match



BiocManager::install("WGCNA")
library(WGCNA)
allowWGCNAThreads()

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering", sub="", xlab="")
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot scale-free topology fit
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit", type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")
abline(h=0.90, col="blue") # Aim for R^2 cutoff ~0.9

softPower <- 6

adjacency <- adjacency(datExpr, power = softPower)

# Turn adjacency into topological overlap matrix (TOM)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1-TOM

# Hierarchical clustering of genes
geneTree <- hclust(as.dist(dissTOM), method="average")
plot(geneTree, xlab="", sub="", main="Gene clustering on TOM", labels=FALSE, hang=0.04)

# Module detection
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

table(dynamicMods)

dynamicColors <- labels2colors(dynamicMods)
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

# Cluster module eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of module eigengenes")

# Merge close modules
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.25, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

# Correlation among module eigengenes
moduleEigengeneCor <- cor(mergedMEs, use = "p")
moduleEigengenePvalue <- corPvalueStudent(moduleEigengeneCor, nrow(datExpr))

# Heatmap of module–eigengene relationships
labeledHeatmap(Matrix = moduleEigengeneCor,
               xLabels = names(mergedMEs),
               yLabels = names(mergedMEs),
               ySymbols = names(mergedMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = signif(moduleEigengeneCor, 2),
               main = "Module eigengene correlations")

# Adjust memory and block size parameters if needed
net <- blockwiseModules(datExpr,
                        power = softPower,
                        TOMType = "unsigned",
                        minModuleSize = 30,
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM",
                        verbose = 5)

# Module labels and colors
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
MEs <- net$MEs
cat("Detected modules (number):", length(unique(moduleColors)), "\n")

# Save module assignment table
geneIDs <- colnames(datExpr)
moduleAssignment <- data.frame(Gene = geneIDs, Module = moduleColors, Label = moduleLabels)
write.csv(moduleAssignment, "module_assignment.csv", row.names = FALSE)
cat("Saved module_assignment.csv\n")

# ---------- 10. Merge modules (if not already merged by blockwiseModules) ----------
# blockwiseModules merges modules; here we ensure we have mergedMEs
mergedMEs <- MEs   # if you used mergeCloseModules you could replace it


dataset_labels <- experiment_traits

if (!is.null(dataset_labels)) {
  # Prepare trait dataframe: dataset membership as binary indicators (one-hot)
  traitDF <- model.matrix(~0 + dataset_labels)
  colnames(traitDF) <- sub("^dataset_labels", "", colnames(traitDF))
  
  # correlation
  moduleTraitCor <- cor(mergedMEs, traitDF, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
  
  # Save matrices
  write.csv(moduleTraitCor, "moduleTraitCor.csv")
  write.csv(moduleTraitPvalue, "moduleTraitPvalue.csv")
  cat("Saved moduleTraitCor.csv and moduleTraitPvalue.csv\n")
  
  # Heatmap
  png("module_trait_heatmap.png", width = 1200, height = 800)
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(traitDF),
                 yLabels = colnames(mergedMEs),
                 ySymbols = colnames(mergedMEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = signif(moduleTraitCor, 2),
                 main = "Module-trait relationships (datasets)")
  dev.off()
  cat("Saved module_trait_heatmap.png\n")
} else {
  cat("dataset_labels not provided: skipping module vs dataset correlations. You can still inspect modules and export gene lists.\n")
  
  # Optionally create module vs module heatmap (eigengene correlations)
  moduleEigengeneCor <- cor(mergedMEs, use = "p")
  write.csv(moduleEigengeneCor, "moduleEigengeneCor.csv")
  png("module_eigengene_heatmap.png", width = 1000, height = 800)
  labeledHeatmap(Matrix = moduleEigengeneCor,
                 xLabels = colnames(mergedMEs),
                 yLabels = colnames(mergedMEs),
                 ySymbols = colnames(mergedMEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = signif(moduleEigengeneCor, 2),
                 main = "Module eigengene correlations")
  dev.off()
  cat("Saved module_eigengene_heatmap.png\n")
}


# ---------- 12. Export gene lists for each module (for enrichment) ----------
modules <- unique(moduleColors)
for (mod in modules) {
  genes_in_module <- geneIDs[moduleColors == mod]
  fname <- paste0("genes_module_", mod, ".txt")
  write.table(genes_in_module, file = fname, quote = FALSE, row.names = FALSE, col.names = FALSE)
}
cat("Exported gene lists for each module (genes_module_<color>.txt)\n")

# ---------- 13. Save workspace and important objects ----------
save(net, moduleAssignment, mergedMEs, file = "WGCNA_results.RData")
cat("Saved WGCNA_results.RData\n")

cat("WGCNA pipeline finished. Review PNGs and CSVs for results.\n")
