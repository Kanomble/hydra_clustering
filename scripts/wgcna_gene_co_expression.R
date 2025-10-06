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
combined_matrix[is.na(combined_matrix)] <- 0

# remove batch effects
counter <- 1
all_datasets <- c()
experiment_traits <- c()
for (exp in experiments) {
  dataset <- paste0("dataset", counter)
  
  for (col in samples[samples$EXP == exp, ]$treatment) {
    all_datasets <- c(all_datasets, dataset)
    experiment_traits <- c(experiment_traits,col)
  }
  
  counter <- counter + 1
}

########################################################
# WGCNA pipeline for your combined VST matrix (62 samples)
# - Assumes combined_matrix has genes in rows, samples in columns
# - Install packages if needed
########################################################

# ---------- 0. Install / load libraries ----------
if (!requireNamespace("WGCNA", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("WGCNA")
}
if (!requireNamespace("limma", quietly = TRUE)) install.packages("limma")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")

library(WGCNA)
library(limma)
library(data.table)

# Allow multithreading for speed
allowWGCNAThreads()

# ---------- 1. Load your data ----------
# Option A: if your matrix is already in the R environment as 'combined_matrix', skip this.
# Option B: read from file (uncomment and edit the file path)
# combined_matrix <- fread("path/to/combined_vst.csv", data.table = FALSE)
# rownames(combined_matrix) <- combined_matrix[,1]; combined_matrix <- combined_matrix[,-1]

# Check dims
cat("combined_matrix dims (genes x samples):", dim(combined_matrix), "\n")
# Expect something like: genes x 62

# ---------- 2. (Optional but recommended) Define dataset labels / batch vector ----------
# Provide a vector of length = number of samples (62), with dataset IDs, e.g. "set1","set2","set3","set4"
# Two options:
#  - Provide as a vector here:
# dataset_labels <- c(rep("A", 15), rep("B", 20), rep("C", 14), rep("D", 13))
#
#  - OR provide a CSV with a column "Sample" and column "Dataset" with sample names matching colnames(combined_matrix)
# traits_df <- read.csv("path/to/sample_metadata.csv", stringsAsFactors = FALSE)
# dataset_labels <- traits_df$Dataset[match(colnames(combined_matrix), traits_df$Sample)]
#
# If you do not have labels, set dataset_labels <- NULL
dataset_labels <- experiment_traits   # <--- REPLACE NULL with your vector or load from file as above
colnames(combined_matrix) <- dataset_labels
combined_matrix <- combined_matrix[, !(colnames(combined_matrix) %in% c("control_B8", "control", "18°C"))]
dataset_labels <- dataset_labels[!dataset_labels %in% c("control_B8", "control", "18°C")]
colnames(combined_matrix) <- dataset_labels

# If user provided labels, validate length
if (!is.null(dataset_labels)) {
  if (length(dataset_labels) != ncol(combined_matrix)) stop("dataset_labels length != number of samples (columns). Fix it and re-run.")
  dataset_labels <- as.factor(dataset_labels)
  cat("Dataset labels provided with levels:", levels(dataset_labels), "\n")
}

# ---------- 3. (Optional) Batch correction ----------
# If dataset_labels is provided and you want to remove dataset-specific batch effects, use removeBatchEffect.
# NOTE: removeBatchEffect works on genes in rows and sample-level batch vector.
combined_matrix_corrected <- combined_matrix
#if (!is.null(dataset_labels)) {
#  cat("Performing batch correction with limma::removeBatchEffect using dataset labels...\n")
#  combined_matrix_corrected <- removeBatchEffect(as.matrix(combined_matrix), batch = dataset_labels)
#} else {
#  cat("No dataset labels provided. Skipping batch correction.\n")
#}
combined_matrix_corrected <- combined_matrix
# ---------- 4. Prepare datExpr for WGCNA (samples in rows) ----------
datExpr0 <- t(as.data.frame(combined_matrix_corrected))  # now samples x genes
# make sure samples are rows and genes are columns
cat("datExpr dims (samples x genes):", dim(datExpr0), "\n")

# ---------- 5. Filter genes by variance (recommended for large combined data) ----------
# Keep top variable genes (customize topN). For very large gene sets, using top 5000-15000 variable genes is common.
topN <- 10000   # change to 5000 or 15000 depending on memory/time
geneVars <- apply(datExpr0, 2, var, na.rm = TRUE)
ord <- order(geneVars, decreasing = TRUE)
topGenes <- ord[1:min(topN, length(ord))]
datExpr <- datExpr0[, topGenes]
cat("Filtered to top variable genes:", ncol(datExpr), "\n")

# ---------- 6. Check for good samples/genes ----------
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodSamples) > 0) {
    cat("Removing bad samples:\n"); print(rownames(datExpr)[!gsg$goodSamples])
    datExpr <- datExpr[gsg$goodSamples, ]
  }
  if (sum(!gsg$goodGenes) > 0) {
    cat("Removing bad genes: ", sum(!gsg$goodGenes), "\n")
    datExpr <- datExpr[, gsg$goodGenes]
  }
}

# ---------- 7. Sample clustering to detect outliers ----------
sampleTree <- hclust(dist(datExpr), method = "average")
# Save tree plot
png("./results/network_analysis/sampleClustering.png", width = 1200, height = 800)
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "")
if (!is.null(dataset_labels)) {
  # Add colored annotation bar under dendrogram
  # create a color vector for dataset labels
  labelColors <- as.character(as.numeric(as.factor(dataset_labels)))
  # simple approach: plot colored rectangles (not perfect if labels not matching sample order)
}
dev.off()
cat("Saved sampleClustering.png\n")

# If you want to manually remove outlier samples, do it here (script leaves all samples in by default)

# ---------- 8. Pick soft-thresholding power ----------
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Print the table similar to your earlier output
sftTable <- data.frame(Power = sft$fitIndices[,1],
                       SFT.R.sq = signif(sft$fitIndices[,2], 4),
                       slope = signif(sft$fitIndices[,3], 4),
                       truncated.R.sq = signif(sft$fitIndices[,5], 4),
                       mean.k. = signif(sft$fitIndices[,6], 4),
                       median.k. = signif(sft$fitIndices[,7], 4),
                       max.k. = signif(sft$fitIndices[,8], 4))
print(sftTable)

# Plot the results
png("./results/network_analysis/scaleFreePlot.png", width = 1400, height = 800)
par(mfrow = c(1,2))
cex1 <- 0.8
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2",
     type="n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")
abline(h=0.80, col="blue", lty=2) # target 0.8

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()
cat("Saved scaleFreePlot.png\n")

# Choose softPower:
# Use automatic rule: pick smallest power with SFT.R.sq >= 0.8 if it exists, otherwise fallback to where curve flattens.
softPower <- NA
idx <- which(sft$fitIndices[,2] >= 0.8)
if (length(idx) > 0) {
  softPower <- sft$fitIndices[min(idx), 1]
  cat("Selected softPower based on R^2 >= 0.8:", softPower, "\n")
} else {
  # fallback: choose power at which mean connectivity drops but not too low
  # as you observed earlier, sometimes no power reaches 0.8. Choose the power with highest SFT.R.sq.
  bestRow <- which.max(sft$fitIndices[,2])
  softPower <- sft$fitIndices[bestRow, 1]
  cat("No power reaches R^2 >= 0.8. Using power with highest SFT.R.sq:", softPower, "\n")
}
cat("softPower =", softPower, "\n")

# ---------- 9. Construct network and detect modules (blockwiseModules for large datasets) ----------
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
write.csv(moduleAssignment, "./results/network_analysis/module_assignment.csv", row.names = FALSE)
cat("Saved module_assignment.csv\n")

# ---------- 10. Merge modules (if not already merged by blockwiseModules) ----------
# blockwiseModules merges modules; here we ensure we have mergedMEs
mergedMEs <- MEs   # if you used mergeCloseModules you could replace it

# ---------- 11. Module-eigengene relationships with dataset labels (if available) ----------
if (!is.null(dataset_labels)) {
  # Prepare trait dataframe: dataset membership as binary indicators (one-hot)
  traitDF <- model.matrix(~0 + dataset_labels)
  colnames(traitDF) <- sub("^dataset_labels", "", colnames(traitDF))
  
  # correlation
  moduleTraitCor <- cor(mergedMEs, traitDF, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
  
  # Save matrices
  write.csv(moduleTraitCor, "./results/network_analysis/moduleTraitCor.csv")
  write.csv(moduleTraitPvalue, "./results/network_analysis/moduleTraitPvalue.csv")
  cat("Saved moduleTraitCor.csv and moduleTraitPvalue.csv\n")
  
  # Heatmap
  png("./results/network_analysis/module_trait_heatmap.png", width = 1200, height = 800)
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
  write.csv(moduleEigengeneCor, "./results/network_analysis/moduleEigengeneCor.csv")
  png("./results/network_analysis/module_eigengene_heatmap.png", width = 1000, height = 800)
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
  fname <- paste0("./results/network_analysis/genes_module_", mod, ".txt")
  write.table(genes_in_module, file = fname, quote = FALSE, row.names = FALSE, col.names = FALSE)
}
cat("Exported gene lists for each module (genes_module_<color>.txt)\n")

# ---------- 13. Save workspace and important objects ----------
save(net, moduleAssignment, mergedMEs, file = "./results/network_analysis/WGCNA_results.RData")
cat("Saved WGCNA_results.RData\n")

cat("WGCNA pipeline finished. Review PNGs and CSVs for results.\n")

# Make sure you have both
head(moduleLabels)   # numeric
head(moduleColors)   # color names
colnames(mergedMEs)  # should be like "ME0", "ME1", ...

# Build mapping table
label2color <- data.frame(
  Label = sort(unique(moduleLabels)),
  Color = labels2colors(sort(unique(moduleLabels)))
)

# Now relabel eigengenes
ME_names <- colnames(mergedMEs)
ME_numbers <- as.numeric(gsub("ME", "", ME_names))  # extract numbers
color_for_ME <- label2color$Color[match(ME_numbers, label2color$Label)]

color2ME <- data.frame(
  ME = ME_names,
  Color = color_for_ME
)

print(color2ME)


moduleGenes <- colnames(datExpr)[moduleColors == "blue"]
datExprModule <- datExpr[, moduleGenes]
# Build TOM similarity matrix
TOM <- TOMsimilarityFromExpr(datExprModule, power = softPower)
dissTOM <- 1 - TOM

# Build igraph network
graph <- graph.adjacency(as.matrix(dissTOM), mode="undirected", weighted=TRUE, diag=FALSE)

# Optionally filter weak edges
graph <- delete_edges(graph, E(graph)[weight < 0.2])

# Set vertex labels to gene names
V(graph)$name <- moduleGenes

# Plot
plot(graph, vertex.label=V(graph)$name, vertex.size=5, edge.width=E(graph)$weight*2,
     main="Gene-level Network for Blue Module (ME1)")

