library(DESeq2)

data_path <- "./data/jacint_data/REG_RNASeq/reg2.gene.counts.matrix"
result_path <- "./results/jacint_results/reg_rnaseq/deseq2/"
count_threshold <- 20

counts_original <- read.table(data_path)
counts <- round(as.matrix(counts_original))

samples <- colnames(counts)
timepoint <- sub("^(T[0-9]+)_.*", "\\1", samples)
condition <- ifelse(grepl("_AS", samples), "AS", "S")

colData <- data.frame(
  sample = samples,
  timepoint = factor(timepoint, levels = c("T0","T1","T2","T3")),
  condition = factor(condition, levels = c("S","AS"))
)

rownames(colData) <- samples
head(colData)

results_list <- list()

for(tp in levels(colData$timepoint)){
  cat("[*] Running DESeq2 for", tp, "\n")
  
  # prepare metadata for specific timepoint experiment
  samples_tp <- rownames(colData)[colData$timepoint == tp]
  sub_counts = counts[,samples_tp]
  sub_colData = colData[samples_tp,]
  
  # read counts and remove low abundant transcripts
  dds <- DESeqDataSetFromMatrix(countData = sub_counts, colData = sub_colData, design = ~condition)
  keep <- rowSums(counts(dds)) >= count_threshold
  dds <- dds[keep,]
  # perform deseq analysis: AS vs S | treatment vs control
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition","AS","S"))
  # save results into results_list
  results_list[[tp]] <- res
  
  write.csv(as.data.frame(res), file = paste0(result_path, "DESeq2_", tp, "_AS_vs_S.csv"))
}

