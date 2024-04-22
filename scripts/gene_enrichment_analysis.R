library(clusterProfiler)
library(enrichplot)
all_genes <- read.csv("results/gene_enrichment/gene_to_KO.csv", row.names = 1)
tentacles <- read.csv("results/gene_enrichment/wild/wild_C_head.csv")

# KO enrichment up-regulated 
kegg_enrich <- enrichKEGG(gene = tentacles$KO,
                          universe = all_genes$KO,
                          organism = "ko",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

dotplot(kegg_enrich)
enrichplot(kegg_enrich)
gseaplot(kegg_enrich,
         title = "KEGG Pathway Enrichment Analysis",
         top_term = 10)

ecokd <- read.csv("results/deseq2_rsem/tables/EcoKD1_Eco1KD_B8_vs_control_B8_results.csv")

ecokd <- ecokd[ecokd$padj <= 0.05, ]
  
extract_function <- function(string) {
  parts <- unlist(strsplit(string, split = "\\."))
  return(parts[2])
}

ecokd$ID <- sapply(ecokd$X, extract_function)
merged_df <- merge(all_genes, ecokd, by="ID", all=TRUE)
merged_df <- na.omit(merged_df)

ecokd_up <- merged_df[merged_df$log2FoldChange >= 1.0, ]
ecokd_down <- merged_df[merged_df$log2FoldChange <= -1.0, ]


# KO enrichment up-regulated 
kegg_enrich <- enrichKEGG(gene = ecokd_up$KO,
                          universe = all_genes$KO,
                          organism = "ko",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

dotplot(kegg_enrich)

# KO enrichment down-regulated 
kegg_enrich <- enrichKEGG(gene = ecokd_down$KO,
                          universe = all_genes$KO,
                          organism = "ko",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

dotplot(kegg_enrich)

ecokd <- read.csv("results/gene_enrichment/ecokd1/ecokd1_L_ECTODERM.csv", row.names = 1)
wild <- read.csv("results/gene_enrichment/wild/wild_L_ECTODERM.csv", row.names = 1)
temp <- read.csv("results/gene_enrichment/temp/temp_L_ECTODERM.csv", row.names = 1)


kegg_enrich <- enrichKEGG(gene = ecokd$KO,
                          universe = all_genes$KO,
                          organism = "ko",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

dotplot(kegg_enrich)


kegg_enrich <- enrichKEGG(gene = wild$KO,
                          universe = all_genes$KO,
                          organism = "ko",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

dotplot(kegg_enrich)


kegg_enrich <- enrichKEGG(gene = temp$KO,
                          universe = all_genes$KO,
                          organism = "ko",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

dotplot(kegg_enrich)