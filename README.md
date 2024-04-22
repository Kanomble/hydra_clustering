# Hydra Clustering
deseq2_rsem_script.R: results/deseq2_rsem/tables/*, results/deseq2_rsem/figures/* produces log2FoldChange and padj result tables and PCA, Volcano-Plot figures
collect_deseq.R: results/deseq2_rsem/: hydra_all_deseq2_lg2.tsv, hydra_all_deseq2_pvalues.tsv, hydra_all_counts.tsv

## Creating Row Means Normalized Celltype Dataframe
dropseq_final_atlas_analysis.R: results/processedData/genes_to_cell_type_table.tsf (rowmeans of each gene from each celltype)

## Assigning Unique Genes For Each Celltype
### Criteria
- Genes expressed over mean(row) for each Cell-Type
create_sc_atlas_cluster_table.ipynb: 
Input:
	../results/processedData/genes_to_cell_type_table.tsf

Output:
	../results/processedData/normalized_mean/: 
	
		hvaep_cell_type_to_gene_cluster_table.tsf, 
		hvaep_uniprot_kegg_go.tsf, 
		hvaep_cell_type_to_gene_cluster_table_just_gos.tsf,
		hvaep_uniprot_go_cleaned_table.tsf,
		hvaep_uniprot_go_full_table.tsf
	
	../results/deseq2_rsem/:
		hydra_all_counts_t_to_g.tsv
		hydra_all_deseq2_lg2_t_to_g.tsv
		hydra_all_deseq2_pvalues_t_to_g.tsv

## Plotting Overlaps
### Criteria
- Genes expressed in more than 1 Cell-Types are assigned to non-unique genes	
get_cluster.R:	
Input: 		
		./results/processedData/normalized_mean/hvaep_cell_type_to_gene_cluster_table.tsf
		./results/processedData/normalized_mean/hvaep_uniprot_kegg_go.tsf

Output:
		./results/processedData/normalized_mean/:
			hvaep_clustering_normalized.tsv
			hvaep_clustering_table.tsv
			hvaep_clustering_table_unique.tsv

		./results/figures/normalized_mean/:
			hvaep_clustering.png
			hvaep_GO_lineage_I_superclasses.png
			hvaep_GO_lineage_C_superclasses.png
			hvaep_GO_celltype.png
			hvaep_ribosomes.png
			hvaep_GO.png


## Clustering Celltypes
analyze_deseq2:

Input:
	./data/rsem_counts/sample_descriptions.csv

	./results/deseq2_rsem/hydra_all_counts_t_to_g.tsv
	./results/deseq2_rsem/hydra_all_deseq2_lg2_t_to_g.tsv
	./results/deseq2_rsem/hydra_all_deseq2_pvalues_t_to_g.tsv

	./results/processedData/normalized_mean/hvaep_clustering_table.tsv
	./results/processedData/normalized_mean/hvaep_uniprot_kegg_go.tsf

Output:
	./results/processedData/normalized_mean/:
		unique_celltypes.table
		ttest_all_statistic.table
		ttest_all_pvalue.table
		ttest_all_target_number.table
		ttest_all_query_number.table
		ttest_sorted_cut_statistic.table
		ttest_sorted_cut_pvalue.table
		ttest_sorted_cut_target_number.table
		ttest_sorted_cut_query_number.table
		kmeans_on_datz.table
		clt.table
		clusterings_gfzero.tsv
		
	./results/figures/normalized_mean/:
		celltype_ttest_updown_cut.png
		celltype_ttest_clustered.png
		celltype_ttest_clustered_both.png
		celltype_ttest_sorted.png
		celltype_wilcox.png
		celltype_wilcox_clustered.png
		cluster_ttest.png
		cormatrix_experiments.png
		cluster_celltype.png
		cluster_celltype_cut.png
		classes.png
		cluster_classes.png
		cluster_classes_sorted.png
		cluster_GO_1e-3.png
		cluster_GO_1e-4.png
		cluster_GO_1e-5.png
		cluster_GO_1e-6.png
		cluster_GO_1e-7.png
		cluster_GO_1e-8.png
		cluster_pca.png
		cluster_umap.png
		experiments_pca.png
		volcano_*.png

# Gene Enrichment Analysis
gene_enrichment_analysis_preparation.ipynb:

Input:
	../results/deseq2_rsem/tables/*
		EcoKD1_Eco1KD_B8_vs_control_B8_results.csv
		...
	../data/annotation_KO_GO.csv
	../results/processedData/normalized_mean/unique_celltypes.table (from analyze_deseq.R)

Output:
	../results/gene_enrichment/gene_to_KO.csv
	../results/gene_enrichment/*
		ecokd1/*
			ecokd1_C_EC_foot.csv
			...
		temp/*
		wild/*
		ahls/*
unknown data:
	- hvaepLRv2_kegg_go