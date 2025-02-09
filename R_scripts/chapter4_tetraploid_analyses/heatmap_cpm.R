#FIGURE 4-6C
#This to create a heatmap of all expression data >1 cpm in at least one sample (~ > 0.05 rpkm) + GSEA of clusters

library(tidyverse)
library(pheatmap)
library(ggplotify)
library(patchwork)
library(gprofiler2)
library(BioVenn)

setwd("/Users/.../R_files/")

#read in file with read counts (cpm)
tetra_dip_cpm <- read.delim("outputs/expression_data/tetra_dip_cpm.txt", header=TRUE, sep="\t")

#create data frame 
tetra_dip_cpm_heatmap <- data.frame("gene" = tetra_dip_cpm$genes, "diploid" = tetra_dip_cpm$diploid_cpm, "tetra_P14" = tetra_dip_cpm$tetra_P14_cpm, "tetra_P30" = tetra_dip_cpm$tetra_P30_cpm)

tetra_dip_cpm_heatmap <- distinct(tetra_dip_cpm_heatmap, gene, .keep_all = TRUE)

#calculate Z-score, row-wise
tetra_dip_cpm_heatmap <- tetra_dip_cpm_heatmap %>% rowwise() %>%
  mutate(mn = mean(c(diploid, tetra_P14, tetra_P30))) %>%
  mutate(stdev = sd(c(diploid, tetra_P14, tetra_P30))) %>%
  mutate(dip_zsc = c((diploid - mn)/stdev)) %>%
  mutate(P14_zsc = c((tetra_P14 - mn)/stdev)) %>%
  mutate(P30_zsc = c((tetra_P30 - mn)/stdev))

#create matrix for plotting, delete NaN rows 
tetra_dip_cpm_heatmap <- dplyr::select(tetra_dip_cpm_heatmap, gene = "gene", diploid = "dip_zsc", tetra_P14 = "P14_zsc", tetra_P30 = "P30_zsc") 

tetra_dip_cpm_heatmap <- as.data.frame(tetra_dip_cpm_heatmap)
.rowNamesDF(tetra_dip_cpm_heatmap, make.names=TRUE) <- tetra_dip_cpm_heatmap[,1]
tetra_dip_cpm_heatmap <- tetra_dip_cpm_heatmap[,2:4]
tetra_dip_cpm_heatmap <- as.numeric(tetra_dip_cpm_heatmap)
tetra_dip_cpm_heatmap <- as.matrix(tetra_dip_cpm_heatmap)

#omit "na" values
tetra_dip_cpm_heatmap <- na.omit(tetra_dip_cpm_heatmap)

#create heatmap
#use kmeans to cluster before constructing heatmap - partions into groups based on minimisation of sum of squares
pheatmap_tetra_dip_cpm_heatmap <- pheatmap(tetra_dip_cpm_heatmap, 
                                      cellwidth = 50,
                                       cluster_cols = F,
                                      cutree_rows = 5,
                                      show_rownames = F)
pheatmap_tetra_dip_cpm_heatmap

ggsave("outputs/heatmap/heatmap_cpm.png", plot = pheatmap_tetra_dip_cpm_heatmap, device = "png", width = 10, height = 10)
#save manually to dictate proportions
#create dataframe of gene names by cluster
cluster = cutree(pheatmap_tetra_dip_cpm_heatmap$tree_row, k = 5)
cluster2 <- as.data.frame(cluster)
cluster3 <- tibble::rownames_to_column(cluster2, var = "gene")
cluster3 <- as_tibble(cluster3)

#convert heatmap matrix to tibble, then merge with cluster numbers
mtrx <- as.data.frame(tetra_dip_cpm_heatmap)
mtrx2 <- tibble::rownames_to_column(mtrx, var = "gene")
mtrx2 <- as_tibble(mtrx2)

cluster4 <- left_join(mtrx2, cluster3)

#GSEA by cluster using gprofiler2
cl1 <- cluster4 %>% filter(cluster == 1)
cl2 <- cluster4 %>% filter(cluster == 2)
cl3 <- cluster4 %>% filter(cluster == 3)
cl4 <- cluster4 %>% filter(cluster == 4)
cl5 <- cluster4 %>% filter(cluster == 5)

cl1_enrich <- gost(query = list("up"= cl1$gene), organism="mdomestica",
                     ordered_query = FALSE, correction_method = "fdr", user_threshold = 0.001, sources = c("GO:BP", "GO:MF", "GO:CC"))

cl2_enrich <- gost(query = list("up"= cl2$gene), organism="mdomestica",
                   ordered_query = FALSE, correction_method = "fdr", user_threshold = 0.001, sources = c("GO:BP", "GO:MF", "GO:CC"))

cl3_enrich <- gost(query = list("up"= cl3$gene), organism="mdomestica",
                   ordered_query = FALSE, correction_method = "fdr", user_threshold = 0.001, sources = c("GO:BP", "GO:MF", "GO:CC"))

cl4_enrich <- gost(query = list("up"= cl4$gene), organism="mdomestica",
                   ordered_query = FALSE, correction_method = "fdr", user_threshold = 0.001, sources = c("GO:BP", "GO:MF", "GO:CC"))

cl5_enrich <- gost(query = list("up"= cl5$gene), organism="mdomestica",
                   ordered_query = FALSE, correction_method = "fdr", user_threshold = 0.001, sources = c("GO:BP", "GO:MF", "GO:CC"))

                 
cl1_enrich_result <- cl1_enrich$result %>% dplyr::select(-query, -significant) %>% filter(term_size < 3000) %>% arrange(desc(p_value))
cl2_enrich_result <- cl2_enrich$result %>% dplyr::select(-query, -significant) %>% filter(term_size < 3000) %>% arrange(desc(p_value))
cl3_enrich_result <- cl3_enrich$result %>% dplyr::select(-query, -significant) %>% filter(term_size < 3000) %>% arrange(desc(p_value))
cl4_enrich_result <- cl4_enrich$result %>% dplyr::select(-query, -significant) %>% filter(term_size < 3000) %>% arrange(desc(p_value))
cl5_enrich_result <- cl5_enrich$result %>% dplyr::select(-query, -significant) %>% filter(term_size < 3000) %>% arrange(desc(p_value))

cl1_enrich_GO <- tibble("GO_term" = cl1_enrich_result$term_name, "GO_id" = cl1_enrich_result$term_id, "source" = cl1_enrich_result$source, "term_size" = cl1_enrich_result$term_size, 
                          "gene_number" = cl1_enrich_result$intersection_size, "p_value"= cl1_enrich_result$p_value) %>% arrange(p_value)
cl2_enrich_GO <- tibble("GO_term" = cl2_enrich_result$term_name, "GO_id" = cl2_enrich_result$term_id, "source" = cl2_enrich_result$source,"term_size" = cl2_enrich_result$term_size, 
                        "gene_number" = cl2_enrich_result$intersection_size, "p_value"= cl2_enrich_result$p_value) %>% arrange(p_value)
cl3_enrich_GO <- tibble("GO_term" = cl3_enrich_result$term_name, "GO_id" = cl3_enrich_result$term_id, "source" = cl3_enrich_result$source,"term_size" = cl3_enrich_result$term_size, 
                        "gene_number" = cl3_enrich_result$intersection_size, "p_value"= cl3_enrich_result$p_value) %>% arrange(p_value)
cl4_enrich_GO <- tibble("GO_term" = cl4_enrich_result$term_name, "GO_id" = cl4_enrich_result$term_id, "source" = cl4_enrich_result$source,"term_size" = cl4_enrich_result$term_size, 
                        "gene_number" = cl4_enrich_result$intersection_size, "p_value"= cl4_enrich_result$p_value) %>% arrange(p_value)
cl5_enrich_GO <- tibble("GO_term" = cl5_enrich_result$term_name, "GO_id" = cl5_enrich_result$term_id, "source" = cl5_enrich_result$source,"term_size" = cl5_enrich_result$term_size, 
                        "gene_number" = cl5_enrich_result$intersection_size, "p_value"= cl5_enrich_result$p_value) %>% arrange(p_value)

setwd("/Users/.../R_files/outputs/heatmap")
write.table(cl1, file = "heatmap_cluster1_cpm.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(cl2, file = "heatmap_cluster2_cpm.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(cl3, file = "heatmap_cluster3_cpm.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(cl4, file = "heatmap_cluster4_cpm.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(cl5, file = "heatmap_cluster5_cpm.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(cl1_enrich_GO, file = "GSEA_heatmap_cluster1_cpm.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(cl2_enrich_GO, file = "GSEA_heatmap_cluster2_cpm.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(cl3_enrich_GO, file = "GSEA_heatmap_cluster3_cpm.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(cl4_enrich_GO, file = "GSEA_heatmap_cluster4_cpm.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(cl5_enrich_GO, file = "GSEA_heatmap_cluster5_cpm.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

