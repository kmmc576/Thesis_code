#FIGURE 4-8A
##this to do Gprofiler enrichment searches for up and down regulated genes.
#generates files to import into cytoscape using Enrichment Map
#run Gprofiler using human annotations (more complete than MDO)
#use column "newgene" which contains the updated HGNC gene name, where it differs from the MonDom5 name


library(gprofiler2)
library(tidyverse)

setwd("/Users/.../R_files/")

#read in normalised cpm counts (cpm>1):
tetra_dip_cpm2 <- read.table("outputs/expression_data/tetra_dip_cpm.txt", header=TRUE, sep="\t")

#separate up/down-regulated genes, sort order
#choose up/down regulated threshold of log2 ratio >1.5, <-1.5
#need to filter for lowly expressed genes, so gene order not biased by low expression genes
upgenes_P30 <- tetra_dip_cpm2 %>% filter(log2_P30_dip>1.5) %>% arrange(desc(log2_P30_dip), desc(tetra_P30_cpm)) %>%
  dplyr::select("chromosome","start","end", "newgene", "tetra_P30_cpm","diploid_cpm", "log2_P30_dip") %>%
  filter(tetra_P30_cpm>5 | diploid_cpm >5)


upgenes_P14 <- tetra_dip_cpm2 %>% filter(log2_P14_dip>1.5) %>% arrange(desc(log2_P14_dip), desc(tetra_P14_cpm)) %>%
  dplyr::select("chromosome","start","end", "newgene", "tetra_P14_cpm","diploid_cpm", "log2_P14_dip") %>%
  filter(tetra_P14_cpm>5 | diploid_cpm >5)


downgenes_P30 <- tetra_dip_cpm2 %>% filter(log2_P30_dip<(-1.5)) %>% arrange((log2_P30_dip), desc(diploid_cpm)) %>%
  dplyr::select("chromosome","start","end", "newgene", "tetra_P30_cpm","diploid_cpm", "log2_P30_dip") %>%
  filter(tetra_P30_cpm>5 | diploid_cpm >5)


downgenes_P14 <- tetra_dip_cpm2 %>% filter(log2_P14_dip<(-1.5)) %>% arrange((log2_P14_dip), desc(diploid_cpm)) %>%
  dplyr::select("chromosome","start","end", "newgene", "tetra_P14_cpm","diploid_cpm", "log2_P14_dip") %>%
  filter(tetra_P14_cpm>5 | diploid_cpm >5)


#RUN using gSCS multiple testing correction
upgenes_enrich_P14 = gost(query = list("up" = upgenes_P14$newgene), organism="hsapiens",
                          evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.01,
                          sources = c("GO:BP"))
P14U = upgenes_enrich_P14$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size")]
colnames(P14U) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size")
P14U$FDR = P14U$p.Val
P14U$Phenotype = "1"
P14U = P14U[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Gene.number", "Genes", "Term.size")]
P14U <- P14U %>% filter(Term.size<2000)
# saving the P14U file
write.table(P14U, file = "outputs/GSEA_for_cytoscape/gProfiler_P14U_g_SCS.txt", sep = "\t", quote = F, row.names = F)


upgenes_enrich_P30 = gost(query = list("up" = upgenes_P30$newgene), organism="hsapiens",
                          evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.01,
                          sources = c("GO:BP"))
P30U = upgenes_enrich_P30$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size")]
colnames(P30U) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size")
P30U$FDR = P30U$p.Val
P30U$Phenotype = "1"
P30U = P30U[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Gene.number", "Genes", "Term.size")]
P30U <- P30U %>% filter(Term.size<2000)
# saving the P30U file
write.table(P30U, file = "outputs/GSEA_for_cytoscape/gProfiler_P30U_g_SCS.txt", sep = "\t", quote = F, row.names = F)


downgenes_enrich_P14 = gost(query = list("down" = downgenes_P14$newgene), organism="hsapiens",
                            evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.01,
                            sources = c("GO:BP"))
P14D = downgenes_enrich_P14$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size")]
colnames(P14D) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size")
P14D$FDR = P14D$p.Val
P14D$Phenotype = "-1"
P14D = P14D[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Gene.number", "Genes", "Term.size")]
P14D <- P14D %>% filter(Term.size<2000)
# saving the P14D file
write.table(P14D, file = "outputs/GSEA_for_cytoscape/gProfiler_P14D_g_SCS.txt", sep = "\t", quote = F, row.names = F)

downgenes_enrich_P30 = gost(query = list("down" = downgenes_P30$newgene), organism="hsapiens",
                            evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.01,
                            sources = c("GO:BP"))
P30D = downgenes_enrich_P30$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size")]
colnames(P30D) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size")
P30D$FDR = P30D$p.Val
P30D$Phenotype = "-1"
P30D = P30D[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Gene.number", "Genes", "Term.size")]
P30D <- P30D %>% filter(Term.size<2000)
# saving the P30D file
write.table(P30D, file = "outputs/GSEA_for_cytoscape/gProfiler_P30D_g_SCS.txt", sep = "\t", quote = F, row.names = F)






