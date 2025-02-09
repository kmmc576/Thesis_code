#FIGURES 4-7 and 4-8
#GSEA using Gprofiler2 - filter for genes up/down regulated in tetraploid relative to diploid by log2 ratio +/-1.5
#use ggvenn to create Venn diagrams for intersection of up/down regulated genes + GSEA terms
#export GSEA term sets for analysis in Cytoscape using Enrichment Map

library(gprofiler2)
library(tidyverse)
library(ggplot2)
library(ggvenn)
library(svglite)


setwd("/Users/.../R_files/")

#read in normalised cpm counts (cpm>1)
tetra_dip_cpm2 <- read.table("outputs/expression_data/tetra_dip_cpm.txt", header=TRUE, sep="\t")
#add in log2 columns
tetra_dip_cpm2 <- tetra_dip_cpm2 %>% mutate("log2_P14_dip"=log2(tetra_P14_cpm/diploid_cpm)) %>%
                                              mutate("log2_P30_dip"=log2(tetra_P30_cpm/diploid_cpm))

#check for duplicates
dup <- tetra_dip_cpm2$genes[duplicated(tetra_dip_cpm2$gene)]

#separate up/down-regulated genes, sort order
#choose up/down regulated threshold of log2 ratio >1.5, <1.5
#filter for lowly expressed genes

upgenes_P30 <- tetra_dip_cpm2 %>% filter(log2_P30_dip>1.5) %>% arrange(desc(log2_P30_dip), desc(tetra_P30_cpm)) %>%
  dplyr::select("chromosome","start","end", genes = "newgene", "tetra_P30_cpm","diploid_cpm", "log2_P30_dip") %>%
  filter(tetra_P30_cpm>5 | diploid_cpm >5) %>%
  arrange(desc(log2_P30_dip), desc(tetra_P30_cpm))

upgenes_P14 <- tetra_dip_cpm2 %>% filter(log2_P14_dip>1.5) %>% arrange(desc(log2_P14_dip), desc(tetra_P14_cpm)) %>%
  dplyr::select("chromosome","start","end", genes = "newgene", "tetra_P14_cpm","diploid_cpm", "log2_P14_dip") %>%
  filter(tetra_P14_cpm>5 | diploid_cpm >5) %>%
  arrange(desc(log2_P14_dip), desc(tetra_P14_cpm))

downgenes_P30 <- tetra_dip_cpm2 %>% filter(log2_P30_dip<(-1.5)) %>% arrange((log2_P30_dip), desc(diploid_cpm)) %>%
  dplyr::select("chromosome","start","end", genes = "newgene", "tetra_P30_cpm","diploid_cpm", "log2_P30_dip") %>%
  filter(tetra_P30_cpm>5 | diploid_cpm >5) %>%
  arrange(log2_P30_dip, desc(diploid_cpm))

downgenes_P14 <- tetra_dip_cpm2 %>% filter(log2_P14_dip<(-1.5)) %>% arrange((log2_P14_dip), desc(diploid_cpm)) %>%
  dplyr::select("chromosome","start","end", genes = "newgene", "tetra_P14_cpm","diploid_cpm", "log2_P14_dip") %>%
  filter(tetra_P14_cpm>5 | diploid_cpm >5) %>%
  arrange(log2_P14_dip, desc(diploid_cpm))


#Venn diagram of genes selected for analysis (+/- 1.5 log2)
ls <- list("Down P14\n(n=1576)" = downgenes_P14$genes, "Up P14\n(n=2373)" = upgenes_P14$genes, "Up P30\n(n=1616)" = upgenes_P30$genes, "Down P30\n(n=1169)" = downgenes_P30$genes)

venn_genesets <- ggvenn(ls, fill_color = c("#4DD0E1", "#F06292", "#FF8A65", "#64B5F6"), set_name_size = 2, text_size=4, show_percentage = FALSE)
venn_genesets

#list genes in each set of Venn
list_venn_genesets <- list_to_data_frame(ls)
colnames(list_venn_genesets) <- c("Gene", "P14D", "P14U", "P30U", "P30D")
write.table(list_venn_genesets, "/Users/kimmcintyre/Library/CloudStorage/OneDrive-UNSW/Kim_Thesis/Chapter_5_tetraploid/R_files/outputs/GSEA/list_venn_genesets.txt", sep = "\t", quote = F, col.names = TRUE, row.names = F)

list_venn_genesets <- read.table("/Users/kimmcintyre/Library/CloudStorage/OneDrive-UNSW/Kim_Thesis/Chapter_5_tetraploid/R_files/outputs/GSEA/list_venn_genesets.txt", header = TRUE, sep = "\t")

setwd("/Users/.../Supp_data4_GSEA/GSEA/")

#GSEA using Gprofiler2 with gSCS multiple testing correction
upgenes_enrich_P14 = gost(query = list("up" = upgenes_P14$genes), organism="hsapiens",
                          evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.001,
                          sources = c("GO:BP"))
P14U = upgenes_enrich_P14$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size")]
colnames(P14U) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size")
P14U$FDR = P14U$p.Val
P14U$Phenotype = "1"
P14U = P14U[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Gene.number", "Genes", "Term.size")]
P14U <- P14U %>% filter(Term.size<2000)

write.table(P14U, file = "gProfiler_P14U_g_SCS.txt", sep = "\t", quote = F, row.names = F)


upgenes_enrich_P30 = gost(query = list("up" = upgenes_P30$genes), organism="hsapiens",
                          evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.001,
                          sources = c("GO:BP"))
P30U = upgenes_enrich_P30$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size")]
colnames(P30U) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size")
P30U$FDR = P30U$p.Val
P30U$Phenotype = "1"
P30U = P30U[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Gene.number", "Genes", "Term.size")]
P30U <- P30U %>% filter(Term.size<2000)

write.table(P30U, file = "gProfiler_P30U_g_SCS.txt", sep = "\t", quote = F, row.names = F)


downgenes_enrich_P14 = gost(query = list("down" = downgenes_P14$genes), organism="hsapiens",
                            evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.001,
                            sources = c("GO:BP"))
P14D = downgenes_enrich_P14$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size")]
colnames(P14D) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size")
P14D$FDR = P14D$p.Val
P14D$Phenotype = "-1"
P14D = P14D[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Gene.number", "Genes", "Term.size")]
P14D <- P14D %>% filter(Term.size<2000)

write.table(P14D, file = "gProfiler_P14D_g_SCS.txt", sep = "\t", quote = F, row.names = F)


downgenes_enrich_P30 = gost(query = list("down" = downgenes_P30$genes), organism="hsapiens",
                            evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.001,
                            sources = c("GO:BP"))
P30D = downgenes_enrich_P30$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size")]
colnames(P30D) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size")
P30D$FDR = P30D$p.Val
P30D$Phenotype = "-1"
P30D = P30D[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Gene.number", "Genes", "Term.size")]
P30D <- P30D %>% filter(Term.size<2000)

write.table(P30D, file = "gProfiler_P30D_g_SCS.txt", sep = "\t", quote = F, row.names = F)


#Venn diagram of genes selected for analysis (+/- 1.5 log2)
ls <- list("Down P14\n(n=102)" = P14D$GO.ID, "Up P14\n(n=404)" = P14U$GO.ID, "Up P30\n(n=194)" = P30U$GO.ID, "Down P30\n(n=182)" = P30D$GO.ID)

venn_GSEAsets <- ggvenn(ls, fill_color = c("#A6E8F0", "#F7B1C9", "#FFC5B2", "#B2DAFB"), set_name_size = 2, text_size=4, show_percentage = FALSE)
venn_GSEAsets

#list GO terms in each set of Venn
list_venn_GSEAsets <- list_to_data_frame(ls)
colnames(list_venn_GSEAsets) <- c("GO.ID", "P14D", "P14U", "P30U", "P30D")
write.table(list_venn_GSEAsets, "/Users/.../R_files/outputs/GSEA/list_venn_GSEAsets.txt", sep = "\t", quote = F, col.names = TRUE, row.names = F)

summary_venn_GSEA <- list_venn_GSEAsets %>% count(P14D, P14U, P30U, P30D)
write.table(summary_venn_GSEA, "/Users/.../R_files/outputs/GSEA/summary_venn_GSEAsets.txt", sep = "\t", quote = F, col.names = TRUE, row.names = F)

#FIGURE 4-8B
#Barplots - most highly enriched GO terms
#filter for top 30, then manually curate for 16 (remove repetitive terms) 
#P14U
P14U_plot <- P14U %>% slice_min(FDR, n=50) %>% arrange(FDR)
P14U_plot <- P14U_plot %>% slice(1,2,5,12,13,14,20,21,22,30,33,36,38,47,50) %>% arrange(desc(FDR))

#use factor to fix order for plotting
P14U_plot$Description <- factor(P14U_plot$Description, levels = unique(P14U_plot$Description))

barplot_P14U <- ggplot(P14U_plot, aes(x= Description, y = -log10(FDR), width = 0.7))+
  geom_bar(aes(), fill = "#F26767", color = "#B8D0ED", stat = "identity")+
  labs(title="P14 upregulated")+
  ylab("")+
  xlab("")+
  scale_y_continuous(limits=c(0,45),
                     breaks=c( 20, 40),
                     labels=c("20","40"))+
  coord_flip()+
  theme_classic(base_size=20)+
  theme(legend.position = "none")
barplot_P14U
ggsave(file = "GSEA_barplot_P14U.svg", plot = barplot_P14U, path = "/outputs", width = 10, height = 7)
ggsave(file = "GSEA_barplot_P14U.png", plot = barplot_P14U, path = "/outputs", width = 10, height = 7)


#P14D
P14D_plot <- P14D %>% slice_min(FDR, n=50) %>% arrange(FDR)
P14D_plot <- P14D_plot %>% slice(1,2,3,7,8,11,14,15,18,19,20,22,24,26,39) %>% arrange(desc(FDR))

#use factor to fix order for plotting
P14D_plot$Description <- factor(P14D_plot$Description, levels = unique(P14D_plot$Description))


barplot_P14D <- ggplot(P14D_plot, aes(x= Description, y = -log10(FDR), width = 0.7))+
  geom_bar(aes(), fill = "#5690CC", color = "#B8D0ED", stat = "identity")+
  labs(title="P14 downregulated")+
  ylab("")+
  xlab("")+
  scale_y_continuous(limits=c(0,45),
                     breaks=c( 20, 40),
                     labels=c("20","40"))+
  coord_flip()+
  theme_classic(base_size=20)+
  theme(legend.position = "none")
barplot_P14D
ggsave(file = "GSEA_barplot_P14D.svg", plot = barplot_P14D, path = "/outputs/GSEA", width = 10, height = 7)
ggsave(file = "GSEA_barplot_P14D.png", plot = barplot_P14D, path = "/outputs/GSEA", width = 10, height = 7)

#P30U
P30U_plot <- P30U %>% slice_min(FDR, n=50) %>% arrange(FDR)
P30U_plot <- P30U_plot %>% arrange(FDR) %>% slice(1,2,6,7,11,12,13,18,21,26,27,30,39,43,44) %>% arrange(desc(FDR))

#use factor to fix order for plotting
P30U_plot$Description <- factor(P30U_plot$Description, levels = unique(P30U_plot$Description))


barplot_P30U <- ggplot(P30U_plot, aes(x= Description, y = -log10(FDR), width = 0.7))+
  geom_bar(aes(), fill = "#FFCD34", color = "#B8D0ED", stat = "identity")+
  labs(title="P30 upregulated")+
  ylab("")+
  xlab("")+
  scale_y_continuous(limits=c(0,45),
                     breaks=c( 20, 40),
                     labels=c("20","40"))+
  coord_flip()+
  theme_classic(base_size=20)+
  theme(legend.position = "none")
barplot_P30U
ggsave(file = "GSEA_barplot_P30U.svg", plot = barplot_P30U, path = "/outputs/GSEA", width = 10, height = 7)
ggsave(file = "GSEA_barplot_P30U.png", plot = barplot_P30U, path = "/outputs/GSEA", width = 10, height = 7)


#P30D
P30D_plot <- P30D %>% slice_min(FDR, n=50) %>% arrange(FDR)
P30D_plot <- P30D_plot %>% arrange(FDR) %>% slice(1,2,4,6,8,9,12,13,17,18,23,26,28,41,43) %>% arrange(desc(FDR))

#use factor to fix order for plotting
P30D_plot$Description <- factor(P30D_plot$Description, levels = unique(P30D_plot$Description))


barplot_P30D <- ggplot(P30D_plot, aes(x= Description, y = -log10(FDR), width = 0.7))+
  geom_bar(aes(), fill = "#41B649", color = "#B8D0ED", stat = "identity")+
  labs(title="P30 downregulated")+
  ylab("")+
  xlab("")+
  scale_y_continuous(limits=c(0,45),
                     breaks=c( 20, 40),
                     labels=c("20","40"))+
  coord_flip()+
  theme_classic(base_size=20)+
  theme(legend.position = "none")
barplot_P30D
ggsave(file = "GSEA_barplot_P30D.svg", plot = barplot_P30D, path = "/outputs/GSEA", width = 10, height = 7)
ggsave(file = "GSEA_barplot_P30D.png", plot = barplot_P30D, path = "/outputs/GSEA", width = 10, height = 7)

#combine into single tibble for facet_wrap
P14U2 <- P14U_plot %>% select(Description, FDR) %>% mutate(plot = "P14 upregulated")
P14D2 <- P14D_plot %>% select(Description, FDR) %>% mutate(plot = "P14 downregulated")
P30U2 <- P30U_plot %>% select(Description, FDR) %>% mutate(plot = "P30 upregulated")
P30D2 <- P30D_plot %>% select(Description, FDR) %>% mutate(plot = "P30 downregulated")
all <- rbind(P14U2,P14D2,P30U2,P30D2)  

barplot_all <- ggplot(all, aes(x= Description, y = -log10(FDR), width = 0.7))+
  geom_bar(aes(), fill = "#41B649", color = "#B8D0ED", stat = "identity")+
  ylab("")+
  xlab("")+
  scale_y_continuous(limits=c(0,45),
                     breaks=c( 20, 40),
                     labels=c("20","40"))+
  coord_flip()+
  theme_classic(base_size=20)+
  theme(legend.position = "none")+
  facet_wrap(~plot)
barplot_all
  
  
