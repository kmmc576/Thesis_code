
library(tidyverse)
library(ggplot2)
library(fabricatr)
library(ggplotify)
library(patchwork)
library(rstatix)
library(ggpubr)
library(ggvenn)
library(rcompanion)
library(splitstackshape)
library(glue)
library(stats)
library(FSA) #Dunn's test
library(readr)
library(ggplot2)
library(RVAideMemoire)
library(coin)
#read in normalised cpm countdata from R script RNA_expression_data_cpm.R
setwd("/Users/.../R_files/")

setwd("/Users/kimmcintyre/Library/CloudStorage/OneDrive-UNSW/Kim_Thesis/Chapter_5_tetraploid/R_files/")

tetra_dip_cpm <- read.table("outputs/expression_data/tetra_dip_cpm.txt", header=TRUE, sep="\t")

#subset genes based on each sample separately
#use split_quantile to assign factor numbers to each row based on order of diploid, tetraploid cpm, then add as column

#add quartiles to dataframe, order by chromomsome then TSS
dip_quartile <- tetra_dip_cpm %>% dplyr::select("cpm" = diploid_cpm, genes, chromosome) %>% filter(cpm > 1) %>%
  mutate("quartile" = split_quantile(x=log2(cpm), type=4)) %>%
  mutate("log2_cpm" = log2(cpm), .before = quartile) %>% mutate("cat" = "diploid") 

P14_quartile <- tetra_dip_cpm %>% dplyr::select("cpm" = tetra_P14_cpm, genes, chromosome) %>% filter(cpm > 1) %>%
  mutate("quartile" = split_quantile(x=log2(cpm), type=4)) %>%
  mutate("log2_cpm" = log2(cpm), .before = quartile) %>% mutate("cat" = "tet_P14") 

P30_quartile <- tetra_dip_cpm %>% dplyr::select("cpm" = tetra_P30_cpm, genes, chromosome) %>% filter(cpm > 1) %>%
  mutate("quartile" = split_quantile(x=log2(cpm), type=4)) %>%
  mutate("log2_cpm" = log2(cpm), .before = quartile) %>% mutate("cat" = "tet_P30") 


all_quart <- rbind(dip_quartile, P14_quartile, P30_quartile) %>% filter(log2_cpm != "-Inf")

#count unique gene name
genecount <- all_quart %>% count(genes)
#13406

#filter out -inf
all_q1 <- all_quart %>% filter(quartile == 1) 
all_q2 <- all_quart %>% filter(quartile == 2) 
all_q3 <- all_quart %>% filter(quartile == 3) 
all_q4 <- all_quart %>% filter(quartile == 4) 

#Define colours, so consistent for each sample:
#Diploid: #F8766D 
#Tetra_P14  #32D824 
#Tetra_P30 #619CFF 

#FIGURE 4-9
#To compare output by chromosome across all samples
setwd("/Users/.../R_files/outputs/expression_by_chromosome_cpm/")

chr = unique(all_quart$chromosome)

chr_median_plots = list()

for(j in chr){
  KW <- rstatix::kruskal_test(data = all_quart %>% filter(chromosome == j), log2_cpm ~ cat)
  KW$p <- format(KW$p, scientific = TRUE, digits = 3)
  
  DT <- dunn_test(all_quart %>% filter(chromosome == j), log2_cpm ~ cat, p.adjust.method = "holm")
  DT$p.adj <- format(DT$p.adj, scientific = TRUE, digits = 3)
  DT <- DT %>% mutate(p.adj2 = (if_else(p.adj.signif == "ns", "", p.adj)))
    
  chr_median_plots[[j]] = ggboxplot(all_quart %>% filter(chromosome == j), x = "cat", y = "log2_cpm", color = "black",
                                      fill = "cat",
                                      palette = c("coral", "chartreuse2", "cornflowerblue"), width = 0.4, alpha = 0.8,
                                      title = paste0("Chr", j),
                                      subtitle = paste0("Kruskal-Wallis p value =", KW$p),
                                      xlab = "",
                                      ylab = "Log2 cpm",
                                      notch = TRUE) +
    theme(legend.position = "none")+  
    scale_x_discrete(breaks = c("diploid", "tet_P14", "tet_P30"), labels = c("Dip", "TetP14", "TetP30"))+
      stat_summary(fun=median, geom = "text", show.legend = FALSE, vjust=-1, size = 3.5, aes( label=round(..y.., digits=2)))+
      stat_pvalue_manual(DT, label = "{p.adj.signif} {p.adj2}", y.position = c(11,12.5,14))+
      geom_text(stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=0.5, size=3.5)+
    print(chr_median_plots[[j]])
    ggsave(chr_median_plots[[j]], file=paste0("Chr", j, "_all_median_KW.svg"))
  }

#chr_median_plots is list, number inside [[]] denotes list item

p1 <- as.ggplot(chr_median_plots[[5]])
p2 <- as.ggplot(chr_median_plots[[8]])
p3 <- as.ggplot(chr_median_plots[[2]])
p4 <- as.ggplot(chr_median_plots[[6]])
p5 <- as.ggplot(chr_median_plots[[3]])
p6 <- as.ggplot(chr_median_plots[[4]])
p7 <-as.ggplot(chr_median_plots[[7]])
p8 <- as.ggplot(chr_median_plots[[1]])
pX <- as.ggplot(chr_median_plots[[10]])
pzzz <- as.ggplot(chr_median_plots[[9]])

# use patchwork to arrange them together

chrs <- p5 + p4 + p8 + p6 + p7 + p3 + p1 + p2 + pX + plot_layout(ncol = 3)
chrs
ggsave(file = "chr_all_median_KW_Dunn_cpm.svg", plot = chrs,  device = "svg", width = 15, height = 15, units = "in")

#FIGURE 4-9B - include only chromosome 2
#chromosomes - by quartiles
#difficult to position stat test automatically, so do loop for chromosomes and define fixed locations for each quartile 
setwd("/Users/.../R_files/outputs/expression_by_chromosome_quartiles_cpm/")

chr = unique(all_quart$chromosome)
quart = unique(all_quart$quartile)

chr_median_plots_quartiles = list()

for(j in chr){
for(i in quart){
  KW <- rstatix::kruskal_test(data = all_quart %>% filter(chromosome == j & quartile == i), log2_cpm ~ cat)
  KW$p <- format(KW$p, scientific = TRUE, digits = 3)
  
  DT <- dunn_test(all_quart %>% filter(chromosome == j & quartile == i), log2_cpm ~ cat, p.adjust.method = "holm")
  DT$p.adj <- format(DT$p.adj, scientific = TRUE, digits = 3)
  DT <- DT %>% mutate(p.adj2 = (if_else(p.adj.signif == "ns", "", p.adj)))
  
  chr_median_plots_quartiles[[j]][[i]] = ggboxplot(all_quart %>% filter(chromosome == j & quartile == i), x = "cat", y = "log2_cpm", color = "black",
                                    fill = "cat",
                                    palette = c("coral", "chartreuse2", "cornflowerblue"), width = 0.4, alpha = 0.8,
                                    title = paste0("Chr", j, " Q",i),
                                    subtitle = paste0("Kruskal-Wallis p value =", KW$p),
                                    xlab = "",
                                    ylab = "Log2 cpm",
                                    notch = TRUE) +
    scale_x_discrete(breaks = c("diploid", "tet_P14", "tet_P30"), labels = c("Dip", "TetP14", "TetP30"))+
    stat_summary(fun=median, geom = "text", show.legend = FALSE, vjust=-1, size = 3.5, aes( label=round(..y.., digits=2)))+
    stat_pvalue_manual(DT, label = "{p.adj.signif} {p.adj2} (Dunn)", y.position = c(11,12.5,14))+
    geom_text(stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=0.5, size=3.5)+
  print(chr_median_plots_quartiles[[j]][[i]])
  ggsave(chr_median_plots_quartiles[[j]][[i]], file=paste0("Chr", j, "_Q", i, "_median_KW_Dunn.tiff"))
  plot_name <- paste0("p",j,i)
  assign(plot_name, chr_median_plots_quartiles[[j]][[i]])
}
  }

KW2 <- rstatix::kruskal_test(data = all_quart %>% filter(chromosome == 2 & quartile == 3), log2_cpm ~ cat)
# assign plot name names each plot in R, so can use with patchwork. 
# use patchwork to arrange them together
p11 <- as.ggplot(p11)
p12 <- as.ggplot(p12)
p13 <- as.ggplot(p13)
p14 <- as.ggplot(p14)

chr1_quartiles <- p11 + p12 + p13 + p14 + plot_layout(ncol = 4)
chr1_quartiles
ggsave("Chr1_quartiles_median.tiff", plot = chr1_quartiles, device = "tiff", width = 35, height = 20, units = "cm")

p21 <- as.ggplot(p21)
p22 <- as.ggplot(p22)
p23 <- as.ggplot(p23)
p24 <- as.ggplot(p24)

chr2_quartiles <- p21 + p22 + p23 + p24 + plot_layout(ncol = 4)
chr2_quartiles
ggsave("Ch2_quartiles_median.tiff", plot = chr2_quartiles, device = "tiff", width = 35, height = 20, units = "cm")


p31 <- as.ggplot(p31)
p32 <- as.ggplot(p32)
p33 <- as.ggplot(p33)
p34 <- as.ggplot(p34)

chr3_quartiles <- p31 + p32 + p33 + p34 + plot_layout(ncol = 4)
chr3_quartiles
ggsave("Ch3_quartiles_median.tiff", plot = chr3_quartiles, device = "tiff", width = 35, height = 20, units = "cm")

p41 <- as.ggplot(p41)
p42 <- as.ggplot(p42)
p43 <- as.ggplot(p43)
p44 <- as.ggplot(p44)

chr4_quartiles <- p41 + p42 + p43 + p44 + plot_layout(ncol = 4)
chr4_quartiles
ggsave("Ch4_quartiles_median.tiff", plot = chr4_quartiles, device = "tiff", width = 35, height = 20, units = "cm")

p51 <- as.ggplot(p51)
p52 <- as.ggplot(p52)
p53 <- as.ggplot(p53)
p54 <- as.ggplot(p54)

chr5_quartiles <- p51 + p52 + p53 + p54 + plot_layout(ncol = 4)
chr5_quartiles
ggsave("Ch5_quartiles_median.tiff", plot = chr5_quartiles, device = "tiff", width = 35, height = 20, units = "cm")


p61 <- as.ggplot(p61)
p62 <- as.ggplot(p62)
p63 <- as.ggplot(p63)
p64 <- as.ggplot(p64)

chr6_quartiles <- p61 + p62 + p63 + p64 + plot_layout(ncol = 4)
chr6_quartiles
ggsave("Ch6_quartiles_median.tiff", plot = chr6_quartiles, device = "tiff", width = 35, height = 20, units = "cm")

p71 <- as.ggplot(p71)
p72 <- as.ggplot(p72)
p73 <- as.ggplot(p73)
p74 <- as.ggplot(p74)

chr7_quartiles <- p71 + p72 + p73 + p74 + plot_layout(ncol = 4)
chr7_quartiles
ggsave("Ch7_quartiles_median.tiff", plot = chr7_quartiles, device = "tiff", width = 35, height = 20, units = "cm")

p81 <- as.ggplot(p81)
p82 <- as.ggplot(p82)
p83 <- as.ggplot(p83)
p84 <- as.ggplot(p84)

chr8_quartiles <- p81 + p82 + p83 + p84 + plot_layout(ncol = 4)
chr8_quartiles
ggsave("Ch8_quartiles_median.tiff", plot = chr8_quartiles, device = "tiff", width = 35, height = 20, units = "cm")


pX1 <- as.ggplot(pX1)
pX2 <- as.ggplot(pX2)
pX3 <- as.ggplot(pX3)
pX4 <- as.ggplot(pX4)

chrX_quartiles <- pX1 + pX2 + pX3 + pX4 + plot_layout(ncol = 4)
chrX_quartiles
ggsave("ChX_quartiles_median.tiff", plot = chrX_quartiles, device = "tiff", width = 35, height = 20, units = "cm")


pzzz1 <- as.ggplot(pzzz1)
pzzz2 <- as.ggplot(pzzz2)
pzzz3 <- as.ggplot(pzzz3)
pzzz4 <- as.ggplot(pzzz4)

chrUnanch_quartiles <- pzzz1 + pzzz2 + pzzz3 + pzzz4 + plot_layout(ncol = 4)
chrUnanch_quartiles
ggsave("ChUnanch_quartiles_median.tiff", plot = chrUnanch_quartiles, device = "tiff", width = 35, height = 20, units = "cm")

#FIGURE 4-9B
#AMEND format for just chromosome 2 quartiles:
setwd("/Users/.../R_files/outputs/expression_by_chromosome_quartiles_cpm/")
quart = unique(all_quart$quartile)
chr_median_plots_quartiles_chr2 = list()
  for(i in quart){
    KW <- rstatix::kruskal_test(data = all_quart %>% filter(chromosome == 2 & quartile == i), log2_cpm ~ cat)
    KW$p <- format(KW$p, scientific = TRUE, digits = 3)
    
    DT <- dunn_test(all_quart %>% filter(chromosome == 2 & quartile == i), log2_cpm ~ cat, p.adjust.method = "holm")
    DT$p.adj <- format(DT$p.adj, scientific = TRUE, digits = 3)
    DT <- DT %>% mutate(p.adj2 = (if_else(p.adj.signif == "ns", "", p.adj)))
    
    chr_median_plots_quartiles_chr2[[i]] = ggboxplot(all_quart %>% filter(chromosome == 2 & quartile == i), x = "cat", y = "log2_cpm", color = "black",
                                                     fill = "cat",
                                                     palette = c("coral", "chartreuse2", "cornflowerblue"), width = 0.3, alpha = 0.8,
                                                     title = paste0("Chr 2", " Q",i),
                                                     subtitle = paste0("Kruskal-Wallis p-value =", KW$p),
                                                     xlab = "",
                                                     ylab = "Log2 cpm",
                                                     ylim = c(0,16),
                                                     notch = TRUE) +
      theme(legend.position = "none")+ 
      scale_x_discrete(breaks = c("diploid", "tet_P14", "tet_P30"), labels = c("Dip", "TetP14", "TetP30"))+
      stat_summary(fun=median, geom = "text", show.legend = FALSE, vjust=-2.5, hjust = -0.1, size = 3.5, aes( label=round(..y.., digits=2)))+
      stat_pvalue_manual(DT, label = "{p.adj.signif} {p.adj2}", y.position = c(12,13,14))+
      geom_text(stat="count", aes(label=..count..), y=-Inf, vjust=-0.1, hjust=0.5, size=3.5)
      print(chr_median_plots_quartiles_chr2[[i]])
    ggsave(chr_median_plots_quartiles_chr2[[i]] , file=paste0("Chr 2_Q", i, "_median_KW_Dunn_amended.svg"), width = 10, height = 20, units = "cm")
  }

#FIGURE 4-19A
#To plot box plots of all autosomal to cf X
setwd("/Users/.../R_files/outputs/expression_by_chromosome_cpm/")

  KW <- rstatix::kruskal_test(data = all_quart %>% filter(chromosome != "zzz" & chromosome != "X"), log2_cpm ~ cat)
  KW$p <- format(KW$p, scientific = TRUE, digits = 3)
  
  DT <- dunn_test(all_quart %>% filter(chromosome != "zzz" & chromosome != "X"), log2_cpm ~ cat, p.adjust.method = "holm")
  DT$p.adj <- format(DT$p.adj, scientific = TRUE, digits = 3)
  DT <- DT %>% mutate(p.adj2 = (if_else(p.adj.signif == "ns", "", p.adj)))
  
  chr_median_plots_auto = ggboxplot(all_quart %>% filter(chromosome != "zzz" & chromosome != "X"), x = "cat", y = "log2_cpm", color = "black",
                                    fill = "cat",
                                    palette = c("coral", "chartreuse2", "cornflowerblue"), width = 0.4, alpha = 0.8,
                                    title = "Autosomes",
                                    subtitle = paste0("Kruskal-Wallis p=", KW$p),
                                    xlab = "",
                                    ylab = "Log2 cpm",
                                    notch = TRUE) +
    theme(legend.position = "none")+  
    scale_x_discrete(breaks = c("diploid", "tet_P14", "tet_P30"), labels = c("Dip", "TetP14", "TetP30"))+
    coord_cartesian(ylim = c(0,15))+
    stat_summary(fun=median, geom = "text", show.legend = FALSE, vjust=-0.5, size = 3.5, aes( label=round(..y.., digits=2)))+
    stat_pvalue_manual(DT, label = "{p.adj.signif} {p.adj}", y.position = c(11,12.5,14))+
    geom_text(stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=0.5, size=3.5)
  
  chr_median_plots_auto
  ggsave("Autosomal_all_median_KW.svg", plot = chr_median_plots_auto, device = "svg", width = 13, height = 15, units = "cm")
  

#FIGURE 4-19
#To plot X:autosomal expression
#Stat test: compare median log2_cpm - use Mood's median test (null hypothesis: medians of two populations are identical)
  #+ Mann Whitney U test (wilcoxon unpaired) (as alternative) - tests similarity of distribution, not just medians


setwd("/Users/.../R_files/outputs/expression_by_chromosome_cpm/")

all_quart2 <- all_quart %>% filter(chromosome != "zzz") %>% mutate("type" = if_else(chromosome == "X", "X", "A"))

cat = unique(all_quart2$cat)

XA_median_plots = list()


for(j in cat){
  MO <- mood.medtest(data = all_quart2 %>% filter(cat == j), log2_cpm ~ type)
  MO$p.value <- format(MO$p.value, scientific = FALSE, digits = 2)
  
  MW <- wilcox.test(data = all_quart2 %>% filter(cat == j), log2_cpm ~ type, alternative = "two.sided")
  MW$p.value <- format(MW$p.value, scientific = FALSE, digits = 2)

  XA_median_plots[[j]] = ggboxplot(all_quart2 %>% filter(cat == j), x = "type", y = "log2_cpm", color = "black",
                                   fill = "type",
                                   palette = c("lightblue3", "orchid2"), width = 0.4, alpha = 0.8,
                                   title = paste0(j),
                                   subtitle = paste0("Mood's median p=", MO$p.value, ", Mann WhitU p=", MW$p.value),
                                   xlab = "",
                                   ylab = "Log2 cpm",
                                   notch = TRUE) +
    theme(legend.position = "none")+ 
    coord_cartesian(ylim = c(0,15))+
    scale_x_discrete(breaks = c("A", "X"), labels = c("A", "X"))+
    stat_summary(fun=median, geom = "text", show.legend = FALSE, vjust=-1, size = 3.5, aes( label=round(..y.., digits=2)))+
    geom_text(stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=0.5, size=3.5)+
    print(XA_median_plots[[j]])
  ggsave(XA_median_plots[[j]], file=paste0(j, "_median_XA.svg"), height = 5.08, width = 4.69)
}

p1 <- as.ggplot(XA_median_plots[["diploid"]])
p2 <- as.ggplot(XA_median_plots[["tet_P14"]])
p3 <- as.ggplot(XA_median_plots[["tet_P30"]])

combined_XA <-p1 + p2 + p3 + plot_layout(ncol = 1)
combined_XA
ggsave("combined_XA.svg", plot = combined_XA, device = "svg", width = 15, height = 30, units = "cm")

#stat tests
#Moods median test - TEST WHETHER medians of INDEPENDENT SAMPLES DIFFERENT
#nonparametric test - 

dip_X <- all_quart2 %>% filter(cat == "diploid") %>% filter(type == "X")
dip_X2 <- as.vector(dip_X[,4])
dip_A <- all_quart2 %>% filter(cat == "diploid") %>% filter(type == "A")
dip_A2 <- as.vector(dip_A[,4])
moods_dip <- mood.test(x = dip_X2, y = dip_A2, alternative = "two.sided")
p.value_dip <- moods_dip$p.value
MannU_dip <- wilcox.test(x = dip_X2, y = dip_A2, alternative = "two.sided")
p.value_MU_dip <- MannU_dip$p.value

p14_X <- all_quart2 %>% filter(cat == "tet_P14") %>% filter(type == "X")
p14_X2 <- as.vector(p14_X[,4])
p14_A <- all_quart2 %>% filter(cat == "tet_P14") %>% filter(type == "A")
p14_A2 <- as.vector(p14_A[,4])
moods_p14 <- mood.test(x = p14_X2, y = p14_A2, alternative = "two.sided")
p.value_p14 <- moods_p14$p.value
MannU_p14 <- wilcox.test(x = p14_X2, y = p14_A2, alternative = "two.sided")
p.value_MU_p14 <- MannU_p14$p.value

p30_X <- all_quart2 %>% filter(cat == "tet_P30") %>% filter(type == "X")
p30_X2 <- as.vector(p30_X[,4])
p30_A <- all_quart2 %>% filter(cat == "tet_P30") %>% filter(type == "A")
p30_A2 <- as.vector(p30_A[,4])
moods_p30 <- mood.test(x = p30_X2, y = p30_A2, alternative = "two.sided")
p.value_p30 <- moods_p30$p.value
MannU_p30 <- wilcox.test(x = p30_X2, y = p30_A2, alternative = "two.sided")
p.value_MU_p30 <- MannU_p30$p.value


#FIGURE 4-10 - ALTERNATIVE APPROACH
#Compare the output of each autosome within samples to see if there is a change of chromosome 2 compared to the other autosomes. 

#To compare output by chromosome within sample
setwd("/Users//R_files/outputs/expression_by_chr_within_sample/separate_A/")
#Arrange and set factor
  
all_quart3 <- as_tibble(all_quart) %>% arrange(chromosome) %>% mutate(chromosome = factor(chromosome, levels = c(1, 2, 3, 4, 5, 6, 7, 8,"X")))
#FOR individual autosomes
#FOR DIPLOID
  KW_dip <- rstatix::kruskal_test(data = all_quart3 %>% filter(chromosome != "zzz") %>% filter(cat == "diploid"), log2_cpm ~ chromosome)
  KW_dip$p <- format(KW_dip$p, scientific = TRUE, digits = 3)
  write.table(KW_dip, file="KW_dip_within_samp.txt", quote=FALSE, sep="\t", row.names = FALSE)
  
  DT_dip <- dunn_test(all_quart3 %>% filter(chromosome != "zzz") %>% filter(cat == "diploid"), log2_cpm ~ chromosome, p.adjust.method = "holm")
  DT_dip$p.adj <- format(DT_dip$p.adj, scientific = TRUE, digits = 3)
  DT_dip <- DT_dip %>% mutate(p.adj2 = (if_else(p.adj.signif == "ns", "", p.adj)))
  write.table(DT_dip, file="DT_dip_within_samp.txt", quote=FALSE, sep="\t", row.names = FALSE)
  DT_dip_filt <- DT_dip %>% filter(p.adj.signif != "ns")
  
  chr_median_plots_samp_dip <- ggboxplot(all_quart3 %>% filter(chromosome != "zzz") %>% filter(cat == "diploid"), x = "chromosome", y = "log2_cpm", color = "black",
                                         fill = "cat",
                                         palette = c("coral"), width = 0.4, alpha = 0.8,
                                         title = paste0("Diploid"),
                                         subtitle = paste0("Kruskal-Wallis p-value =", KW_dip$p),
                                         xlab = "Chromosome",
                                         ylab = "Log2 cpm",
                                         notch = TRUE) +
    theme(legend.position = "none")+  
    scale_x_discrete(breaks = c(1, 2, 3, 4, 5, 6, 7, 8,"X"), labels = c("1", "2", "3", "4", "5", "6", "7", "8", "X"))+
    stat_summary(fun=median, geom = "text", show.legend = FALSE, vjust=-0.5, size = 3.5, aes( label=round(..y.., digits=2)))+
    stat_pvalue_manual(DT_dip_filt, label = "{p.adj.signif}", y.position = c(13,14,15,16,17,18,19,20,21))+
    geom_text(stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=0.5, size=3.5)

  chr_median_plots_samp_dip
  ggsave(chr_median_plots_samp_dip, file = "chr_median_plots_samp_dip.svg", width = 7.85, height = 4.5, units = "in")
  
  #FOR P14
  KW_P14 <- rstatix::kruskal_test(data = all_quart3 %>% filter(chromosome != "zzz") %>% filter(cat == "tet_P14"), log2_cpm ~ chromosome)
  KW_P14$p <- format(KW_P14$p, scientific = TRUE, digits = 3)
  write.table(KW_P14, file="KW_P14_within_samp.txt", quote=FALSE, sep="\t", row.names = FALSE)
  
  
  DT_P14 <- dunn_test(all_quart3 %>% filter(chromosome != "zzz") %>% filter(cat == "tet_P14"), log2_cpm ~ chromosome, p.adjust.method = "holm")
  DT_P14$p.adj <- format(DT_P14$p.adj, scientific = TRUE, digits = 3)
  DT_P14 <- DT_P14 %>% mutate(p.adj2 = (if_else(p.adj.signif == "ns", "", p.adj)))
  write.table(DT_P14, file="DT_P14_within_samp.txt", quote=FALSE, sep="\t", row.names = FALSE)
  DT_P14_filt <- DT_P14 %>% filter(p.adj.signif != "ns")
  
  
  chr_median_plots_samp_P14 <- ggboxplot(all_quart3 %>% filter(chromosome != "zzz") %>% filter(cat == "tet_P14"), x = "chromosome", y = "log2_cpm", color = "black",
                                         fill = "cat",
                                         palette = c("chartreuse2"), width = 0.4, alpha = 0.8,
                                         title = paste0("Tetraploid P14"),
                                         subtitle = paste0("Kruskal-Wallis p-value =", KW_P14$p),
                                         xlab = "Chromosome",
                                         ylab = "Log2 cpm",
                                         notch = TRUE) +
    theme(legend.position = "none")+  
    scale_x_discrete(breaks = c(1, 2, 3, 4, 5, 6, 7, 8,"X"), labels = c("1", "2", "3", "4", "5", "6", "7", "8", "X"))+
    stat_summary(fun=median, geom = "text", show.legend = FALSE, vjust=-0.5, size = 3.5, aes( label=round(..y.., digits=2)))+
    stat_pvalue_manual(DT_P14_filt, label = "{p.adj.signif}", y.position = c(13,14,15,16))+
    geom_text(stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=0.5, size=3.5)
  
  chr_median_plots_samp_P14
  ggsave(chr_median_plots_samp_P14, file = "chr_median_plots_samp_P14.svg", width = 7.85, height = 4.4, units = "in")
  
  #FOR P30
  KW_P30 <- rstatix::kruskal_test(data = all_quart3 %>% filter(chromosome != "zzz") %>% filter(cat == "tet_P30"), log2_cpm ~ chromosome)
  KW_P30$p <- format(KW_P30$p, scientific = TRUE, digits = 3)
  write.table(KW_P30, file="KW_P30_within_samp.txt", quote=FALSE, sep="\t", row.names = FALSE)
  
  
  DT_P30 <- dunn_test(all_quart3 %>% filter(chromosome != "zzz") %>% filter(cat == "tet_P30"), log2_cpm ~ chromosome, p.adjust.method = "holm")
  DT_P30$p.adj <- format(DT_P30$p.adj, scientific = TRUE, digits = 3)
  DT_P30 <- DT_P30 %>% mutate(p.adj2 = (if_else(p.adj.signif == "ns", "", p.adj)))
  write.table(DT_P30, file="DT_P30_within_samp.txt", quote=FALSE, sep="\t", row.names = FALSE)
  DT_P30_filt <- DT_P30 %>% filter(p.adj.signif != "ns")
  
  chr_median_plots_samp_P30 <- ggboxplot(all_quart3 %>% filter(chromosome != "zzz") %>% filter(cat == "tet_P30"), x = "chromosome", y = "log2_cpm", color = "black",
                                         fill = "cat",
                                         palette = c("cornflowerblue"), width = 0.4, alpha = 0.8,
                                         title = paste0("Tetraploid P30"),
                                         subtitle = paste0("Kruskal-Wallis p-value =", KW_P30$p),
                                         xlab = "Chromosome",
                                         ylab = "Log2 cpm",
                                         notch = TRUE) +
    theme(legend.position = "none")+  
    scale_x_discrete(breaks = c(1, 2, 3, 4, 5, 6, 7, 8,"X"), labels = c("1", "2", "3", "4", "5", "6", "7", "8", "X"))+
    stat_summary(fun=median, geom = "text", show.legend = FALSE, vjust=-0.5, size = 3.5, aes( label=round(..y.., digits=2)))+
    stat_pvalue_manual(DT_P30_filt, label = "{p.adj.signif}", y.position = c(13,14,15,16,17,18,19,20,21,22,23))+
    geom_text(stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=0.5, size=3.5)
  
  chr_median_plots_samp_P30
  ggsave(chr_median_plots_samp_P30, file = "chr_median_plots_samp_P30.svg", width = 7.85, height = 4.5, units = "in")
  
  p1 <- as.ggplot(chr_median_plots_samp_dip)
  p2 <- as.ggplot(chr_median_plots_samp_P14)
  p3 <- as.ggplot(chr_median_plots_samp_P30)
  
  combined_XA_within_samp <-p1+ p2 + p3 + plot_layout(ncol = 1)
  combined_XA_within_samp
  ggsave("combined_XA_within_samp.svg", plot = combined_XA_within_samp, device = "svg", width = 25, height = 30, units = "cm")
  
#FIGURE 4-10A - FOR GROUPED autosomes - excluding 2
setwd("/Users/.../R_files/outputs/expression_by_chr_within_sample/combined_A/")  
all_quart4 <- all_quart3 %>% filter(chromosome != "zzz") %>% mutate(chromosome_sum = case_when(chromosome != "2" & chromosome != "X" ~ "other_autosomes",
                                                            chromosome == "2" ~ "2",
                                                            chromosome == "X" ~ "X"))
all_quart4 <- all_quart4 %>% arrange(chromosome_sum) %>% mutate(chromosome_sum = factor(chromosome_sum, levels = c(2, "other_autosomes","X")))
  
  
  #FOR DIPLOID - combined autosomes
  KW_dip_sum <- rstatix::kruskal_test(data = all_quart4 %>% filter(cat == "diploid"), log2_cpm ~ chromosome_sum)
  KW_dip_sum$p <- format(KW_dip_sum$p, scientific = TRUE, digits = 3)
  write.table(KW_dip_sum, file="KW_dip_sum_within_samp.txt", quote=FALSE, sep="\t", row.names = FALSE)
  
  DT_dip_sum <- dunn_test(all_quart4 %>% filter(cat == "diploid"), log2_cpm ~ chromosome_sum, p.adjust.method = "holm")
  DT_dip_sum$p.adj <- format(DT_dip_sum$p.adj, scientific = TRUE, digits = 3)
  DT_dip_sum <- DT_dip_sum %>% mutate(p.adj2 = (if_else(p.adj.signif == "ns", "", p.adj)))
  write.table(DT_dip_sum, file="DT_dip_sum_within_samp.txt", quote=FALSE, sep="\t", row.names = FALSE)
  
  chr_median_plots_samp_dip_sum <- ggboxplot(all_quart4 %>% filter(cat == "diploid"), x = "chromosome_sum", y = "log2_cpm", color = "black",
                                         fill = "cat",
                                         palette = c("coral"), width = 0.4, alpha = 0.8,
                                         title = paste0("Diploid"),
                                         subtitle = paste0("Kruskal-Wallis p-value =", KW_dip_sum$p),
                                         xlab = "",
                                         ylab = "Log2 cpm",
                                         notch = TRUE) +
    theme(legend.position = "none")+  
    scale_x_discrete(breaks = c("2", "other_autosomes","X"), labels = c("chr 2", "other autosomes","X chr"))+
    stat_summary(fun=median, geom = "text", show.legend = FALSE, vjust=-0.5, size = 3.5, aes( label=round(..y.., digits=2)))+
    stat_pvalue_manual(DT_dip_sum, label = "{p.adj.signif}", y.position = c(13.8,14.8,15.8))+
    geom_text(stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=0.5, size=3.5)
  
  chr_median_plots_samp_dip_sum
  ggsave(chr_median_plots_samp_dip_sum, file = "chr_median_plots_samp_dip_sum.svg", width = 4.5, height = 5.51, units = "in")
  
  #FOR P14 - combined autosomes
  KW_P14_sum <- rstatix::kruskal_test(data = all_quart4 %>% filter(cat == "tet_P14"), log2_cpm ~ chromosome_sum)
  KW_P14_sum$p <- format(KW_P14_sum$p, scientific = TRUE, digits = 3)
  write.table(KW_P14_sum, file="KW_P14_sum_within_samp.txt", quote=FALSE, sep="\t", row.names = FALSE)
  
  DT_P14_sum <- dunn_test(all_quart4 %>% filter(cat == "tet_P14"), log2_cpm ~ chromosome_sum, p.adjust.method = "holm")
  DT_P14_sum$p.adj <- format(DT_P14_sum$p.adj, scientific = TRUE, digits = 3)
  DT_P14_sum <- DT_P14_sum %>% mutate(p.adj2 = (if_else(p.adj.signif == "ns", "", p.adj)))
  write.table(DT_P14_sum, file="DT_P14_sum_within_samp.txt", quote=FALSE, sep="\t", row.names = FALSE)
  
  chr_median_plots_samp_P14_sum <- ggboxplot(all_quart4 %>% filter(cat == "tet_P14"), x = "chromosome_sum", y = "log2_cpm", color = "black",
                                             fill = "cat",
                                             palette = c("chartreuse2"), width = 0.4, alpha = 0.8,
                                             title = paste0("Tetraploid P14"),
                                             subtitle = paste0("Kruskal-Wallis p-value =", KW_P14_sum$p),
                                             xlab = "",
                                             ylab = "Log2 cpm",
                                             notch = TRUE) +
    theme(legend.position = "none")+  
    scale_x_discrete(breaks = c("2", "other_autosomes","X"), labels = c("chr 2", "other autosomes","X chr"))+
    stat_summary(fun=median, geom = "text", show.legend = FALSE, vjust=-0.5, size = 3.5, aes( label=round(..y.., digits=2)))+
    stat_pvalue_manual(DT_P14_sum, label = "{p.adj.signif}", y.position = c(13.8,14.8,15.8))+
    geom_text(stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=0.5, size=3.5)
  
  chr_median_plots_samp_P14_sum
  ggsave(chr_median_plots_samp_P14_sum, file = "chr_median_plots_samp_P14_sum.svg", width = 4.5, height = 5.51, units = "in")
  
  #FOR P30 - combined autosomes
  KW_P30_sum <- rstatix::kruskal_test(data = all_quart4 %>% filter(cat == "tet_P30"), log2_cpm ~ chromosome_sum)
  KW_P30_sum$p <- format(KW_P30_sum$p, scientific = TRUE, digits = 3)
  write.table(KW_P30_sum, file="KW_P30_sum_within_samp.txt", quote=FALSE, sep="\t", row.names = FALSE)
  
  DT_P30_sum <- dunn_test(all_quart4 %>% filter(cat == "tet_P30"), log2_cpm ~ chromosome_sum, p.adjust.method = "holm")
  DT_P30_sum$p.adj <- format(DT_P30_sum$p.adj, scientific = TRUE, digits = 3)
  DT_P30_sum <- DT_P30_sum %>% mutate(p.adj2 = (if_else(p.adj.signif == "ns", "", p.adj)))
  write.table(DT_P30_sum, file="DT_P30_sum_within_samp.txt", quote=FALSE, sep="\t", row.names = FALSE)
  
  chr_median_plots_samp_P30_sum <- ggboxplot(all_quart4 %>% filter(cat == "tet_P30"), x = "chromosome_sum", y = "log2_cpm", color = "black",
                                             fill = "cat",
                                             palette = c("cornflowerblue"), width = 0.4, alpha = 0.8,
                                             title = paste0("Tetraploid P30"),
                                             subtitle = paste0("Kruskal-Wallis p-value =", KW_P30_sum$p),
                                             xlab = "",
                                             ylab = "Log2 cpm",
                                             notch = TRUE) +
    theme(legend.position = "none")+  
    scale_x_discrete(breaks = c("2", "other_autosomes","X"), labels = c("chr 2", "other autosomes","X chr"))+
    stat_summary(fun=median, geom = "text", show.legend = FALSE, vjust=-0.5, size = 3.5, aes( label=round(..y.., digits=2)))+
    stat_pvalue_manual(DT_P30_sum, label = "{p.adj.signif}", y.position = c(13.8,14.8,15.8))+
    geom_text(stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=0.5, size=3.5)
  
  chr_median_plots_samp_P30_sum
  ggsave(chr_median_plots_samp_P30_sum, file = "chr_median_plots_samp_P30_sum.svg", width = 4.5, height = 5.51, units = "in")
  
  p1 <- as.ggplot(chr_median_plots_samp_dip_sum)
  p2 <- as.ggplot(chr_median_plots_samp_P14_sum)
  p3 <- as.ggplot(chr_median_plots_samp_P30_sum)
  
  combined_XA_within_samp <-p1 + p2 + p3 + plot_layout(ncol = 3)
  combined_XA_within_samp
  ggsave("combined_XA_within_samp_sum.svg", plot = combined_XA_within_samp, device = "svg", width = 35, height = 12, units = "cm")
  
