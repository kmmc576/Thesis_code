#FIGURE 4-5C
library(tidyverse)
library(fabricatr)
library(ggplotify)
library(patchwork)
library(viridis)
library(plotly)
setwd("/Users/.../MDO_tetra_5mC/")
tetra_spreads <- read.delim2("/Users/.../MDO_tetra_5mC/tetraploidy_summary.txt", header=TRUE, sep="\t")

tetra_spreads <- tetra_spreads %>% arrange(chr_num, Xchr_count)
tetra_spreads$Xchr_count <- as.character(tetra_spreads$Xchr_count)
tetra_spreads$nucleus <- factor(tetra_spreads$nucleus, levels = unique(tetra_spreads$nucleus))

plot_tetra2 <- ggplot(tetra_spreads, aes(chr_num))+
  geom_histogram(alpha = 0.7, binwidth = 0.5, colour = "blue")+
 # scale_fill_manual(values = c("#0072B2", "#00FF33", "#FF9000", "#FF66CC", "#33CCFF", "#FF0000","#999999"))+
  scale_x_continuous(limits = c(10,80), 
                     breaks = c(18,27,36,45,72),
                     labels = c("18(2N)","27", "36(4N)", "45", "72(8N)"))+
  scale_y_continuous(limits = c(0,100), 
                     breaks = c(0,20,40,60,80),
                     labels = c("0","20","40","60","80"))+
  labs(title="Tetraploid nuclei (P21 - P28)")+
  labs(subtitle="Chromosome number per nucleus (n = 193)")+
  theme_classic()+
  ylab("number of nuclei")+
  xlab("chromosome number")
plot_tetra2

ggsave("/Users/kimmcintyre/Library/CloudStorage/OneDrive-UNSW/Tetraploid chapter/Metaphase_spread_plots/histogram_metaphase_spread_counts.tiff",
       width = 10, height = 12, units = "cm")
chr_count <- tetra_spreads %>% count(chr_num)
write.table(chr_count, file="/Users/kimmcintyre/Library/CloudStorage/OneDrive-UNSW/Tetraploid chapter/Metaphase_spread_plots/chromosome_count.txt", quote=FALSE, sep="\t", row.names = FALSE)

tetra_spreads2 <- tetra_spreads %>% filter(chr_num >=31)  %>% filter(chr_num <=40) 
count2 <- tetra_spreads2 %>% count(missing_c2)
count3 <- tetra_spreads2 %>% count(Xchr_count)
count4 <- tetra_spreads2 %>% count(extra.chr)
