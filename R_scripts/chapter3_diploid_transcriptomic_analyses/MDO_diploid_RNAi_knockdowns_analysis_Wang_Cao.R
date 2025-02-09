

library("edgeR")
library(tidyverse)
library(dplyr)
library(splitstackshape)
library(plotly)
library(rcompanion)
library(ggplot2)
library(stats)
library(ggpubr)
library(gprofiler2)
library(rstatix)
library(svglite)
library(GSA)


#comparison: CAO/WANG MDO SNP (2014) major allele proportions
#LIMIT ONLY TO RECIPROCAL F1 HYBRID CROSSES:A0571_1, A0571_4, A0579_3, A0579_4
#Autosomal data published with Cao(2023)

setwd("/Users/.../Wang_Cao_reanalysis/")

#Analyse separately - data from brain and placenta
#filter out intergenic and intronic SNVs

#separate brain
brain <- read.table("/Users/.../Wang_Cao_reanalysis/Cao_MDO_brain_SNVs.txt", header=TRUE, sep="\t")

brain <- brain %>% filter(chr != "chrUn") %>% filter(Effect != "INTRON") %>% filter(Effect != "INTER")

brain_effect <- brain %>% dplyr::count(Effect)

#filter for minimum counts, then calculate major allele proportion
#individual SNVs
brain2 <- brain %>% filter(Gene_name != "") %>% 
  filter(A0579_b3_ref_c >5 | A0579_b3_ref_c >5 | A0579_b4_alter_c >5 | A0579_b4_ref_c >5 |
           A0571_b1_alter_c >5 | A0571_b1_ref_c > 5 | A0571_b4_ref_c >5 | A0571_b4_alter_c >5 | A0580_b1_alter_c >5 | A0580_b1_ref_c >5 |
           A0580_b5_ref_c >5 | A0580_b5_alter_c >5 | A0572_b1_ref_c >5 | A0572_b1_alter_c >5 | A0572_b3_ref_c >5 | A0572_b3_alter_c >5) %>%
  mutate(A0579b3 = case_when(A0579_b3_ref_c > A0579_b3_alter_c ~ A0579_b3_ref_c/(A0579_b3_ref_c+A0579_b3_alter_c), A0579_b3_ref_c < A0579_b3_alter_c ~ A0579_b3_alter_c/(A0579_b3_ref_c+A0579_b3_alter_c), A0579_b3_ref_c == A0579_b3_alter_c ~ 0.5)) %>%
  mutate(A0579b4 = case_when(A0579_b4_ref_c > A0579_b4_alter_c ~ A0579_b4_ref_c/(A0579_b4_ref_c+A0579_b4_alter_c), A0579_b4_ref_c < A0579_b4_alter_c ~ A0579_b4_alter_c/(A0579_b4_ref_c+A0579_b4_alter_c), A0579_b4_ref_c == A0579_b4_alter_c ~ 0.5)) %>%
  mutate(A0571b1 = case_when(A0571_b1_ref_c > A0571_b1_alter_c ~ A0571_b1_ref_c/(A0571_b1_ref_c+A0571_b1_alter_c), A0571_b1_ref_c < A0571_b1_alter_c ~ A0571_b1_alter_c/(A0571_b1_ref_c+A0571_b1_alter_c), A0571_b1_ref_c == A0571_b1_alter_c ~ 0.5)) %>%
  mutate(A0571b4 = case_when(A0571_b4_ref_c > A0571_b4_alter_c ~ A0571_b4_ref_c/(A0571_b4_ref_c+A0571_b4_alter_c), A0571_b4_ref_c < A0571_b4_alter_c ~ A0571_b4_alter_c/(A0571_b4_ref_c+A0571_b4_alter_c), A0571_b4_ref_c == A0571_b4_alter_c ~ 0.5)) %>%
  mutate(A0580b1 = case_when(A0580_b1_ref_c > A0580_b1_alter_c ~ A0580_b1_ref_c/(A0580_b1_ref_c+A0580_b1_alter_c), A0580_b1_ref_c < A0580_b1_alter_c ~ A0580_b1_alter_c/(A0580_b1_ref_c+A0580_b1_alter_c), A0580_b1_ref_c == A0580_b1_alter_c ~ 0.5)) %>%
  mutate(A0580b5 = case_when(A0580_b5_ref_c > A0580_b5_alter_c ~ A0580_b5_ref_c/(A0580_b5_ref_c+A0580_b5_alter_c), A0580_b5_ref_c < A0580_b5_alter_c ~ A0580_b5_alter_c/(A0580_b5_ref_c+A0580_b5_alter_c), A0580_b5_ref_c == A0580_b5_alter_c ~ 0.5)) %>%
  mutate(A0572b1 = case_when(A0572_b1_ref_c > A0572_b1_alter_c ~ A0572_b1_ref_c/(A0572_b1_ref_c+A0572_b1_alter_c), A0572_b1_ref_c < A0572_b1_alter_c ~ A0572_b1_alter_c/(A0572_b1_ref_c+A0572_b1_alter_c), A0572_b1_ref_c == A0572_b1_alter_c ~ 0.5)) %>%
  mutate(A0572b3 = case_when(A0572_b3_ref_c > A0572_b3_alter_c ~ A0572_b3_ref_c/(A0572_b3_ref_c+A0572_b3_alter_c), A0572_b3_ref_c < A0572_b3_alter_c ~ A0572_b3_alter_c/(A0572_b3_ref_c+A0572_b3_alter_c), A0572_b3_ref_c == A0572_b3_alter_c ~ 0.5)) %>%
  dplyr::select(chr, position, Gene_name, A0579b3, A0579b4, A0571b1, A0571b4, A0580b1, A0580b5, A0572b1, A0572b3)

brain3 <- brain2 %>% pivot_longer(cols = c(A0579b3, A0579b4, A0571b1, A0571b4, A0580b1, A0580b5, A0572b1, A0572b3), names_to = "sample", values_to = "map")

brain_plot <- ggplot(brain3, aes(x = sample, y = map, fill = sample), alpha = 0.2)+
  geom_boxplot(notch = TRUE, outliers = FALSE)+
  xlab("")+
  ylab("")+
  theme_classic(base_size = 18)+
  theme(legend.position = "none")+
  geom_text(data = . %>% dplyr::count(sample), aes(y = 0.4, label = n), size = 5, position = "dodge")+
  stat_summary(data = brain3, fun = "median", geom = "text", aes(label = round(after_stat(y),3), vjust = -1))+
  theme(legend.title = element_blank())
brain_plot
ggsave(file = "Cao_brain_A_separate_snvs.svg", plot = brain_plot,  device = "svg", width = 9, height = 5, units = "in")

pal <- hue_pal()(8)
pal #"#F8766D" "#CD9600" "#7CAE00" "#00BE67" "#00BFC4" "#00A9FF" "#C77CFF" "#FF61CC"

brain_plot_density <- ggplot()+
  #geom_histogram(bins = 60, fill = "dodgerblue3", alpha = 0.7)+
  geom_density(data = brain3 %>% filter(sample == "A0579b3"), aes(x = map, y = ..density..), colour = "#F8766D")+
  geom_density(data = brain3 %>% filter(sample == "A0579b4"), aes(x = map, y = ..density..), colour = "#CD9600")+
  geom_density(data = brain3 %>% filter(sample == "A0571b1"), aes(x = map, y = ..density..), colour = "#7CAE00")+
  geom_density(data = brain3 %>% filter(sample == "A0571b4"), aes(x = map, y = ..density..), colour = "#00BE67")+
  geom_density(data = brain3 %>% filter(sample == "A0580b1"), aes(x = map, y = ..density..), colour = "#00BFC4")+
  geom_density(data = brain3 %>% filter(sample == "A0580b5"), aes(x = map, y = ..density..), colour = "#00A9FF")+
  geom_density(data = brain3 %>% filter(sample == "A0572b1"), aes(x = map, y = ..density..), colour = "#C77CFF")+
  geom_density(data = brain3 %>% filter(sample == "A0572b3"), aes(x = map, y = ..density..), colour = "#FF61CC")+
  theme_classic(base_size = 18)+
  xlab("")+
  ylab("")+ #gene number
  theme(legend.position = "none")
brain_plot_density
ggsave(file = "Cao_brain_A_separate_density.svg", plot = brain_plot_density,  device = "svg", width = 4, height = 4, units = "in")

#PREFER to use gene median major allele proportion:
#gene median map

brain4 <- brain2 %>% group_by(Gene_name) %>% mutate(A0579b3_m = median(A0579b3)) %>% mutate(A0579b4_m = median(A0579b4)) %>% mutate(A0571b1_m = median(A0571b1)) %>%
  mutate(A0571b4_m = median(A0571b4)) %>% mutate(A0580b1_m = median(A0580b1)) %>% mutate(A0580b5_m = median(A0580b5)) %>% mutate(A0572b1_m = median(A0572b1)) %>% mutate(A0572b3_m = median(A0572b3)) %>% ungroup 

brain5 <- brain4 %>% dplyr::select(-4, -5, -6, -7, -8, -9, -10, -11) %>% distinct(chr, Gene_name, .keep_all = TRUE) %>%
  pivot_longer(cols = c(A0579b3_m, A0579b4_m, A0571b1_m, A0571b4_m, A0580b1_m, A0580b5_m, A0572b1_m, A0572b3_m), names_to = "sample", values_to = "map")

brain_plot_gene <- ggplot(brain5, aes(x = sample, y = map, fill = sample), alpha = 0.2)+
  geom_boxplot(notch = TRUE, outliers = FALSE)+
  xlab("")+
  ylab("")+
  theme_classic(base_size = 18)+
  theme(legend.position = "none")+
  geom_text(data = . %>% dplyr::count(sample), aes(y = 0.4, label = n), size = 5, position = "dodge")+
  stat_summary(data = brain5, fun = "median", geom = "text", aes(label = round(after_stat(y),3), vjust = -1))+
  theme(legend.title = element_blank())
brain_plot_gene
ggsave(file = "Cao_brain_A_separate_gene.svg", plot = brain_plot_gene,  device = "svg", width = 10, height = 5, units = "in")

brain_plot_gene_density <- ggplot()+
  #geom_histogram(bins = 60, fill = "dodgerblue3", alpha = 0.7)+
  geom_density(data = brain5 %>% filter(sample == "A0579b3_m"), aes(x = map, y = ..density..), colour = "#F8766D")+
  geom_density(data = brain5 %>% filter(sample == "A0579b4_m"), aes(x = map, y = ..density..), colour = "#CD9600")+
  geom_density(data = brain5 %>% filter(sample == "A0571b1_m"), aes(x = map, y = ..density..), colour = "#7CAE00")+
  geom_density(data = brain5 %>% filter(sample == "A0571b4_m"), aes(x = map, y = ..density..), colour = "#00BE67")+
  geom_density(data = brain5 %>% filter(sample == "A0580b1_m"), aes(x = map, y = ..density..), colour = "#00BFC4")+
  geom_density(data = brain5 %>% filter(sample == "A0580b5_m"), aes(x = map, y = ..density..), colour = "#00A9FF")+
  geom_density(data = brain5 %>% filter(sample == "A0572b1_m"), aes(x = map, y = ..density..), colour = "#C77CFF")+
  geom_density(data = brain5 %>% filter(sample == "A0572b3_m"), aes(x = map, y = ..density..), colour = "#FF61CC")+
  theme_classic(base_size = 18)+
  xlab("")+
  ylab("")+ #gene number
  theme(legend.position = "none")
brain_plot_gene_density
ggsave(file = "Cao_brain_A_separate_gene_density.svg", plot = brain_plot_gene_density,  device = "svg", width = 5, height = 4, units = "in")

#separate placenta
plac <- read.table("/Users/.../Wang_Cao_reanalysis/Cao_MDO_placenta_SNVs.txt", header=TRUE, sep="\t")

plac <- plac %>% filter(chr != "chrUn") %>% filter(Effect != "INTRON") %>% filter(Effect != "INTER")

plac_effect <- plac %>% dplyr::count(Effect)

plac2 <- plac %>% filter(Gene_name != "") %>% 
  filter(A0579_p3_ref_c >5 | A0579_p3_ref_c >5 | A0579_p4_alter_c >5 | A0579_p4_ref_c >5 |
           A0571_p1_alter_c >5 | A0571_p1_ref_c > 5 | A0571_p4_ref_c >5 | A0571_p4_alter_c >5 | A0580_p1_alter_c >5 | A0580_p1_ref_c >5 |
           A0580_p5_ref_c >5 | A0580_p5_alter_c >5 | A0572_p1_ref_c >5 | A0572_p1_alter_c >5 | A0572_p3_ref_c >5 | A0572_p3_alter_c >5) %>%
  mutate(A0579p3 = case_when(A0579_p3_ref_c > A0579_p3_alter_c ~ A0579_p3_ref_c/(A0579_p3_ref_c+A0579_p3_alter_c), A0579_p3_ref_c < A0579_p3_alter_c ~ A0579_p3_alter_c/(A0579_p3_ref_c+A0579_p3_alter_c), A0579_p3_ref_c == A0579_p3_alter_c ~ 0.5)) %>%
  mutate(A0579p4 = case_when(A0579_p4_ref_c > A0579_p4_alter_c ~ A0579_p4_ref_c/(A0579_p4_ref_c+A0579_p4_alter_c), A0579_p4_ref_c < A0579_p4_alter_c ~ A0579_p4_alter_c/(A0579_p4_ref_c+A0579_p4_alter_c), A0579_p4_ref_c == A0579_p4_alter_c ~ 0.5)) %>%
  mutate(A0571p1 = case_when(A0571_p1_ref_c > A0571_p1_alter_c ~ A0571_p1_ref_c/(A0571_p1_ref_c+A0571_p1_alter_c), A0571_p1_ref_c < A0571_p1_alter_c ~ A0571_p1_alter_c/(A0571_p1_ref_c+A0571_p1_alter_c), A0571_p1_ref_c == A0571_p1_alter_c ~ 0.5)) %>%
  mutate(A0571p4 = case_when(A0571_p4_ref_c > A0571_p4_alter_c ~ A0571_p4_ref_c/(A0571_p4_ref_c+A0571_p4_alter_c), A0571_p4_ref_c < A0571_p4_alter_c ~ A0571_p4_alter_c/(A0571_p4_ref_c+A0571_p4_alter_c), A0571_p4_ref_c == A0571_p4_alter_c ~ 0.5)) %>%
  mutate(A0580p1 = case_when(A0580_p1_ref_c > A0580_p1_alter_c ~ A0580_p1_ref_c/(A0580_p1_ref_c+A0580_p1_alter_c), A0580_p1_ref_c < A0580_p1_alter_c ~ A0580_p1_alter_c/(A0580_p1_ref_c+A0580_p1_alter_c), A0580_p1_ref_c == A0580_p1_alter_c ~ 0.5)) %>%
  mutate(A0580p5 = case_when(A0580_p5_ref_c > A0580_p5_alter_c ~ A0580_p5_ref_c/(A0580_p5_ref_c+A0580_p5_alter_c), A0580_p5_ref_c < A0580_p5_alter_c ~ A0580_p5_alter_c/(A0580_p5_ref_c+A0580_p5_alter_c), A0580_p5_ref_c == A0580_p5_alter_c ~ 0.5)) %>%
  mutate(A0572p1 = case_when(A0572_p1_ref_c > A0572_p1_alter_c ~ A0572_p1_ref_c/(A0572_p1_ref_c+A0572_p1_alter_c), A0572_p1_ref_c < A0572_p1_alter_c ~ A0572_p1_alter_c/(A0572_p1_ref_c+A0572_p1_alter_c), A0572_p1_ref_c == A0572_p1_alter_c ~ 0.5)) %>%
  mutate(A0572p3 = case_when(A0572_p3_ref_c > A0572_p3_alter_c ~ A0572_p3_ref_c/(A0572_p3_ref_c+A0572_p3_alter_c), A0572_p3_ref_c < A0572_p3_alter_c ~ A0572_p3_alter_c/(A0572_p3_ref_c+A0572_p3_alter_c), A0572_p3_ref_c == A0572_p3_alter_c ~ 0.5)) %>%
  dplyr::select(chr, position, Gene_name, A0579p3, A0579p4, A0571p1, A0571p4, A0580p1, A0580p5, A0572p1, A0572p3)

plac3 <- plac2 %>% pivot_longer(cols = c(A0579p3, A0579p4, A0571p1, A0571p4, A0580p1, A0580p5, A0572p1, A0572p3), names_to = "sample", values_to = "map")

plac_plot <- ggplot(plac3, aes(x = sample, y = map, fill = sample), alpha = 0.2)+
  geom_boxplot(notch = TRUE, outliers = FALSE)+
  xlab("")+
  ylab("")+
  theme_classic(base_size = 18)+
  theme(legend.position = "none")+
  geom_text(data = . %>% dplyr::count(sample), aes(y = 0.4, label = n), size = 5, position = "dodge")+
  stat_summary(data = plac3, fun = "median", geom = "text", aes(label = round(after_stat(y),3), vjust = -1))+
  theme(legend.title = element_blank())
plac_plot
ggsave(file = "Cao_plac_A_separate_snvs.svg", plot = plac_plot,  device = "svg", width = 10, height = 5, units = "in")


plac_plot_density <- ggplot()+
  #geom_histogram(bins = 60, fill = "dodgerblue3", alpha = 0.7)+
  geom_density(data = plac3 %>% filter(sample == "A0579p3"), aes(x = map, y = ..density..), colour = "#F8766D")+
  geom_density(data = plac3 %>% filter(sample == "A0579p4"), aes(x = map, y = ..density..), colour = "#CD9600")+
  geom_density(data = plac3 %>% filter(sample == "A0571p1"), aes(x = map, y = ..density..), colour = "#7CAE00")+
  geom_density(data = plac3 %>% filter(sample == "A0571p4"), aes(x = map, y = ..density..), colour = "#00BE67")+
  geom_density(data = plac3 %>% filter(sample == "A0580p1"), aes(x = map, y = ..density..), colour = "#00BFC4")+
  geom_density(data = plac3 %>% filter(sample == "A0580p5"), aes(x = map, y = ..density..), colour = "#00A9FF")+
  geom_density(data = plac3 %>% filter(sample == "A0572p1"), aes(x = map, y = ..density..), colour = "#C77CFF")+
  geom_density(data = plac3 %>% filter(sample == "A0572p3"), aes(x = map, y = ..density..), colour = "#FF61CC")+
  theme_classic(base_size = 18)+
  xlab("")+
  ylab("")+ #gene number
  theme(legend.position = "none")
plac_plot_density
ggsave(file = "Cao_plac_A_separate_density.svg", plot = plac_plot_density,  device = "svg", width = 4, height = 4, units = "in")

#As above, prefer to use gene median major allel proportion

plac4 <- plac2 %>% group_by(Gene_name) %>% mutate(A0579p3_m = median(A0579p3)) %>% mutate(A0579p4_m = median(A0579p4)) %>% mutate(A0571p1_m = median(A0571p1)) %>%
  mutate(A0571p4_m = median(A0571p4)) %>% mutate(A0580p1_m = median(A0580p1)) %>% mutate(A0580p5_m = median(A0580p5)) %>% mutate(A0572p1_m = median(A0572p1)) %>% mutate(A0572p3_m = median(A0572p3)) %>% ungroup 

plac5 <- plac4 %>% dplyr::select(-4, -5, -6, -7, -8, -9, -10, -11) %>% distinct(chr, Gene_name, .keep_all = TRUE) %>%
  pivot_longer(cols = c(A0579p3_m, A0579p4_m, A0571p1_m, A0571p4_m, A0580p1_m, A0580p5_m, A0572p1_m, A0572p3_m), names_to = "sample", values_to = "map")

plac_plot_gene <- ggplot(plac5, aes(x = sample, y = map, fill = sample), alpha = 0.2)+
  geom_boxplot(notch = TRUE, outliers = FALSE)+
  xlab("")+
  ylab("")+
  theme_classic(base_size = 18)+
  theme(legend.position = "none")+
  geom_text(data = . %>% dplyr::count(sample), aes(y = 0.4, label = n), size = 5, position = "dodge")+
  stat_summary(data = plac5, fun = "median", geom = "text", aes(label = round(after_stat(y),3), vjust = -1))+
  theme(legend.title = element_blank())
plac_plot_gene
ggsave(file = "Cao_plac_A_separate_gene.svg", plot = plac_plot_gene,  device = "svg", width = 10, height = 5, units = "in")

plac_plot_gene_density <- ggplot()+
  #geom_histogram(bins = 60, fill = "dodgerblue3", alpha = 0.7)+
  geom_density(data = plac5 %>% filter(sample == "A0579p3_m"), aes(x = map, y = ..density..), colour = "#F8766D")+
  geom_density(data = plac5 %>% filter(sample == "A0579p4_m"), aes(x = map, y = ..density..), colour = "#CD9600")+
  geom_density(data = plac5 %>% filter(sample == "A0571p1_m"), aes(x = map, y = ..density..), colour = "#7CAE00")+
  geom_density(data = plac5 %>% filter(sample == "A0571p4_m"), aes(x = map, y = ..density..), colour = "#00BE67")+
  geom_density(data = plac5 %>% filter(sample == "A0580p1_m"), aes(x = map, y = ..density..), colour = "#00BFC4")+
  geom_density(data = plac5 %>% filter(sample == "A0580p5_m"), aes(x = map, y = ..density..), colour = "#00A9FF")+
  geom_density(data = plac5 %>% filter(sample == "A0572p1_m"), aes(x = map, y = ..density..), colour = "#C77CFF")+
  geom_density(data = plac5 %>% filter(sample == "A0572p3_m"), aes(x = map, y = ..density..), colour = "#FF61CC")+
  theme_classic(base_size = 18)+
  xlab("")+
  ylab("")+ #gene number
  theme(legend.position = "none")
plac_plot_gene_density
ggsave(file = "Cao_plac_A_separate_gene_density.svg", plot = plac_plot_gene_density,  device = "svg", width = 4, height = 4, units = "in")


#Proportion skewed - based on gene median major allele proportion:
skew_brain <- brain5 %>% group_by(Gene_name, sample) %>% mutate(med_map = median(map)) %>% ungroup() %>%
  mutate(cat = (case_when(med_map >= 0.50 & med_map < 0.70 ~ "balanced", 
                          (med_map >= 0.70 & med_map <= 0.85) ~ "inter",
                          med_map > 0.85 ~ "skewed"))) %>% distinct(Gene_name, sample, .keep_all = TRUE)

skew_brain_count <- skew_brain %>% dplyr::count(sample, cat) %>% group_by(sample) %>% mutate(prop = n/sum(n))

skew_brain_plot <- ggplot(skew_brain_count, aes(x=sample, y=n, fill=cat))+
  geom_col() +
  coord_flip() +
  ylab("number of genes") +
  xlab("") +
  theme_classic(base_size = 16) +
  theme(legend.title = element_blank())
skew_brain_plot

ggsave(file = "skew_brain_plot.svg", plot = skew_brain_plot,  device = "svg", width = 5, height = 3, units = "in")


