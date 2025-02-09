#FIGURE 4-18
##read in excel doc: 5mC MDO tetra July 2021 (Olympus BX53) - saved as tab delimited text
#X intensity normalised by mean autosomal intensity for that spread.
#mean autosomal intensity calculated by deducting area, and total intensity for all Xs from area and total intensity for whole spread 

library(ggplot2)
library(tidyverse)
library(gt)
library(rstatix)
library(ggpubr)

setwd("/Users/.../R_files/")


data <- read.csv("5mC_intensity/5mC_MDOtetra_Jul2021.txt", header = T, stringsAsFactors = F, sep = '\t')
P22 <- data %>% filter(Passage==22)
P28 <- data %>% filter(Passage==28)
P29 <- data %>% filter(Passage==29)

##pulls file names down in data
name <- ''
for(i in 1:nrow(data)){
  if(data[i, "File_name"] == ''){
    data[i, "File_name"] <-  name
  } else {
    name <- data[i, "File_name"]
  }
} 

##add column log2(X_int_norm)
lg2X_int_norm <- c(log2(data$X_int_norm))
data <- cbind(data, lg2X_int_norm)

##SUB DATA FRAMES
#create sub data frame, rank entries for each image, each sorted on log2_X_int_norm
data_lg2Xnorm_auto <- data.frame()
for(i in unique(data$File_name)){
  
  data_lg2Xnorm_auto <- rbind(data_lg2Xnorm_auto, data.frame(name = i, rank = c(1:4), lg2X_int_norm = sort(data[data$File_name==i, "lg2X_int_norm"])))
  
}
data_lg2Xnorm_auto$rank <- as.character(data_lg2Xnorm_auto$rank)
#add in passage number
pass <- data %>% dplyr::select("name" = "File_name", "Passage")
data_lg2Xnorm_auto <- left_join(data_lg2Xnorm_auto, pass, by = "name", keep = FALSE)
data_lg2Xnorm_auto <- unique(data_lg2Xnorm_auto)


##BOX PLOTS
#Hinges at 1st and 3rd quartiles. Whiskers default to 1.5X interquartile range. Bar is median
##1) makes a box plot with data_Xnorm_auto

#Notched boxplot
#Notches are used to compare groups; if the notches of two boxes do not overlap, this suggests that the medians are significantly different.
#delete legend.position=none for full legend (cell names)
#delete notch=TRUE for normal box plot

#Use ggboxplot to plot p values
#Do anova, then  TukeyHSD


stat.test <- aov(lg2X_int_norm ~ rank, data = data_lg2Xnorm_auto) %>% tukey_hsd() 
stat.test$p.scientif <- format(stat.test$p.adj, scientific = TRUE, digits = 3)

ggboxplot_all <- ggboxplot(data_lg2Xnorm_auto, x = "rank", y = "lg2X_int_norm", color = "black",
                           notch = TRUE,
                           fill = "rank",
                         palette = c("#80BEE2", "#80BEE2", "#80BEE2", "#80BEE2"),
                           title = "",
                         legend = "none",
                           xlab = "Rank (X intensity)",
                           ylab = "Log2(X intensity)") +
  scale_x_discrete(breaks = c("1", "2", "3", "4"), labels = c("1", "2", "3", "4"))+
  stat_pvalue_manual(
    stat.test, label = "p.adj.signif",
    y.position = c(0.6, 0.8, 1, 1.2, 1.4, 1.6))
ggboxplot_all

ggsave(file = "ggboxplot_all.png", plot = ggboxplot_all, path = "outputs/5mC", width = 10, height = 7)
ggsave(file = "ggboxplot_all.svg", plot = ggboxplot_all, path = "outputs/5mC", width = 10, height = 7)


##ANOVA anova analysis (to assess whether diff btwn means of 3 or more groups)
##do one-way anova for single variable 
##NB if there is difference between means, this test will not tell which mean is different

anova_lg2Xnorm_auto <- aov(data_lg2Xnorm_auto$lg2X_int_norm ~ data_lg2Xnorm_auto$rank, data = data)
anova_lg2Xnorm_sum <- summary(anova_lg2Xnorm_auto)
anova_lg2Xnorm_sum <- as.data.frame(anova_lg2Xnorm_sum[[1]])

##do Tukey test to see which group has different mean, ie pairwise analysis
##gives P value adjusted for multiple groups
Tukey_lg2Xnorm <- TukeyHSD(anova_lg2Xnorm_auto)
Tukey_lg2Xnorm_sum <- as.data.frame(Tukey_lg2Xnorm[1])
colnames(Tukey_lg2Xnorm_sum) <- c("diff", "lower", "upper", "adj_pvalue")
Tukey_lg2Xnorm_sum <- Tukey_lg2Xnorm_sum %>% mutate("pair" = row.names(Tukey_lg2Xnorm_sum), .before = diff)

#write in tables
  Tukey_lg2Xnorm_sum_table <- Tukey_lg2Xnorm_sum %>% gt() %>%
    tab_header(title = "TukeyHSD log2 norm X intensity (aov p = 2.15e-38)") %>%
    fmt_scientific(diff, decimals = 2) %>%
    fmt_scientific(lower, decimals = 2) %>%
    fmt_scientific(upper, decimals = 2) %>%
    fmt_scientific(adj_pvalue, decimals = 2) %>%
    cols_align(align = c("left")) %>%
    tab_options(heading.align = "left")
  Tukey_lg2Xnorm_sum_table

gtsave(Tukey_lg2Xnorm_sum_table, "outputs/5mC/Tukey_lg2Xnorm_sum_table.png")


