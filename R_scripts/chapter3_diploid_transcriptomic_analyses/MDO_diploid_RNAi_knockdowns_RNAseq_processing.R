

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

setwd("/Users/.../Featurecounts_subread_mMonDom1/")

#read in counts from featurecounts output:

fc_u6m2 <- read.table("subread_fc_U6_diploid.txt",header=TRUE, sep="\t")
fc_hnrnpk <- read.table("subread_fc_Hnrnpk_diploid.txt",header=TRUE, sep="\t")
fc_ckap4 <- read.table("subread_fc_Ckap4_diploid.txt",header=TRUE, sep="\t")
fc_syncrip <- read.table("subread_fc_Syncrip_diploid.txt",header=TRUE, sep="\t")
fc_caprin1 <- read.table("subread_fc_Caprin_diploid.txt",header=TRUE, sep="\t")
fc_nono <- read.table("subread_fc_Nono_diploid.txt",header=TRUE, sep="\t")
fc_ncl <- read.table("subread_fc_Ncl_diploid.txt",header=TRUE, sep="\t")
fc_nude <- read.table("subread_fc_nude_diploid.txt",header=TRUE, sep="\t")

#combine all into single matrix
all_counts <- cbind(fc_caprin1[1:29466,], fc_ckap4[1:29466,7], fc_hnrnpk[1:29466,7], fc_ncl[1:29466,7], fc_nono[1:29466,7], fc_syncrip[1:29466,7], fc_u6m2[1:29466,7], fc_nude[1:29466,7])

colnames(all_counts) <- c("gene", "chromosome", "start", "end", "strand", "length", "CAPRIN1", "CKAP4", "HNRNPK", "NCL", "NONO", "SYNCRIP", "U6M2", "NUDE")

#convert to DGEList object
#rows = features, columns = samples
all_counts_DGE <- DGEList(counts=all_counts[ ,7:14], gene=all_counts[ ,1])

#normalise TMM (normalises to calculate effective library sizes). Does not change count values 
all_counts_DGE_norm <- calcNormFactors(all_counts_DGE)
head(all_counts_DGE_norm)
#group        lib.size      norm.factors
#CAPRIN1     1 124760386    0.9895571
#CKAP4       1  93197944    1.0042095
#HNRNPK      1  81126171    0.9843531
#NCL         1  82641814    1.0117585
#NONO        1  91071005    0.9920451
#SYNCRIP     1  79285005    0.9949529
#U6M2        1  73597729    0.9977603
#NUDE        1  82055005    1.0259989

#normalise using cpm 
keep <- rowSums(cpm(all_counts_DGE_norm)>1) >= 1
summary(keep)
#Mode   FALSE    TRUE 
#logical   15608   13858

#don't keep library sizes, ie so re-calculated after filtering for lowly expressed genes
all_counts_DGE_filtered <- all_counts_DGE_norm[keep, , keep.lib.sizes=FALSE]
head(all_counts_DGE_filtered)

all_counts_DGE_filtered_a <- cpm(all_counts_DGE_filtered, normalized.lib.sizes = TRUE)
#counts_filtered remains a DGEList object
#counts_normal_cpm is matrix of CPM

#add in gene names
diploid_kd_all_cpm <- cbind(all_counts_DGE_filtered$gene, all_counts_DGE_filtered_a)
diploid_kd_all_cpm <- diploid_kd_all_cpm %>% dplyr::rename(gene = genes)

#strip leading/trailing spaces from gene name column
diploid_kd_all_cpm$gene <- trimws(diploid_kd_all_cpm$gene, which = c("both"))

#read in ncbi GTF to extract TSS
#calculate total exon length for X and autosomal genes
mdo_GTF <- read.table("/Users/.../MonDom_GTF_mMomDom1.pri/ncbi_dataset/data/GCF_027887165.1/genomic.gtf", header=FALSE, sep="\t")
tot_exons <- mdo_GTF %>% filter(V3 == "exon")
tot_exons_A <- mdo_GTF %>% filter(V3 == "exon") %>% filter(V1 == "NC_077227.1" | V1 == "NC_077228.1" | V1 == "NC_077229.1" | V1 == "NC_077230.1" | V1 == "NC_077231.1" | V1 == "NC_077232.1" | V1 == "NC_077233.1" | V1 == "NC_077234.1")
tot_exons_X <- mdo_GTF %>% filter(V3 == "exon") %>% filter(V1 == "NC_077235.1")
tot_ex_length <- sum(tot_exons$V5) - sum(tot_exons$V4)
tot_ex_length_A <- sum(tot_exons_A$V5) - sum(tot_exons_A$V4)
tot_ex_length_X <- sum(tot_exons_X$V5) - sum(tot_exons_X$V4)

tot_gene_number_A <- mdo_GTF %>% filter(V3 == "gene") %>% filter(V1 == "NC_077227.1" | V1 == "NC_077228.1" | V1 == "NC_077229.1" | V1 == "NC_077230.1" | V1 == "NC_077231.1" | V1 == "NC_077232.1" | V1 == "NC_077233.1" | V1 == "NC_077234.1")
tot_gene_number_X <- mdo_GTF %>% filter(V3 == "gene") %>% filter(V1 == "NC_077235.1")
#total exon number 932,468
#total exon length 320,726,330
#1,024,473,833 = gene length
#autosomal_exon_length = 307,252,840
#X_exon_length  = 13,167,455

#filter for exons
ncbi_exon <- mdo_GTF %>% filter(V3=="exon")
ncbi_exon <- concat.split(ncbi_exon, 9, sep = ";")
ncbi_exon <- tibble(chromosome = ncbi_exon$V1, start = ncbi_exon$V4, end = ncbi_exon$V5, gene = ncbi_exon$V9_01)
ncbi_exon <- ncbi_exon %>% mutate(gene = sub("gene_id ","", gene)) %>% mutate("exonlength" = end - start)
ncbi_exon <- ncbi_exon %>% group_by(gene) %>% mutate(tot_exonlength = sum(exonlength))
ncbi_exon$gene <- trimws(ncbi_exon$gene, which = c("both"))
ncbi_exon2 <- ncbi_exon %>% dplyr::select(chromosome, gene, tot_exonlength) %>% distinct(chromosome, gene, .keep_all = TRUE)

#filter for gene entry
ncbi_gene <- mdo_GTF %>% filter(V3=="gene")
#gene number 32,665

#split out info column
ncbi_gene <- concat.split(ncbi_gene, 9, sep = ";")

#select columns
ncbi_gene <- tibble(chromosome = ncbi_gene$V1, start = ncbi_gene$V4, end = ncbi_gene$V5, strand = ncbi_gene$V7, gene = ncbi_gene$V9_01)

#remove text from column
ncbi_gene <- ncbi_gene %>% mutate(gene = sub("gene_id ","", gene)) %>% mutate("genelength" = end - start)

#filter gene list for unique entries (chose longest transcript for each)
ncbi_gene <- ncbi_gene %>% group_by(gene) %>% filter((end - start) == max(end - start))
#NB none lost with filtering

#strip leading/trailing spaces from gene name column
#space after gene name in ncbi_names file means no match.
ncbi_gene$gene <- trimws(ncbi_gene$gene, which = c("both"))
diploid_kd_all_cpm2 <- left_join(diploid_kd_all_cpm, ncbi_gene, keep=FALSE)

#reorder tibble
diploid_kd_all_cpm2 <- diploid_kd_all_cpm2 %>% relocate(chromosome, .before=gene) %>% relocate(start, .after=chromosome) %>% relocate(end, .after=start)

#add in exon length
diploid_kd_all_cpm_exon <- left_join(diploid_kd_all_cpm2, ncbi_exon2, by = c("gene"), keep=FALSE)  %>% dplyr::select(-chromosome.y) %>% dplyr::rename("chromosome" = "chromosome.x")

#rename chromosomes
#annotate unanchored genes with "zz"
diploid_kd_all_cpm2 <- diploid_kd_all_cpm2 %>% mutate(chromosome = sub("NC_077227.1","1", chromosome)) %>%
  mutate(chromosome = sub("NC_077228.1","2", chromosome)) %>%
  mutate(chromosome = sub("NC_077229.1","3", chromosome)) %>%
  mutate(chromosome = sub("NC_077230.1","4", chromosome)) %>%
  mutate(chromosome = sub("NC_077231.1","5", chromosome)) %>%
  mutate(chromosome = sub("NC_077232.1","6", chromosome)) %>%
  mutate(chromosome = sub("NC_077233.1","7", chromosome)) %>%
  mutate(chromosome = sub("NC_077234.1","8", chromosome)) %>%
  mutate(chromosome = sub("NC_077235.1","X", chromosome)) %>%
  mutate(chromosome = replace(chromosome, str_detect(chromosome, "NW_"),"zz")) %>%
  mutate(chromosome = sub("NC_077236.1","Y", chromosome)) 


#FILTER FOR EXPRESSION >1 CPM IN AT LEAST ONE KD OR CONTROL
diploid_kd_all_cpm2 <- diploid_kd_all_cpm2 %>% dplyr::filter(CAPRIN1 > 1 | CKAP4 > 1 | HNRNPK > 1 | NCL > 1 | NONO > 1 | SYNCRIP > 1 | U6M2 > 1)

diploid_kd_all_cpm2_A <- diploid_kd_all_cpm2 %>% filter(XA == "A")
diploid_kd_all_cpm2_X <- diploid_kd_all_cpm2 %>% filter(XA == "X")

#Add expression ratio relative to U6M2 control
a <- diploid_kd_all_cpm2 %>% mutate("expression_ratio" = log2(HNRNPK/U6M2)) %>% mutate("type" = "U6M2") %>% mutate("kd" = "HNRNPK") %>% dplyr::select("chromosome", "start", "end", "gene", "expression_ratio","type", "kd") 
c <- diploid_kd_all_cpm2 %>% mutate("expression_ratio" = log2(CKAP4/U6M2)) %>% mutate("type" = "U6M2") %>% mutate("kd" = "CKAP4") %>% dplyr::select("chromosome", "start", "end", "gene", "expression_ratio","type", "kd") 
e <- diploid_kd_all_cpm2 %>% mutate("expression_ratio" = log2(SYNCRIP/U6M2)) %>% mutate("type" = "U6M2") %>% mutate("kd" = "SYNCRIP") %>% dplyr::select("chromosome", "start", "end", "gene", "expression_ratio","type", "kd") 
g <- diploid_kd_all_cpm2 %>% mutate("expression_ratio" = log2(CAPRIN1/U6M2)) %>% mutate("type" = "U6M2") %>% mutate("kd" = "CAPRIN1") %>% dplyr::select("chromosome", "start", "end", "gene", "expression_ratio","type", "kd") 
i <- diploid_kd_all_cpm2 %>% mutate("expression_ratio" = log2(NONO/U6M2)) %>% mutate("type" = "U6M2") %>% mutate("kd" = "NONO") %>% dplyr::select("chromosome", "start", "end", "gene", "expression_ratio","type", "kd") 
k <- diploid_kd_all_cpm2 %>% mutate("expression_ratio" = log2(NCL/U6M2)) %>% mutate("type" = "U6M2") %>% mutate("kd" = "NCL") %>% dplyr::select("chromosome", "start", "end", "gene", "expression_ratio","type", "kd") 
m <- diploid_kd_all_cpm2 %>% mutate("expression_ratio" = log2(NUDE/U6M2)) %>% mutate("type" = "U6M2") %>% mutate("kd" = "zNUDE") %>% dplyr::select("chromosome", "start", "end", "gene", "expression_ratio","type", "kd") 
dip_kd_u6m2 <- rbind(a,c,e,g,i,k,m)
dip_kd_u6m2 <- as_tibble(dip_kd_u6m2)

#add in column for A v X - filter out unplaced genes (chromosome == "zz")
dip_kd_u6m2 <- dip_kd_u6m2 %>% filter(chromosome != "zz" & chromosome != "Y") %>% mutate(chr = if_else(chromosome == "X", "X", "A"))

#compile tibbles of counts
a <- diploid_kd_all_cpm2 %>% mutate("normalised_cpm" = HNRNPK) %>% mutate("kd" = "HNRNPK") %>% select("chromosome", "start", "end", "gene", "normalised_cpm","kd") 
b <- diploid_kd_all_cpm2 %>% mutate("normalised_cpm" = CKAP4) %>% mutate("kd" = "CKAP4") %>% select("chromosome", "start", "end", "gene", "normalised_cpm","kd") 
c <- diploid_kd_all_cpm2 %>% mutate("normalised_cpm" = SYNCRIP) %>% mutate("kd" = "SYNCRIP") %>% select("chromosome", "start", "end", "gene", "normalised_cpm","kd") 
d <- diploid_kd_all_cpm2 %>% mutate("normalised_cpm" = CAPRIN1) %>% mutate("kd" = "CAPRIN1") %>% select("chromosome", "start", "end", "gene", "normalised_cpm","kd") 
e <- diploid_kd_all_cpm2 %>% mutate("normalised_cpm" = NONO) %>% mutate("kd" = "NONO") %>% select("chromosome", "start", "end", "gene", "normalised_cpm","kd") 
f <- diploid_kd_all_cpm2 %>% mutate("normalised_cpm" = NCL) %>% mutate("kd" = "NCL") %>% select("chromosome", "start", "end", "gene", "normalised_cpm","kd") 
g <- diploid_kd_all_cpm2 %>% mutate("normalised_cpm" = U6M2) %>% mutate("kd" = "U6M2") %>% select("chromosome", "start", "end", "gene", "normalised_cpm","kd") 
h <- diploid_kd_all_cpm2 %>% mutate("normalised_cpm" = NUDE) %>% mutate("kd" = "NUDE") %>% select("chromosome", "start", "end", "gene", "normalised_cpm","kd") 

dip_kd_u6m2a <- rbind(a,b,c,d,e,f,g,h)
dip_kd_u6m2a <- as_tibble(dip_kd_u6m2a)
dip_kd_u6m2a <- dip_kd_u6m2a %>% filter (chromosome != "zz" & chromosome != "Y") %>% mutate(chr = if_else(chromosome == "X", "X", "A"))

tallies <- dip_kd_u6m2a %>% count(kd,chromosome)
dc <- tallies %>% filter(chromosome != "X")


#**FIGURE2
#CHECK knockdowns expression level cf RT-qPCR (Figure)
setwd("/Users/.../MDO_diploid_kd_normalised_counts_fc_mMonDom1/")

diploid_kd_all_cpm2 <- read.table("/Users/.../diploid_kd_all_cpm2_subread_mMonDom1_1cpm.txt", header=TRUE, sep="\t")

diploid_kd_all_cpm2 <- as_tibble(diploid_kd_all_cpm2) 

#replace NONO gene name (renamed in MonDom1 cf MonDom5)
diploid_kd_all_cpm2 <- diploid_kd_all_cpm2 %>% mutate(gene = sub("LOC100013374", "NONO", gene))

kds_fc <- diploid_kd_all_cpm2 %>% filter(gene == "CAPRIN1" | gene == "CKAP4" | gene == "HNRNPK" | gene == "NCL" | gene == "NONO" | gene == "SYNCRIP") %>%
  select(-chromosome) %>% select(-start) %>% select(-end) %>% select(-strand) %>% arrange(gene)

write.table(kds_fc, file="/Users/.../knockdowns_fc_mMonDom1_cpm1.txt", quote=FALSE, sep="\t", row.names = FALSE)

kds <- c("CAPRIN1", "CKAP4", "HNRNPK", "NCL", "NONO", "SYNCRIP")

#calculate proportions in order
exp_prop_fc <- as.numeric(c(kds_fc[1,2]/kds_fc[1,8], kds_fc[2,3]/kds_fc[2,8], kds_fc[3,4]/kds_fc[3,8], kds_fc[4,5]/kds_fc[4,8], kds_fc[5,6]/kds_fc[5,8], kds_fc[6,7]/kds_fc[6,8]))
exp_prop_fc
#0.0.6121682 0.7422934 0.7931676 0.5621149 0.8306200 0.6097362 0.6097362

#input RT-qPCR results
qpcr <- as.numeric(c(0.2929, 0.3502, 0.7673, 0.3802, 0.4647, 0.4102))

kd_table1 <- tibble(target = kds, exp_prop = exp_prop_fc) %>% mutate(type = "RNA-seq")
kd_table2 <- tibble(target = kds, exp_prop = qpcr) %>% mutate(type = "qPCR")

#PLOT separately to add in error bars for RT-qPCR:
kd_plot_qpcr <- ggplot(kd_table2, aes(x = target, y = exp_prop))+
  geom_col(fill = "skyblue3")+
  geom_errorbar(aes(x = target, ymin = c(0.271, 0.292, 0.726, 0.360, 0.437, 0.375),
                    ymax = c(0.316, 0.401, 0.811, 0.420, 0.495, 0.449)), width = 0.4, alpha = 0.9)+
  geom_hline(yintercept = 1, size=0.4, linetype=1)+
  xlab("")+
  ylab("Expression relative to U6M2 control")+
  theme_classic(base_size = 18) +
  theme(legend.position="none")
kd_plot_qpcr
ggsave("colplot_fc_mMonDom1_qPCR_cpm1.svg", width = 8, height = 6, units = "in")

kd_plot_rnaseq <- ggplot(kd_table1, aes(x = target, y = exp_prop))+
  geom_col(fill = "orange1")+
  geom_hline(yintercept = 1, size=0.4, linetype=1)+
  xlab("")+
  ylab("Expression relative to U6M2 control")+
  theme_classic(base_size = 18) +
  theme(legend.position="none")
kd_plot_rnaseq
ggsave("colplot_fc_mMonDom1_RNAseq_cpm1.svg", width = 8, height = 6, units = "in")

#Calculate p values for RT-qPCR - use Wilcoxon two sample - no assumption of Gaussian distribution (as required for t-test)
setwd("/Users/.../kd_effect/")
qpcr <- read.table("RT-qPCR_pvalue_calcs.txt",header=TRUE, sep="\t")

#Mann-Whitney U test - equivalent to Kruskall-Wallis, but for 2 groups
#The two-sample Mann–Whitney U test is a rank-based test that compares values for two groups.  A significant result suggests that the values for the two groups are different.  
#It is equivalent to a two-sample Wilcoxon rank-sum test.

results_WX <- data.frame(id = character(0), p.value = numeric(0))
id = unique(qpcr$target_kd)

for(j in id){
  WX <- rstatix::wilcox_test(data = qpcr %>% filter(target_kd == j), deltaCT_kd ~ cat, paired = FALSE, p.adjust.method = "holm")
  p.value <- WX$p
  res <- c(j, p.value)
  results_WX[nrow(results_WX)+1,] <- res
}  
results_WX <- as_tibble(results_WX)

setwd("/Users/.../MDO_diploid_kd_normalised_counts_fc_mMonDom1/")
write.table(results_WX, file="RT-qPCR_WX_test_pvalue.txt", quote=FALSE, sep="\t", row.names = FALSE)
# p=0.1 for all


#**FIGURE 3
#Boxplot of log2(kd/u6m2) for all chromosomes - autosomes combined
boxplot_c_htseq <- ggplot(dip_kd_u6m2 %>% filter (type == "U6M2"), aes(kd, expression_ratio))+
  geom_boxplot(notch = TRUE, aes(fill = chr), outlier.shape=NA) +
  scale_fill_manual(values = c("A" = "skyblue2", "X" = "orange2"),
                    labels = c("A" = "Autosomes", "X chromosome"))+
  # labs(title = "X chromosome") +
  geom_hline(yintercept = 0, size=0.4, linetype=1)+
  xlab("")+
  ylab("log2 expression ratio (knockdown/U6M2 control)")+
  theme_classic(base_size = 16) +
  scale_y_continuous(breaks = c(-0.4, 0, 0.4), labels = c("-0.4", "0", "0.4"))+
  coord_cartesian(ylim = c(-0.8, 0.8))+
  theme(legend.title = element_blank())
boxplot_c_htseq  
ggsave(file = "boxplot_fc_mMonDom1_cpm1.svg",
       plot = boxplot_c_htseq, device = "svg", width = 12, height = 8, units = "in")

#count gene number for plot
dat <- dip_kd_u6m2 %>% filter (type == "U6M2") %>% count(chr, kd)
#autosomal = 13,487
#X = 369

#LOOP to calculate 
setwd("/Users/.../MDO_diploid_kd_normalised_counts_fc_mMonDom1/")
#Wilcox test (non-paired) - test whether median differs from zero
#Grouped autosomes
results_gr <- data.frame(cat = character(0), kd = character(0), medi = character(0), p.value = numeric(0))
chr = unique(dip_kd_u6m2$chr)
kd = unique(dip_kd_u6m2$kd)

for(i in chr){
  for(j in kd){
      data <- dip_kd_u6m2 %>% dplyr::filter(type == "U6M2") %>% dplyr::filter(chr == i) %>% dplyr::filter(kd == j) 
      data2 <- as.vector(data$expression_ratio)
      data2 <- unlist(data2)
      wilcox <- wilcox.test(data2, alternative = "two.sided", mu = 0, paired = FALSE)
      p.value <- wilcox$p.value
      medi <- median(data2)
      res <- c(i, j, medi, p.value)
      results_gr[nrow(results_gr)+1,] <- res
    }
  } 


results_gr <- as_tibble(results_gr)
results_gr$p.value <- as.numeric(results_gr$p.value)
results_gr$medi <- as.numeric(results_gr$medi)
results_gr <- results_gr %>% mutate("p.signif" = case_when(p.value <= 0.0001 ~ "****", p.value <= 0.001 ~ "***", p.value <= 0.01 ~ "**", p.value <= 0.05 ~ "*", p.value > 0.05 ~ "NS"))
results_gr$p.value <- format(results_gr$p.value, scientific = TRUE, digits = 3)
results_gr$medi <- format(results_gr$medi, scientific = TRUE, digits = 3)

write.table(results_gr, file="Wilcox_test_fc_mMonDom1_cpm1.txt", quote=FALSE, sep="\t", row.names = FALSE)


#MOOD'S MEDIAN TEST - TEST WHETHER INDEPENDENT SAMPLES HAVE DIFFERENT MEDIANS
setwd("/Users/.../MDO_diploid_kd_normalised_counts_fc_mMonDom1/")
results_moods <- data.frame(kd = character(0), p.value = numeric(0))
kd = unique(dip_kd_u6m2$kd)

    for(k in kd){
      moods <- rcompanion::pairwiseMedianTest(expression_ratio ~ chr, data = dip_kd_u6m2 %>% dplyr::filter(kd == k) %>% dplyr::filter(expression_ratio != "-Inf"))
      p.value <- moods$p.value
      res <- c(k, p.value)
      results_moods[nrow(results_moods)+1,] <- res
    }


results_moods <- as_tibble(results_moods)
results_moods$p.value <- as.numeric(results_moods$p.value)
results_moods <- results_moods %>% mutate("p.signif" = case_when(p.value <= 0.0001 ~ "****", p.value <= 0.001 ~ "***", p.value <= 0.01 ~ "**", p.value <= 0.05 ~ "*", p.value > 0.05 ~ "NS"))
results_moods$p.value <- format(results_moods$p.value, scientific = TRUE, digits = 3)

write.table(results_moods, file="Moods_median_fc_mMonDom1_cpm1.txt", quote=FALSE, sep="\t", row.names = FALSE)

#WILCOX/MANN WHITNEY U TEST - TEST WHETHER INDEPENDENT SAMPLES DIFFERENT
#NB: stats::wilcox.test does not do adjustment - prefer to use rstatix::wilcox_test (Holm adjustment)
results_mw <- data.frame(kd = character(0), p.value = numeric(0))
kd = unique(dip_kd_u6m2$kd)

for(k in kd){
  mw <- rstatix::wilcox_test(data = dip_kd_u6m2 %>% dplyr::filter(kd == k) %>% dplyr::filter(expression_ratio != "-Inf"), expression_ratio ~ chr)
  p.value <- mw$p
  res <- c(k, p.value)
  results_mw[nrow(results_mw)+1,] <- res
}

results_mw <- as_tibble(results_mw)
results_mw$p.value <- as.numeric(results_mw$p.value)
results_mw <- results_mw %>% mutate("p.signif" = case_when(p.value <= 0.0001 ~ "****", p.value <= 0.001 ~ "***", p.value <= 0.01 ~ "**", p.value <= 0.05 ~ "*", p.value > 0.05 ~ "NS"))
results_mw$p.value <- format(results_mw$p.value, scientific = TRUE, digits = 3)

write.table(results_mw, file="Mann_WhitneyU_fc_mMonDom1_cpm1.txt", quote=FALSE, sep="\t", row.names = FALSE)


#FIGURE 3B - X:MEDIAN AUTOSOMAL
# upregulation/downregulation - normalised by autosomal medians
#X to autosomal expression ratio
#convert cpm to fpkm to adjust for different gene lengths (use fpkm for paired end data rather than rpkm)
#FPKM = (CPM x 1000) / length (bp)

setwd("/Users/.../MDO_diploid_kd_normalised_counts_fc_mMonDom1/")

diploid_kd_all_cpm_exon <- read.table("/Users/diploid_kd_all_cpm2_subread_mMonDom1_1cpm_exon.txt", header=TRUE, sep="\t")

diploid_kd_all_cpm_exon <- diploid_kd_all_cpm_exon %>% mutate(chromosome = sub("NC_077227.1","1", chromosome)) %>%
  mutate(chromosome = sub("NC_077228.1","2", chromosome)) %>%
  mutate(chromosome = sub("NC_077229.1","3", chromosome)) %>%
  mutate(chromosome = sub("NC_077230.1","4", chromosome)) %>%
  mutate(chromosome = sub("NC_077231.1","5", chromosome)) %>%
  mutate(chromosome = sub("NC_077232.1","6", chromosome)) %>%
  mutate(chromosome = sub("NC_077233.1","7", chromosome)) %>%
  mutate(chromosome = sub("NC_077234.1","8", chromosome)) %>%
  mutate(chromosome = sub("NC_077235.1","X", chromosome)) %>%
  mutate(chromosome = replace(chromosome, str_detect(chromosome, "NW_"),"zz")) %>%
  mutate(chromosome = sub("NC_077236.1","Y", chromosome)) 

diploid_kd_all_cpm_exon <- diploid_kd_all_cpm_exon %>% dplyr::filter(chromosome != "Y" & chromosome != "zz") %>% mutate(XA = if_else(chromosome =="X", "X", "A"), .after = chromosome)

diploid_kd_all_cpm_exon2 <- diploid_kd_all_cpm_exon %>% dplyr::filter(CAPRIN1 > 1 | CKAP4 > 1 | HNRNPK > 1 | NCL > 1 | NONO > 1 | SYNCRIP > 1 | U6M2 > 1)

diploid_kd_all_fpkm <- diploid_kd_all_cpm_exon2 %>% mutate(CAPRIN1_fpkm = (CAPRIN1*1000)/tot_exonlength) %>% mutate(CKAP4_fpkm = (CKAP4*1000)/tot_exonlength) %>%
  mutate(HNRNPK_fpkm = (HNRNPK*1000)/tot_exonlength) %>% mutate(NCL_fpkm = (NCL*1000)/tot_exonlength) %>% mutate(NONO_fpkm = (NONO*1000)/tot_exonlength) %>%
  mutate(SYNCRIP_fpkm = (SYNCRIP*1000)/tot_exonlength) %>% mutate(U6M2_fpkm = (U6M2*1000)/tot_exonlength) %>% mutate(NUDE_fpkm = (NUDE*1000)/tot_exonlength)

diploid_kd_all_fpkmX <- diploid_kd_all_fpkm %>% filter(XA == "X")
diploid_kd_all_fpkmA <- diploid_kd_all_fpkm %>% filter(XA == "A")

#MEDIAN AUTOSOMAL EXPRESSION IN EACH SAMPLE
dip_auto <- diploid_kd_all_fpkm %>% dplyr::filter(XA == "A") 

dip_med_cap <- median(dip_auto[,17])
dip_med_ckap <- median(dip_auto[,18])
dip_med_hk <- median(dip_auto[,19])
dip_med_ncl <- median(dip_auto[,20])
dip_med_nono <- median(dip_auto[,21])
dip_med_syn <- median(dip_auto[,22])
dip_med_u6m2 <- median(dip_auto[,23])
dip_med_nude <- median(dip_auto[,24])

#MEDIAN X EXPRESSION IN EACH SAMPLE
dip_autoX <- diploid_kd_all_fpkm %>% dplyr::filter(XA == "X") 

dip_med_capX <- median(dip_autoX[,17])
dip_med_ckapX <- median(dip_autoX[,18])
dip_med_hkX <- median(dip_autoX[,19])
dip_med_nclX <- median(dip_autoX[,20])
dip_med_nonoX <- median(dip_autoX[,21])
dip_med_synX <- median(dip_autoX[,22])
dip_med_u6m2X <- median(dip_autoX[,23])
dip_med_nudex <- median(dip_autoX[,24])

#NORMALISE EXPRESSION BY AUTOSOMAL MEDIAN FOR EACH SAMPLE
dip_kd_norm <- diploid_kd_all_fpkm %>% mutate(CAPRIN1_fpkm = CAPRIN1_fpkm/dip_med_cap) %>% mutate(CKAP4_fpkm = CKAP4_fpkm/dip_med_ckap) %>% mutate(HNRNPK_fpkm = HNRNPK_fpkm/dip_med_hk) %>%
  mutate(NCL_fpkm = NCL_fpkm/dip_med_ncl) %>% mutate(NONO_fpkm = NONO_fpkm/dip_med_nono) %>% mutate(SYNCRIP_fpkm = SYNCRIP_fpkm/dip_med_syn) %>% mutate(U6M2_fpkm = U6M2_fpkm/dip_med_u6m2) %>% mutate(NUDE_fpkm = NUDE_fpkm/dip_med_nude)

a <- dip_kd_norm %>% dplyr::select(chromosome, XA, gene, "fpkm" = CAPRIN1_fpkm) %>% mutate(id = "CAPRIN1")
b <- dip_kd_norm %>% dplyr::select(chromosome, XA, gene, "fpkm" = CKAP4_fpkm) %>% mutate(id = "CKAP4")
c <- dip_kd_norm %>% dplyr::select(chromosome, XA, gene, "fpkm" = HNRNPK_fpkm) %>% mutate(id = "HNRNPK")
d <- dip_kd_norm %>% dplyr::select(chromosome, XA, gene, "fpkm" = NCL_fpkm) %>% mutate(id = "NCL")
e <- dip_kd_norm %>% dplyr::select(chromosome, XA, gene, "fpkm" = NONO_fpkm) %>% mutate(id = "NONO")
f <- dip_kd_norm %>% dplyr::select(chromosome, XA, gene, "fpkm" = SYNCRIP_fpkm) %>% mutate(id = "SYNCRIP")
g <- dip_kd_norm %>% dplyr::select(chromosome, XA, gene, "fpkm" = U6M2_fpkm) %>% mutate(id = "U6M2")
h <- dip_kd_norm %>% dplyr::select(chromosome, XA, gene, "fpkm" = NUDE_fpkm) %>% mutate(id = "zNO_VECTOR")

all_norm <- rbind(a, b, c, d, e, f, g, h)

#NB median for NUDE is not zero due to genes not expressed in Nude sample only (ie log2 = undefined, so median different)
all_norm_plot <- ggplot(all_norm, aes(x = id, y = fpkm, fill = XA))+
  geom_boxplot(notch = TRUE, outliers = FALSE)+
  scale_fill_manual(values = c("A" = "skyblue2", "X" = "orange2"),
                    labels = c("A" = "Autosomes", "X chromosome"))+
  #geom_hline(yintercept = 0, size=0.4, linetype=1)+
  ylab("FPKM expression ratio (gene / autosomal median)")+
  xlab("")+
  theme_classic(base_size = 16)+
  stat_summary(data = all_norm %>% filter(XA == "A"), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = 2, vjust = -1))+
  stat_summary(data = all_norm %>% filter(XA == "X"), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = -0.5, vjust = -1))+
  theme(legend.title = element_blank())
all_norm_plot

ggsave(file = "FPKM_boxplot_norm.svg", plot = all_norm_plot,  device = "svg", width = 12, height = 8, units = "in")

#MOODS MEDIAN
all_norm2 <- all_norm %>% filter(fpkm != 0) %>% mutate(log2_fpkm = log2(fpkm))
results_moods_fpkm <- data.frame(id = character(0), p.value = numeric(0))
id = unique(all_norm2$id)

for(k in id){
  moods <- rcompanion::pairwiseMedianTest(fpkm ~ XA, data = all_norm2 %>% dplyr::filter(id == k))
  p.value <- moods$p.value
  res <- c(k, p.value)
  results_moods_fpkm[nrow(results_moods_fpkm)+1,] <- res
}


results_moods_fpkm <- as_tibble(results_moods_fpkm)
results_moods_fpkm$p.value <- as.numeric(results_moods_fpkm$p.value)
results_moods_fpkm <- results_moods_fpkm %>% mutate("p.signif" = case_when(p.value <= 0.0001 ~ "****", p.value <= 0.001 ~ "***", p.value <= 0.01 ~ "**", p.value <= 0.05 ~ "*", p.value > 0.05 ~ "NS"))
results_moods_fpkm$p.value <- format(results_moods_fpkm$p.value, scientific = TRUE, digits = 3)

write.table(results_moods_fpkm, file="Moods_median_fc_mMonDom1_cpm1_fpkm.txt", quote=FALSE, sep="\t", row.names = FALSE)

#WILCOX/MANN WHITNEY U TEST - TEST WHETHER INDEPENDENT SAMPLES DIFFERENT
#NB: stats::wilcox.test does not seem to do adjustment - prefer to use rstatix::wilcox_test (Holm adjustment)
results_mw <- data.frame(kd = character(0), p.value = numeric(0))
kd = c("CAPRIN1", "CKAP4", "HNRNPK", "NCL", "NONO", "SYNCRIP", "zNUDE")

for(k in kd){
  mw <- rstatix::wilcox_test(data = all_norm2 %>% dplyr::filter(id == k | id == "U6M2"), fpkm ~ id)
  p.value <- mw$p
  res <- c(k, p.value)
  results_mw[nrow(results_mw)+1,] <- res
}


results_mw <- as_tibble(results_mw)
results_mw$p.value <- as.numeric(results_mw$p.value)
results_mw <- results_mw %>% mutate("p.signif" = case_when(p.value <= 0.0001 ~ "****", p.value <= 0.001 ~ "***", p.value <= 0.01 ~ "**", p.value <= 0.05 ~ "*", p.value > 0.05 ~ "NS"))
results_mw$p.value <- format(results_mw$p.value, scientific = TRUE, digits = 3)

write.table(results_mw, file="Mann_WhitneyU_fc_mMonDom1_cpm1_fpkm.txt", quote=FALSE, sep="\t", row.names = FALSE)

#COMPARE X ACROSS ALL SAMPLES
#USE PAIRWISE MOODS MEDIAN TEST
moods_cap <- rcompanion::pairwiseMedianTest(fpkm ~ id, data = all_norm2 %>% dplyr::filter(XA == "X") %>% dplyr::filter(id == "CAPRIN1" | id == "U6M2"))
moods_ckap <- rcompanion::pairwiseMedianTest(fpkm ~ id, data = all_norm2 %>% dplyr::filter(XA == "X") %>% dplyr::filter(id == "CKAP4" | id == "U6M2"))
moods_hk <- rcompanion::pairwiseMedianTest(fpkm ~ id, data = all_norm2 %>% dplyr::filter(XA == "X") %>% dplyr::filter(id == "HNRNPK" | id == "U6M2"))
moods_ncl <- rcompanion::pairwiseMedianTest(fpkm ~ id, data = all_norm2 %>% dplyr::filter(XA == "X") %>% dplyr::filter(id == "NCL" | id == "U6M2"))
moods_nono <- rcompanion::pairwiseMedianTest(fpkm ~ id, data = all_norm2 %>% dplyr::filter(XA == "X") %>% dplyr::filter(id == "NONO" | id == "U6M2"))
moods_syncrip <- rcompanion::pairwiseMedianTest(fpkm ~ id, data = all_norm2 %>% dplyr::filter(XA == "X") %>% dplyr::filter(id == "SYNCRIP" | id == "U6M2"))
moods_nude <- rcompanion::pairwiseMedianTest(fpkm ~ id, data = all_norm2 %>% dplyr::filter(XA == "X") %>% dplyr::filter(id == "zNUDE" | id == "U6M2"))

#Kruskall-Wallis test
KW <- rstatix::kruskal_test(data = all_norm2 %>% dplyr::filter(XA == "X"), fpkm ~ id)
write.table(KW, file="KW_u6m2_FPKM_median_fc_mMonDom1.txt", quote=FALSE, sep="\t", row.names = FALSE)

#DUNN Test
DT <- rstatix::dunn_test(data = all_norm2 %>% dplyr::filter(XA == "X"), fpkm ~ id, p.adjust.method = "holm")
write.table(DT, file="DT_u6m2_FPKM_median_fc_mMonDom1.txt", quote=FALSE, sep="\t", row.names = FALSE)

#WILCOXON test
WX <- rstatix::wilcox_test(data = all_norm2 %>% dplyr::filter(XA == "X"), fpkm ~ id)
write.table(WX, file="WX_u6m2_FPKM_median_fc_mMonDom1.txt", quote=FALSE, sep="\t", row.names = FALSE)

#counts
counts_fpkm <- all_norm2 %>% count(XA, id)


#GSEA: functional enrichments + chromosomal locations of genes +/- 0.5 log2 up/down regulated for nude cf U6M2 – ie, what genes are up/down regulated with transfection

setwd("/Users/.../MDO_diploid_kd_normalised_counts_fc_mMonDom1/")

nude_genes <- dip_kd_u6m2 %>% filter(kd == "zNUDE")
nude_genesA <- nude_genes %>% filter(chromosome != "X") #for autosomal
nude_genesX <- nude_genes %>% filter(chromosome == "X") #for X

upgenes_list <- nude_genes %>% filter(expression_ratio > 0.5 | expression_ratio == "Inf") %>% arrange(desc(expression_ratio))#1406
downgenes_list <- nude_genes %>% filter(expression_ratio < -0.5 | expression_ratio == "-Inf") %>% arrange(expression_ratio)#1504
upgenes_listA <- nude_genesA %>% filter(expression_ratio > 0.5 | expression_ratio == "Inf") %>% arrange(desc(expression_ratio))#1360
downgenes_listA <- nude_genesA %>% filter(expression_ratio < -0.5 | expression_ratio == "-Inf") %>% arrange(expression_ratio)#1485
upgenes_listX <- nude_genesX %>% filter(expression_ratio > 0.5 | expression_ratio == "Inf") %>% arrange(desc(expression_ratio))#46
downgenes_listX <- nude_genesX %>% filter(expression_ratio < -0.5 | expression_ratio == "-Inf") %>% arrange(expression_ratio)#19

#GSEA using Gprofiler2 with gSCS multiple testing correction, ordered query
#For all combined
upgenes = gost(query = list("up" = upgenes_list$gene), organism="hsapiens",
               evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.05, ordered_query = TRUE,
               sources = c("GO:BP", "GO:MF", "GO:CC"))
Up = upgenes$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size", "source")]
colnames(Up) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size", "Source")
Up$FDR = Up$p.Val
Up$Phenotype = "1"
Up = Up[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Gene.number", "Genes", "Term.size", "Source")]
Up <- Up %>% filter(Term.size<4000)

write.table(Up, file = "gProfiler_Up_nude_SCS.txt", sep = "\t", quote = F, row.names = F)

downgenes = gost(query = list("down" = downgenes_list$gene), organism="hsapiens",
                 evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.05, ordered_query = TRUE,
                 sources = c("GO:BP", "GO:MF", "GO:CC"))
Down = downgenes$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size", "source")]
colnames(Down) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size", "Source")
Down$FDR = Down$p.Val
Down$Phenotype = "1"
Down = Down[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Gene.number", "Genes", "Term.size", "Source")]
Down <- Down %>% filter(Term.size<4000)

write.table(Down, file = "gProfiler_Down_nude_SCS.txt", sep = "\t", quote = F, row.names = F)

#For autosomal only and X only - separately
upgenes = gost(query = list("up" = upgenes_listA$gene), organism="hsapiens",
                          evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.05, ordered_query = TRUE,
                          sources = c("GO:BP", "GO:MF", "GO:CC"))
Up = upgenes$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size", "source")]
colnames(Up) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size", "Source")
Up$FDR = Up$p.Val
Up$Phenotype = "1"
Up = Up[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Gene.number", "Genes", "Term.size", "Source")]
Up <- Up %>% filter(Term.size<4000)

write.table(Up, file = "gProfiler_Up_nudeA_SCS.txt", sep = "\t", quote = F, row.names = F)


upgenesX = gost(query = list("up" = upgenes_listX$gene), organism="hsapiens",
                  evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.05, ordered_query = TRUE,
                  sources = c("GO:BP", "GO:MF", "GO:CC"))
upX = upgenesX$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size", "source")]
colnames(upX) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size", "Source")
upX$FDR = upX$p.Val
upX = upX[,c("GO.ID", "Description", "p.Val", "FDR", "Gene.number", "Genes", "Term.size", "Source")]
write.table(upX, file = "gProfiler_Up_nudeX_SCS.txt", sep = "\t", quote = F, row.names = F)


downgenes = gost(query = list("down" = downgenes_listA$gene), organism="hsapiens",
               evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.05, ordered_query = TRUE,
               sources = c("GO:BP", "GO:MF", "GO:CC"))
Down = downgenes$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size", "source")]
colnames(Down) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size", "Source")
Down$FDR = Down$p.Val
Down$Phenotype = "1"
Down = Down[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Gene.number", "Genes", "Term.size", "Source")]
Down <- Down %>% filter(Term.size<4000)

write.table(Down, file = "gProfiler_Down_nudeA_SCS.txt", sep = "\t", quote = F, row.names = F)

downgenesX = gost(query = list("down" = downgenes_listX$gene), organism="hsapiens",
                 evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.05, ordered_query = TRUE,
                 sources = c("GO:BP", "GO:MF", "GO:CC"))
DownX = downgenesX$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size", "source")]
colnames(DownX) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size", "Source")
DownX$FDR = DownX$p.Val
DownX = DownX[,c("GO.ID", "Description", "p.Val", "FDR", "Gene.number", "Genes", "Term.size", "Source")]
write.table(DownX, file = "gProfiler_Down_nudeX_SCS.txt", sep = "\t", quote = F, row.names = F)

#GSEA of AUTOSOMAL genes UPregulated on RNAi knockdown. Writes tables for input to Cytoscape Enrichment map
#Figure 13
A_dip_caprin <- dip_kd_u6m2 %>% filter(kd == "CAPRIN1") %>% filter(chr == "A") %>% filter(expression_ratio > 0.5) %>% arrange(desc(expression_ratio))
A_dip_syncrip <- dip_kd_u6m2 %>% filter(kd == "SYNCRIP") %>% filter(chr == "A") %>% filter(expression_ratio > 0.5) %>% arrange(desc(expression_ratio))
A_dip_hnrnpk <- dip_kd_u6m2 %>% filter(kd == "HNRNPK") %>% filter(chr == "A") %>% filter(expression_ratio > 0.5) %>% arrange(desc(expression_ratio))
A_dip_ckap <- dip_kd_u6m2 %>% filter(kd == "CKAP4") %>% filter(chr == "A") %>% filter(expression_ratio > 0.5) %>% arrange(desc(expression_ratio))

#caprin
upgenes = gost(query = list("up" = A_dip_caprin$gene), organism="hsapiens",
               evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.05, ordered_query = TRUE,
               sources = c("GO:BP", "GO:MF", "GO:CC"))
Up = upgenes$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size", "source")]
colnames(Up) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size", "Source")
Up$FDR = Up$p.Val
Up$Phenotype = "1"
Up = Up[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Gene.number", "Genes", "Term.size", "Source")]
Up_caprin <- Up %>% filter(Term.size<=3000) #%>% slice_min(p.Val, n = 50)
write.table(Up_caprin, file = "gProfiler_A_up_caprin.txt", sep = "\t", quote = F, row.names = F)

#hnrnpk
upgenes = gost(query = list("up" = A_dip_hnrnpk$gene), organism="hsapiens",
               evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.05, ordered_query = TRUE,
               sources = c("GO:BP", "GO:MF", "GO:CC"))
Up = upgenes$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size", "source")]
colnames(Up) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size", "Source")
Up$FDR = Up$p.Val
Up$Phenotype = "1"
Up = Up[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Gene.number", "Genes", "Term.size", "Source")]
Up_hnrnpk <- Up %>% filter(Term.size<=3000) #%>% slice_min(p.Val, n = 50)
write.table(Up_hnrnpk, file = "gProfiler_A_up_hnrnpk.txt", sep = "\t", quote = F, row.names = F)

#syncrip
upgenes = gost(query = list("up" = A_dip_syncrip$gene), organism="hsapiens",
               evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.05, ordered_query = TRUE,
               sources = c("GO:BP", "GO:MF", "GO:CC"))
Up = upgenes$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size", "source")]
colnames(Up) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size", "Source")
Up$FDR = Up$p.Val
Up$Phenotype = "1"
Up = Up[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Gene.number", "Genes", "Term.size", "Source")]
Up_syncrip <- Up %>% filter(Term.size<=3000) #%>% slice_min(p.Val, n = 50)
write.table(Up_syncrip, file = "gProfiler_A_up_syncrip.txt", sep = "\t", quote = F, row.names = F)

#CKAP4
#ckap
upgenes = gost(query = list("up" = A_dip_ckap$gene), organism="hsapiens",
               evcodes = TRUE, multi_query = FALSE, correction_method = "g_SCS", user_threshold = 0.05, ordered_query = TRUE,
               sources = c("GO:BP", "GO:MF", "GO:CC"))
Up = upgenes$result[,c("term_id", "term_name", "p_value", "intersection_size", "intersection", "term_size", "source")]
colnames(Up) = c("GO.ID", "Description", "p.Val", "Gene.number", "Genes", "Term.size", "Source")
Up$FDR = Up$p.Val
Up$Phenotype = "1"
Up = Up[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Gene.number", "Genes", "Term.size", "Source")]
Up_ckap <- Up %>% filter(Term.size<=3000) #%>% slice_min(p.Val, n = 50)
write.table(Up_ckap, file = "gProfiler_A_up_ckap.txt", sep = "\t", quote = F, row.names = F)

