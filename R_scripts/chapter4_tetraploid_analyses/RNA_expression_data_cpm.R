#this compares expression diploid v tetraploid in terms of log2(tetraploid_cpm/diploid_cpm)
#consider genes expressed foldchange different between diploid and tetraploid
#this uses MonDom5 GTF from ncbi

library("edgeR")
library(scales)
library(fabricatr) #use for split_quantile
library(tidyverse)
library(splitstackshape)
library(plotly)

setwd("/Users/.../R_files/")

#read in normalised htseq count files
htseq_tetra_P30 <- read.table("htseq_count_files/MDO_htseq_tetra_nude_P30_ncbi.txt", header=FALSE, sep="\t")
htseq_tetra_P14 <- read.table("htseq_count_files/MDO_htseq_tetra_nude_P14_ncbi.txt", header=FALSE, sep="\t")
htseq_diploid <- read.table("htseq_count_files/MDO_diploid_htseq_ncbi.txt", header=FALSE, sep="\t")

htseq_tetra_P30 <- htseq_tetra_P30[2:33749,]
htseq_tetra_P14 <- htseq_tetra_P14[2:33749,]
htseq_diploid <- htseq_diploid[2:33749,]

names(htseq_tetra_P30) <- c("gene", "tetra_P30")
names(htseq_tetra_P14) <- c("gene", "tetra_P14")
names(htseq_diploid) <- c("gene", "diploid")

#combine tetraploid and diploid 
#first check all names align
htseq_all <- left_join(htseq_tetra_P14, htseq_tetra_P30)
htseq_all <- left_join(htseq_all, htseq_diploid)

#convert to DGEList object
group = c("tetra_P14", "tetra_P30", "diploid")
htseq_all_DGE <- DGEList(counts=htseq_all[ ,2:4], group=group, genes=htseq_all[ ,1])

#normalise TMM (normalises to calculate effective library sizes). Does not change count values 
htseq_all_norm <- calcNormFactors(htseq_all_DGE)
head(htseq_all_norm)

keep <- rowSums(cpm(htseq_all_norm)>1) >= 1
summary(keep)
#Mode   FALSE    TRUE 
#logical   20342   13406  

htseq_all_filtered <- htseq_all_norm[keep, , keep.lib.sizes=FALSE]
head(htseq_all_filtered)

htseq_all_cpm <- cpm(htseq_all_filtered, normalized.lib.sizes = TRUE)
#counts_filtered remains a DGEList object
#counts_normal_cpm is matrix of CPM

#add in gene names
htseq_all_cpm <- data.frame("gene" = htseq_all_filtered$genes, "tetra_P14_cpm" = htseq_all_cpm[,1], "tetra_P30_cpm" = htseq_all_cpm[,2], "diploid_cpm" = htseq_all_cpm[,3])

#strip leading/trailing spaces from gene name column
htseq_all_cpm$gene <- trimws(htseq_all_cpm$gene, which = c("both"))

#read in ncbi GTF to extract TSS
ncbi_GTF <- read.table("htseq_count_files/MDO_GTF_ncbi.gtf", header=FALSE, sep="\t")

#filter for gene entry
ncbi_genes <- ncbi_GTF %>% filter(V3=="gene")

#split out info column
ncbi_genes <- concat.split(ncbi_genes, 9, sep = ";")

#select columns
ncbi_genes <- tibble(chromosome = ncbi_genes$V1, start = ncbi_genes$V4, end = ncbi_genes$V5, strand = ncbi_genes$V7, gene = ncbi_genes$V9_01)

#remove text from column
ncbi_genes <- ncbi_genes %>% mutate(gene = sub("gene_id ","", gene)) %>% mutate("genelength" = end - start)

#filter gene list for unique entries (chose longest transcript for each)
ncbi_genes <- ncbi_genes %>% group_by(gene) %>% filter((end - start) == max(end - start))

#strip leading/trailing spaces from gene name column
#space after gene name in ncbi_names file means no match.
ncbi_genes$gene <- trimws(ncbi_genes$gene, which = c("both"))
tetra_dip_cpm <- left_join(htseq_all_cpm, ncbi_genes, keep=FALSE)

#reorder
tetra_dip_cpm <- tetra_dip_cpm %>% relocate(chromosome, .before=genes) %>% relocate(start, .after=chromosome) %>% relocate(end, .after=start)

#rename chromosomes
#annotate unanchored genes with "zzz
tetra_dip_cpm <- tetra_dip_cpm %>% mutate(chromosome = sub("NC_008801.1","1", chromosome)) %>%
  mutate(chromosome = sub("NC_008802.1","2", chromosome)) %>%
  mutate(chromosome = sub("NC_008803.1","3", chromosome)) %>%
  mutate(chromosome = sub("NC_008804.1","4", chromosome)) %>%
  mutate(chromosome = sub("NC_008805.1","5", chromosome)) %>%
  mutate(chromosome = sub("NC_008806.1","6", chromosome)) %>%
  mutate(chromosome = sub("NC_008807.1","7", chromosome)) %>%
  mutate(chromosome = sub("NC_008808.1","8", chromosome)) %>%
  mutate(chromosome = sub("NC_008809.1","X", chromosome)) %>%
  mutate(chromosome = replace(chromosome, str_detect(chromosome, "NW_"),"zzz"))

#search human ortholog names [uniprot, ensembl], import
tetra_dip_cpm_genes <- tetra_dip_cpm %>% dplyr::select(genes)
write.table(tetra_dip_cpm_genes, file="outputs/expression_data/tetra_dip_cpm_genes.txt", quote=FALSE, sep="\t", row.names = FALSE)

#from ensembl 104 biomart
#delete FBXL3, ZFP36L2, MPP2, CFAP94,  - as these gene names already exist in MDO
MDO_gene_annot1 <- read.table("htseq_count_files/MDO_GENE_ANNOT1.txt", header=FALSE, sep="\t")
names(MDO_gene_annot1)[1] <- "genes"
names(MDO_gene_annot1)[2] <- "hum_orth"

MDO_gene_annot2 <- read.table("htseq_count_files/MDO_GENE_ANNOT2.txt", header=FALSE, sep="\t")
names(MDO_gene_annot2)[1] <- "genes"
names(MDO_gene_annot2)[2] <- "hum_orth"

#combine, remove duplication
MDO_gene_annot <- rbind(MDO_gene_annot1, MDO_gene_annot2) %>% distinct()

tetra_dip_cpm <- left_join(tetra_dip_cpm, MDO_gene_annot) 

#where human ortholog gene name differs from MDO, replace with human ortholog gene name
tetra_dip_cpm <-  tetra_dip_cpm %>% mutate(newgene = case_when(is.na(hum_orth) ~ genes,
                                                                  genes != hum_orth ~ hum_orth))

#this adds new column "newgene" with updated gene name (compatible with human orthologs)                                                                                                                 

#write with unanchored included
write.table(tetra_dip_cpm, file="outputs/expression_data/tetra_dip_cpm.txt", quote=FALSE, sep="\t", row.names = FALSE)
tetra_dip_cpm <- read.table("outputs/expression_data/tetra_dip_cpm.txt", header = TRUE, sep="\t")

#write lists of genes for Venn (use BioVenn online)
gene_dip_cpm <- tetra_dip_cpm %>% dplyr::select("cpm" = diploid_cpm, genes) %>% filter(cpm > 1) %>% dplyr::select(genes)
gene_P14_cpm <- tetra_dip_cpm %>% dplyr::select("cpm" = tetra_P14_cpm, genes) %>% filter(cpm > 1) %>% dplyr::select(genes)
gene_P30_cpm <- tetra_dip_cpm %>% dplyr::select("cpm" = tetra_P30_cpm, genes) %>% filter(cpm > 1) %>% dplyr::select(genes)

setwd("/Users/kimmcintyre/Library/CloudStorage/OneDrive-UNSW/Kim_Thesis/Chapter_5_tetraploid/R_files/outputs/venn_expression_genelists/")
write_delim(gene_dip_cpm, file = "gene_dip_cpm", delim= "\t")
write_delim(gene_P14_cpm, file = "gene_P14_cpm", delim= "\t")
write_delim(gene_P30_cpm, file = "gene_P30_cpm", delim= "\t")



