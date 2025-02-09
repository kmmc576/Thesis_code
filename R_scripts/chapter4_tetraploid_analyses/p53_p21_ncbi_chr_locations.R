
library(tidyverse)
library(dplyr)
library(readr)
p53_list <- read.delim2("/Users/.../R_files/outputs/TP53_p21_locations/p53_gene_result.txt", header=TRUE, sep="\t")
p21_list <- read.delim2("/Users/.../R_files/outputs/TP53_p21_locations/p21_gene_result.txt", header=TRUE, sep="\t")

p53_list <- p53_list %>% dplyr::select("tax_id", "Org_name", "GeneID", "Symbol", "chromosome", "exon_count") 
colnames(p53_list) <- c("tax_id", "Org_name", "GeneID", "gene", "chromosome", "exon_count")

p53_list <- dplyr::filter(p53_list, gene=="TP53"|gene=="tp53") 
p53_list <- dplyr::filter(p53_list, chromosome!="Un") 

p21_list <- p21_list %>% dplyr::select("tax_id", "Org_name", "GeneID", "Symbol", "chromosome", "exon_count") 
colnames(p21_list) <- c("tax_id", "Org_name", "GeneID", "gene", "chromosome", "exon_count")

p21_list <- dplyr::filter(p21_list, gene == "CDKN1A" | gene == "cdkn1a" | gene == "Cdkn1a") 
p21_list <- dplyr::filter(p21_list, chromosome!="Un") 

joint <- full_join(p53_list, p21_list, by = "tax_id", keep = TRUE)
joint <- joint %>% select(tax_id.x, Org_name.x, gene.x, gene.y, chromosome.x, chromosome.y)
joint <- joint %>% filter(!is.na(gene.x)) %>% filter(!is.na(gene.y))
joint <- joint %>% mutate("match" = if_else(chromosome.x == chromosome.y, "true", "false"))

write_delim(joint, file="/Users/.../p53_p21_chromosomal_locations.txt", sep = "\t")

