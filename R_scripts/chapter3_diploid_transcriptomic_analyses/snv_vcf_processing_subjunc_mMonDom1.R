#import VCF files from bcftools mpileup, mcall
#vcf codes
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of raw reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=MQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Mapping Quality Bias (closer to 0 is better)">
##INFO=<ID=BQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Base Quality Bias (closer to 0 is better)">
##INFO=<ID=MQSBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Mapping Quality vs Strand Bias (closer to 0 is better)">
##INFO=<ID=SCBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Soft-Clip Length Bias (closer to 0 is better)">
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)">
##INFO=<ID=AD,Number=R,Type=Integer,Description="Total allelic depths (high-quality bases)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">
##INFO=<ID=AF2,Number=1,Type=Float,Description="Max-likelihood estimate of the first and second group ALT allele frequency (assuming HWE)">
##INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">
##INFO=<ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same">
##INFO=<ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">
##INFO=<ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies">
##INFO=<ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
#SAMPLE infor fields: GT:PL:DP:AD

library(vcfR)
library(splitstackshape)
library(tidyverse)
library(fuzzyjoin)
library(GenomicRanges)
library(dplyr)
library(reshape2)
library(ggpubr) #stat_cor
library(ggplotify)
library(patchwork) #to combine ggplots
library(RVAideMemoire)#mood.medhurst
library(FSA) #Dunn's test
library(rstatix) #Kruskall-Wallis
library(fabricatr)
library(patchwork)
library("edgeR")
library(gprofiler2)
library(viridis)
library(ggpointdensity)
library(GGally)
library(boot)


#extract gene/exon info
gtf <- read.table("/Users/.../ncbi_dataset/data/GCF_027887165.1/genomic.gtf", sep="\t", quote = "")
gtf <- dplyr::as_tibble(gtf)

genes <- gtf %>% dplyr::filter(V3 == "gene" | V3 == "exon") %>% dplyr::filter(str_detect(V1, "NC_0772.*")) %>% dplyr::select (-V2, -V6, -V8)

names(genes) <- c("chromosome", "type", "start", "end", "strand", "gene")

#split out INFO, separate ALTs
genes <- concat.split(genes, 6, sep = ";")

#remove text from columns, filter for gene name only - NB this keeps in micro RNAs
genes <- genes %>% mutate(gene_01 = sub("gene_id ","", gene_01)) %>%
  mutate(gene = gene_01) %>%
  dplyr::select(chromosome, type, start, end, strand, gene)

#remove quotations from gene name
genes$gene <- trimws(genes$gene, whitespace = '"')

exons <- genes %>% filter(type == "exon")
exons <- exons %>% mutate(length = end - start)

genes2 <- genes %>% filter(type == "gene") %>% dplyr::select("CHROM" = chromosome, start, end, gene)
exons2 <- exons %>% dplyr::select("CHROM" = chromosome, start, end, gene)

#read in vcf files
setwd("/Users/.../mpileup_vcf/")

mpu <- read.vcfR("MDO_diploid_mpileup_combined_subjunc_mMonDom1.vcf")
codes <- as_tibble(mpu@meta)
as <- mpu@fix

#convert vcf meta-data to tibble #use @ to specify component of S4 object; filter for identified chromosomes only
dip <- as_tibble(mpu@fix) 
dip2 <- as_tibble(mpu@gt)
meta <- as_tibble(mpu@meta)
dip3 <- cbind(dip, dip2)

kd_dip <- dip

kd_dip <- concat.split(kd_dip, 8, sep = ";")

kd_dip$INFO <- str_extract_all(kd_dip$INFO, "DP4=.*")
kd_dip <- concat.split(kd_dip, 8, sep = ";")
colnames(kd_dip)[9] <- "DPinfo"
colnames(kd_dip)[10] <- "ADinfo"
colnames(kd_dip)[28] <- "DP4"
colnames(kd_dip)[29] <- "MQ"
kd_dip <- dplyr::select(kd_dip, CHROM, POS, "ref" = REF, "alt" = ALT,
                             "qual" = QUAL, "MQ", "DPinfo", "ADinfo", "DP4") 

#remove text from columns
kd_dip$DPinfo <- sub("DP=","", kd_dip$DPinfo) 
kd_dip$ADinfo <- sub("AD=","", kd_dip$ADinfo)
kd_dip$DP4 <- sub("DP4=","", kd_dip$DP4)
kd_dip$MQ <- sub("MQ=","", kd_dip$MQ)

#split counts columns
kd_dip <- concat.split(kd_dip, "alt", sep = ",")
kd_dip <- concat.split(kd_dip, "ADinfo", sep = ",")
kd_dip <- concat.split(kd_dip, "DP4", sep = ",")
colnames(kd_dip)[15:18] <- c("DP4_F_ref", "DP4_R_ref", "DP4_F_alt", "DP4_R_alt")

kd_dip <- kd_dip %>% select(CHROM, POS, ref, alt_1, alt_2, qual, MQ, DPinfo, ADinfo_1, ADinfo_2, DP4_F_ref,DP4_R_ref, DP4_F_alt, DP4_R_alt) %>%
  mutate("DPqual" = ADinfo_1 + ADinfo_2, .after = "DPinfo")

#DPqual is read depth of quality reads, DPinfo is raw read depth

kd_dip$POS <- as.numeric(kd_dip$POS)
kd_dip$qual <- as.numeric(kd_dip$qual)
kd_dip$MQ <- as.numeric(kd_dip$MQ)
kd_dip$DPinfo <- as.numeric(kd_dip$DPinfo)

#samp
#_1 is genotype (1/1 = homo ref, 0/0 = homo alt, 0/1 = hetero) 
#_3 (DP) is total reads 
#_2 is genotype probability (phred) 
#_4 is allelic depth (high quality)
#all biallelic - but few ALT + ALT2 (i.e. no REF - so different from reference genome)

samp <- dip3 %>% select(-ID, -REF, -ALT, -QUAL, -FILTER,-INFO)
samp$POS <- as.numeric(samp$POS)

colnames(samp)[4:11] <- c("CAPRIN1", "CKAP4", "HNRNPK", "NCL", "NONO", "NO_VECTOR", "SYNCRIP", "U6M2")
samp <- concat.split(samp, "CAPRIN1", sep = ":", drop = TRUE)
samp <- concat.split(samp, "CKAP4", sep = ":", drop = TRUE)
samp <- concat.split(samp, "HNRNPK", sep = ":", drop = TRUE)
samp <- concat.split(samp, "NCL", sep = ":", drop = TRUE)
samp <- concat.split(samp, "NONO", sep = ":", drop = TRUE)
samp <- concat.split(samp, "NO_VECTOR", sep = ":", drop = TRUE)
samp <- concat.split(samp, "SYNCRIP", sep = ":", drop = TRUE)
samp <- concat.split(samp, "U6M2", sep = ":", drop = TRUE)

samp <- concat.split(samp, "CAPRIN1_4", sep = ",", drop = TRUE)
samp <- concat.split(samp, "CKAP4_4", sep = ",", drop = TRUE)
samp <- concat.split(samp, "HNRNPK_4", sep = ",", drop = TRUE)
samp <- concat.split(samp, "NCL_4", sep = ",", drop = TRUE)
samp <- concat.split(samp, "NONO_4", sep = ",", drop = TRUE)
samp <- concat.split(samp, "NO_VECTOR_4", sep = ",", drop = TRUE)
samp <- concat.split(samp, "SYNCRIP_4", sep = ",", drop = TRUE)
samp <- concat.split(samp, "U6M2_4", sep = ",", drop = TRUE)

samp <- samp %>% select(CHROM, POS,    
                        CAPRIN1_1, CAPRIN1_3, "CAPRIN1_REF" = CAPRIN1_4_1, "CAPRIN1_ALT" = CAPRIN1_4_2, "CAPRIN1_ALT2" = CAPRIN1_4_3,
                        CKAP4_1, CKAP4_3, "CKAP4_REF" = CKAP4_4_1, "CKAP4_ALT" = CKAP4_4_2, "CKAP4_ALT2" = CKAP4_4_3,
                        HNRNPK_1, HNRNPK_3, "HNRNPK_REF" = HNRNPK_4_1, "HNRNPK_ALT" = HNRNPK_4_2, "HNRNPK_ALT2" = HNRNPK_4_3,
                        NCL_1, NCL_3, "NCL_REF" = NCL_4_1, "NCL_ALT" = NCL_4_2, "NCL_ALT2" = NCL_4_3,
                        NONO_1, NONO_3, "NONO_REF" = NONO_4_1, "NONO_ALT" = NONO_4_2, "NONO_ALT2" = NONO_4_3,
                        NO_VECTOR_1, NO_VECTOR_3, "NO_VECTOR_REF" = NO_VECTOR_4_1, "NO_VECTOR_ALT" = NO_VECTOR_4_2, "NO_VECTOR_ALT2" = NO_VECTOR_4_3,
                        SYNCRIP_1, SYNCRIP_3, "SYNCRIP_REF" = SYNCRIP_4_1, "SYNCRIP_ALT" = SYNCRIP_4_2, "SYNCRIP_ALT2" = SYNCRIP_4_3,
                        U6M2_1, U6M2_3, "U6M2_REF" = U6M2_4_1, "U6M2_ALT" = U6M2_4_2, "U6M2_ALT2" = U6M2_4_3)
#replace all NA with 0
samp <- samp %>% replace(is.na(.), 0)


#major allele proportion (map)
samp2 <- samp %>% dplyr::mutate("U6M2_qr" = U6M2_REF + U6M2_ALT + U6M2_ALT2, .after = U6M2_3) %>%
 mutate("U6M2_ALT2p" = U6M2_ALT2/U6M2_qr, .after = U6M2_ALT2) %>% mutate("U6M2_ALTp" = U6M2_ALT/U6M2_qr, .after = U6M2_ALT2) %>% 
  mutate("U6M2_REFp" = U6M2_REF/U6M2_qr, .after = U6M2_ALT2) %>% 
  replace(is.na(.), 0) %>%
  mutate("U6M2_map" = pmax(U6M2_REFp, U6M2_ALTp, U6M2_ALT2p), .after = POS)

samp2 <- samp2 %>% dplyr::mutate("SYNCRIP_qr" = SYNCRIP_REF + SYNCRIP_ALT + SYNCRIP_ALT2, .after = SYNCRIP_3) %>%
  mutate("SYNCRIP_ALT2p" = SYNCRIP_ALT2/SYNCRIP_qr, .after = SYNCRIP_ALT2) %>% mutate("SYNCRIP_ALTp" = SYNCRIP_ALT/SYNCRIP_qr, .after = SYNCRIP_ALT2) %>% 
  mutate("SYNCRIP_REFp" = SYNCRIP_REF/SYNCRIP_qr, .after = SYNCRIP_ALT2) %>% 
  replace(is.na(.), 0) %>%
  mutate("SYNCRIP_map" = pmax(SYNCRIP_REFp, SYNCRIP_ALTp, SYNCRIP_ALT2p), .after = POS)

samp2 <- samp2 %>% dplyr::mutate("NO_VECTOR_qr" = NO_VECTOR_REF + NO_VECTOR_ALT + NO_VECTOR_ALT2, .after = NO_VECTOR_3) %>%
  mutate("NO_VECTOR_ALT2p" = NO_VECTOR_ALT2/NO_VECTOR_qr, .after = NO_VECTOR_ALT2) %>% mutate("NO_VECTOR_ALTp" = NO_VECTOR_ALT/NO_VECTOR_qr, .after = NO_VECTOR_ALT2) %>% 
  mutate("NO_VECTOR_REFp" = NO_VECTOR_REF/NO_VECTOR_qr, .after = NO_VECTOR_ALT2) %>% 
  replace(is.na(.), 0) %>%
  mutate("NO_VECTOR_map" = pmax(NO_VECTOR_REFp, NO_VECTOR_ALTp, NO_VECTOR_ALT2p), .after = POS)

samp2 <- samp2 %>% dplyr::mutate("NONO_qr" = NONO_REF + NONO_ALT + NONO_ALT2, .after = NONO_3) %>%
  mutate("NONO_ALT2p" = NONO_ALT2/NONO_qr, .after = NONO_ALT2) %>% mutate("NONO_ALTp" = NONO_ALT/NONO_qr, .after = NONO_ALT2) %>% 
  mutate("NONO_REFp" = NONO_REF/NONO_qr, .after = NONO_ALT2) %>% 
  replace(is.na(.), 0) %>%
  mutate("NONO_map" = pmax(NONO_REFp, NONO_ALTp, NONO_ALT2p), .after = POS)

samp2 <- samp2 %>% dplyr::mutate("NCL_qr" = NCL_REF + NCL_ALT + NCL_ALT2, .after = NCL_3) %>%
  mutate("NCL_ALT2p" = NCL_ALT2/NCL_qr, .after = NCL_ALT2) %>% mutate("NCL_ALTp" = NCL_ALT/NCL_qr, .after = NCL_ALT2) %>% 
  mutate("NCL_REFp" = NCL_REF/NCL_qr, .after = NCL_ALT2) %>% 
  replace(is.na(.), 0) %>%
  mutate("NCL_map" = pmax(NCL_REFp, NCL_ALTp, NCL_ALT2p), .after = POS)

samp2 <- samp2 %>% dplyr::mutate("HNRNPK_qr" = HNRNPK_REF + HNRNPK_ALT + HNRNPK_ALT2, .after = HNRNPK_3) %>%
  mutate("HNRNPK_ALT2p" = HNRNPK_ALT2/HNRNPK_qr, .after = HNRNPK_ALT2) %>% mutate("HNRNPK_ALTp" = HNRNPK_ALT/HNRNPK_qr, .after = HNRNPK_ALT2) %>% 
  mutate("HNRNPK_REFp" = HNRNPK_REF/HNRNPK_qr, .after = HNRNPK_ALT2) %>% 
  replace(is.na(.), 0) %>%
  mutate("HNRNPK_map" = pmax(HNRNPK_REFp, HNRNPK_ALTp, HNRNPK_ALT2p), .after = POS)

samp2 <- samp2 %>% dplyr::mutate("CKAP4_qr" = CKAP4_REF + CKAP4_ALT + CKAP4_ALT2, .after = CKAP4_3) %>%
  mutate("CKAP4_ALT2p" = CKAP4_ALT2/CKAP4_qr, .after = CKAP4_ALT2) %>% mutate("CKAP4_ALTp" = CKAP4_ALT/CKAP4_qr, .after = CKAP4_ALT2) %>% 
  mutate("CKAP4_REFp" = CKAP4_REF/CKAP4_qr, .after = CKAP4_ALT2) %>% 
  replace(is.na(.), 0) %>%
  mutate("CKAP4_map" = pmax(CKAP4_REFp, CKAP4_ALTp, CKAP4_ALT2p), .after = POS)

samp2 <- samp2 %>% dplyr::mutate("CAPRIN1_qr" = CAPRIN1_REF + CAPRIN1_ALT + CAPRIN1_ALT2, .after = CAPRIN1_3) %>%
  mutate("CAPRIN1_ALT2p" = CAPRIN1_ALT2/CAPRIN1_qr, .after = CAPRIN1_ALT2) %>% mutate("CAPRIN1_ALTp" = CAPRIN1_ALT/CAPRIN1_qr, .after = CAPRIN1_ALT2) %>% 
  mutate("CAPRIN1_REFp" = CAPRIN1_REF/CAPRIN1_qr, .after = CAPRIN1_ALT2) %>% 
  replace(is.na(.), 0) %>%
  mutate("CAPRIN1_map" = pmax(CAPRIN1_REFp, CAPRIN1_ALTp, CAPRIN1_ALT2p), .after = POS)
  
kd_dip2 <- full_join(kd_dip, samp2, by = c("CHROM", "POS"))


#use genome_join from fuzzyjoin package to combine genes dataframe with kd_dip2 + keep other columns
#NB genome_intersect (GenomicRanges) also joins, but doesn't keep other columns
#first re-arrange kd_dip2 to match genes tibble
kd_dip2 <- kd_dip2 %>% mutate(start = as.numeric(POS), end = as.numeric(POS)) %>% relocate(start, .after=CHROM) %>%
  relocate(end, .after=start) 

kd_dip2 <- genome_join(kd_dip2, genes2, by=c("CHROM", "start", "end"), mode = "left") %>% relocate(gene, .after = POS) %>% relocate("gene_start" = start.y, .after = gene) %>%
  relocate("gene_end" = end.y, .before = ref) %>% select(-CHROM.y)

kd_dip2 <- genome_join(kd_dip2, exons2, by=c("CHROM.x" = "CHROM", "start.x" = "start", "end.x" = "end"), mode = "left") %>% relocate("gene" = gene.x, .after = POS) %>% relocate("exon_start" = start, .after = gene_end) %>%
  relocate("exon_end" = end, .before = ref) %>% select(-start.x, -end.x, -gene.y)

#Normalise total mapped read counts using EdgeR. Based on kd_dip2 file (ie before additional filtering)
qcounts <- kd_dip2 %>% select(CHROM.x, POS, CAPRIN1_qr, CKAP4_qr, HNRNPK_qr, NCL_qr, NONO_qr, SYNCRIP_qr, U6M2_qr, NO_VECTOR_qr)
qcounts <- qcounts %>% unite(locus, c("CHROM.x", "POS"))

#convert to DGEList object
#rows = features, columns = samples
qcounts_DGE <- DGEList(counts=qcounts[ ,2:9], genes=qcounts[ ,1])

#normalise TMM (normalises to calculate effective library sizes). Does not change count values 
qcounts_DGE_norm <- calcNormFactors(qcounts_DGE)
head(qcounts_DGE_norm)
#$samples
#group lib.size norm.factors
#CAPRIN1_qr       1 55897568    1.0147760
#CKAP4_qr         1 41406274    1.0279583
#HNRNPK_qr        1 36907676    0.9921610
#NCL_qr           1 37502281    1.0036259
#NONO_qr          1 41957344    0.9992046
#SYNCRIP_qr       1 36650080    0.9794332
#U6M2_qr          1 33367612    0.9665680
#NO_VECTOR_qr     1 37486304    1.0177442

#$genes
#genes
#1 NC_077227.1_499587
#2 NC_077227.1_528114
#3 NC_077227.1_530428
#4 NC_077227.1_540858
#5 NC_077227.1_540877
#6 NC_077227.1_540890

#For sample with lowest coverage (U6M2), 6 quality reads = 0.186 (SNV normalised cpm); 7 quality read = 0.217 --> use 0.2 cpm as minimium
#For sample with highest coverage (CAPRIN1), 11 quality reads = -.19 cpm; 12 q reads = 0.211 --

#normalise using cpm 
keep <- rowSums(cpm(qcounts_DGE_norm)>=0) >= 0
summary(keep)
#Mode    TRUE 
#logical 1193106  

#don't keep library sizes, ie so re-calculated after filtering for lowly expressed genes
qcounts_DGE_filtered <- qcounts_DGE_norm[keep, , keep.lib.sizes=FALSE]
head(qcounts_DGE_filtered)

qcounts_cpm <- cpm(qcounts_DGE_filtered, normalized.lib.sizes = TRUE)
head(qcounts_cpm)
#counts_filtered remains a DGEList object
#counts_normal_cpm is matrix of CPM

colnames(qcounts_cpm) <- c("caprin1_cpm", "ckap4_cpm", "hnrnpk_cpm", "ncl_cpm", "nono_cpm", "syncrip_cpm", "u6m2_cpm", "no_vector_cpm")

kd_dip2 <- cbind(kd_dip2, qcounts_cpm)

#Add in gene counts for genes
#same GTF used for genes2 annotation and for Htseq-count - so just join on gene
diploid_counts_1cpm <- read.table("/Users/.../diploid_kd_all_cpm2_subread_mMonDom1_1cpm.txt", header=TRUE, sep="\t")

counts <- diploid_counts_1cpm %>% select(-chromosome, -start, -end, -strand, -genelength)

kd_dip3 <- left_join(kd_dip2, counts, by = "gene", keep = TRUE) 

#filter out Y chrom = NC_077236.1
kd_dip3 <- kd_dip3 %>% filter(str_detect(CHROM.x, "NC_0772.*")) %>% filter(CHROM.x != "NC_077236.1")

#filter for call quality
kd_dip3 <- kd_dip3 %>% filter(qual >=10)

#filter for depth of quality reads
kd_dip3 <- kd_dip3 %>% filter(DPqual >=20)

#passing quality = 534,553

#filter for biallelic + EXPRESSED > 0.2 CPM - in at least one of kd samples or u6m2 control
kd_dip4 <- kd_dip3 %>% filter((CAPRIN1_map <=0.98 & caprin1_cpm > 0.2) | (CKAP4_map <=0.98 & ckap4_cpm > 0.2) | (HNRNPK_map <=0.98 & hnrnpk_cpm > 0.2) | 
                                (NCL_map <=0.98 & ncl_cpm > 0.2) | (NONO_map <=0.98 & nono_cpm > 0.2) | (SYNCRIP_map <=0.98 & syncrip_cpm > 0.2) | (U6M2_map <=0.98 & u6m2_cpm > 0.2)) 

#biallelic in at least one - with expression filter = 184,994

#check how many biallelic in only one kd sample + control - WITH EXPRESSION MINIMUM
kd_dip4_u6 <- kd_dip3 %>% filter(U6M2_map <=0.98 & u6m2_cpm > 0.2) %>% filter(U6M2_map >=0.5) #131279
kd_dip4_cap <- kd_dip3 %>% filter(CAPRIN1_map <=0.98 & caprin1_cpm > 0.2) %>% filter(CAPRIN1_map >=0.5) #134031
kd_dip4_ckap <- kd_dip3 %>% filter(CKAP4_map <=0.98 & ckap4_cpm > 0.2) %>% filter(CKAP4_map >=0.5) #132957
kd_dip4_hk <- kd_dip3 %>% filter(HNRNPK_map <=0.98 & hnrnpk_cpm > 0.2) %>% filter(HNRNPK_map >=0.5) #131687
kd_dip4_ncl <- kd_dip3 %>% filter(NCL_map <=0.98 & ncl_cpm > 0.2) %>% filter(NCL_map >=0.5) #131687
kd_dip4_nono <- kd_dip3 %>% filter(NONO_map <=0.98 & nono_cpm > 0.2) %>% filter(NONO_map >=0.5) #130552
kd_dip4_syn <- kd_dip3 %>% filter(SYNCRIP_map <=0.98 & syncrip_cpm > 0.2) %>% filter(SYNCRIP_map >=0.5) #130812
kd_dip4_nude <- kd_dip3 %>% filter(NO_VECTOR_map <=0.98 & no_vector_cpm > 0.2) %>% filter(NO_VECTOR_map >=0.5) #131452

setwd("/Users/.../new_combined_vcf_mMonDom1/")

write.table(kd_dip4, file = "kd_dip4.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

kd_dip4 <- read.table("kd_dip4.txt", header=TRUE, sep="\t")

#Figure 4 - FILTERING SUMMARIES  
setwd("/Users/.../new_combined_vcf_mMonDom1/")
kd_dip3a <- kd_dip3 %>% filter(caprin1_cpm > 0.2 | ckap4_cpm > 0.2 | hnrnpk_cpm > 0.2 | ncl_cpm > 0.2 | nono_cpm > 0.2 | syncrip_cpm > 0.2 | u6m2_cpm > 0.2)
#passing quality + minimum exp = 433,144

summary_kd_dip3 <-  kd_dip3a %>% count(CHROM.x)
summary_gene_dip3 <- kd_dip3a %>% filter(gene.x != "NA") %>% count(CHROM.x)
a <- sum(summary_gene_dip3$n) #397903 snvs in genes - after quality + expression filter
summary_gene_dip3exon <- kd_dip3a %>% filter(exon_start != "NA") %>% count(CHROM.x)
b <- sum(summary_gene_dip3exon$n) #302842 snvs in exons - after quality filter
summary_kd_dip4 <- kd_dip4 %>% count(CHROM.x)

summary_gene <- kd_dip4 %>% filter(gene.x != "NA") %>% count(CHROM.x)
c <- sum(summary_gene$n) #170212 snvs in genes - after biallelic +expression filter
summary_exon <- kd_dip4 %>% filter(exon_start != "NA") %>% count(CHROM.x)
d <- sum(summary_exon$n) #129930 snvs in genes - after biallelic filter
summary_sum <- left_join(summary_kd_dip3, summary_kd_dip4, by = "CHROM.x", keep = FALSE)
summary_sum <- left_join(summary_sum, summary_gene, by = "CHROM.x", keep = FALSE)
summary_sum <- left_join(summary_sum, summary_exon, by = "CHROM.x", keep = FALSE)
colnames(summary_sum) <- c("chromosome", "quality", "biallelic", "genic", "exonic") #NB all these subject to minimum SNV cpm > 0.2 in at least one sample

summary_sum <- summary_sum %>% mutate(chromosome = sub("NC_077227.1","1", chromosome)) %>%
  mutate(chromosome = sub("NC_077228.1","2", chromosome)) %>%
  mutate(chromosome = sub("NC_077229.1","3", chromosome)) %>%
  mutate(chromosome = sub("NC_077230.1","4", chromosome)) %>%
  mutate(chromosome = sub("NC_077231.1","5", chromosome)) %>%
  mutate(chromosome = sub("NC_077232.1","6", chromosome)) %>%
  mutate(chromosome = sub("NC_077233.1","7", chromosome)) %>%
  mutate(chromosome = sub("NC_077234.1","8", chromosome)) %>%
  mutate(chromosome = sub("NC_077235.1","X", chromosome)) %>%
  mutate(XA = if_else(chromosome =="X", "X", "A"), .after = chromosome) 

write.table(summary_sum, file = "filtering_summary.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

sum_sum <- summary_sum %>% pivot_longer(cols = c(quality, biallelic, genic, exonic), names_to = "filter")
sum_sum2 <- sum_sum %>% group_by(XA, filter) %>% tally(value)
sum_sum3 <- sum_sum2 %>% group_by(filter) %>% tally(n)

filters <- c("quality", "biallelic", "genic", "exonic")

all_summ_plot_X <- ggplot(sum_sum %>% filter(XA == "X")) +
  geom_col(aes(x = factor(filter, filters), y = value), fill = "turquoise3", alpha = 0.8) +
  ylab("") +
  xlab("") +
  ggtitle("X chromosome")+
  theme_classic(base_size = 18) +
  geom_text(data = . %>% filter(XA == "X"), aes(y = value, x = filter, label = value), vjust = -0.5, size = 5)+
  ylim(0,4500)+
  theme(legend.position = "none")
all_summ_plot_X
ggsave("all_summ_plot_X.svg", width = 7, height = 4, units = "in")

all_summ_plot_A <- ggplot(sum_sum %>% filter(XA == "A")) +
  geom_col(aes(x = factor(filter, filters), y = value), fill = "coral2", alpha = 0.8) +
  ylab("") +
  xlab("") +
  ggtitle("Autosomes")+
  theme_classic(base_size = 18) +
  geom_text(data = sum_sum2 %>% filter(XA == "A"), aes(y = n, x = filter, label = n), vjust = -0.5, size = 5)+
 ylim(0,450000)+
  theme(legend.position = "none")
all_summ_plot_A
ggsave("all_summ_plot_A.svg", width = 7, height = 4, units = "in")


#SNVs with biallelic expression in at least one sample
kd_dip5 <- kd_dip4 %>% dplyr::select("chromosome" = CHROM.x, "pos" = POS, "gene" = gene.x, gene_start, gene_end, exon_start, exon_end,
                              qual, "CAPRIN1" = CAPRIN1_map, "CKAP4" = CKAP4_map, "HNRNPK" = HNRNPK_map, "NCL" = NCL_map, "NONO" = NONO_map,
                              "SYNCRIP" = SYNCRIP_map, "U6M2" = U6M2_map, "NO_VECTOR" = NO_VECTOR_map, CAPRIN1_qr, CKAP4_qr, HNRNPK_qr, NCL_qr, NONO_qr, SYNCRIP_qr, U6M2_qr, NO_VECTOR_qr, "CAPRIN1_genecpm" = CAPRIN1, "CKAP4_genecpm" = CKAP4, "HNRNPK_genecpm" = HNRNPK, "NCL_genecpm" = NCL, "NONO_genecpm" = NONO, "SYNCRIP_genecpm" = SYNCRIP, "U6M2_genecpm" = U6M2, "NO_VECTOR_genecpm" = NUDE,
                              caprin1_cpm, ckap4_cpm, hnrnpk_cpm, ncl_cpm, nono_cpm, syncrip_cpm, u6m2_cpm, no_vector_cpm)
kd_dip5 <- kd_dip5 %>% mutate(chromosome = sub("NC_077227.1","1", chromosome)) %>%
  mutate(chromosome = sub("NC_077228.1","2", chromosome)) %>%
  mutate(chromosome = sub("NC_077229.1","3", chromosome)) %>%
  mutate(chromosome = sub("NC_077230.1","4", chromosome)) %>%
  mutate(chromosome = sub("NC_077231.1","5", chromosome)) %>%
  mutate(chromosome = sub("NC_077232.1","6", chromosome)) %>%
  mutate(chromosome = sub("NC_077233.1","7", chromosome)) %>%
  mutate(chromosome = sub("NC_077234.1","8", chromosome)) %>%
  mutate(chromosome = sub("NC_077235.1","X", chromosome)) %>%
  mutate(XA = if_else(chromosome =="X", "X", "A"), .after = chromosome) 


write.table(kd_dip5, file = "kd_dip5.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
kd_dip5 <- read.table("kd_dip5.txt", header=TRUE, sep="\t")

#check how many biallelic in only one kd sample + control - WITH EXPRESSION MINIMUM
kd_dip5_u6 <- kd_dip5 %>% filter(U6M2 <=0.98 & u6m2_cpm > 0.2) %>% filter(U6M2 >=0.5) #131279
kd_dip5_cap <- kd_dip5 %>% filter(CAPRIN1 <=0.98 & caprin1_cpm > 0.2) %>% filter(CAPRIN1 >=0.5) #134031
kd_dip5_ckap <- kd_dip5 %>% filter(CKAP4 <=0.98 & ckap4_cpm > 0.2) %>% filter(CKAP4 >=0.5) #132957
kd_dip5_hk <- kd_dip5 %>% filter(HNRNPK <=0.98 & hnrnpk_cpm > 0.2) %>% filter(HNRNPK >=0.5) #131687
kd_dip5_ncl <- kd_dip5 %>% filter(NCL <=0.98 & ncl_cpm > 0.2) %>% filter(NCL >=0.5) #131687
kd_dip5_nono <- kd_dip5 %>% filter(NONO <=0.98 & nono_cpm > 0.2) %>% filter(NONO >=0.5) #130552
kd_dip5_syn <- kd_dip5 %>% filter(SYNCRIP <=0.98 & syncrip_cpm > 0.2) %>% filter(SYNCRIP >=0.5) #130812
kd_dip5_nude <- kd_dip5 %>% filter(NO_VECTOR <=0.98 & no_vector_cpm > 0.2) %>% filter(NO_VECTOR >=0.5) #131452

kd_dip5_u6 <- kd_dip5 %>% filter(U6M2 <=0.98 & u6m2_cpm > 0.5) %>% filter(U6M2 >=0.5) #92249
kd_dip5_cap <- kd_dip5 %>% filter(CAPRIN1 <=0.98 & caprin1_cpm > 0.5) %>% filter(CAPRIN1 >=0.5) #93480
kd_dip5_ckap <- kd_dip5 %>% filter(CKAP4 <=0.98 & ckap4_cpm > 0.5) %>% filter(CKAP4 >=0.5) #92769
kd_dip5_hk <- kd_dip5 %>% filter(HNRNPK <=0.98 & hnrnpk_cpm > 0.5) %>% filter(HNRNPK >=0.5) #92915
kd_dip5_ncl <- kd_dip5 %>% filter(NCL <=0.98 & ncl_cpm > 0.5) %>% filter(NCL >=0.5) #92916
kd_dip5_nono <- kd_dip5 %>% filter(NONO <=0.98 & nono_cpm > 0.5) %>% filter(NONO >=0.5) #93165
kd_dip5_syn <- kd_dip5 %>% filter(SYNCRIP <=0.98 & syncrip_cpm > 0.5) %>% filter(SYNCRIP >=0.5) #93404
kd_dip5_nude <- kd_dip5 %>% filter(NO_VECTOR <=0.98 & no_vector_cpm > 0.5) %>% filter(NO_VECTOR >=0.5) #91530


#NUMBER OF SNVS PER GENE - Figure 3B
#USE kd_dip5 (ie, snvs with biallelic expression in at least one sample)
tall <- kd_dip5 %>% filter(gene != "NA" & exon_start != "NA") %>% distinct(chromosome, pos, .keep_all = TRUE) 
tall <- tall %>% count(XA, gene)
tall2 <- tall %>% count(XA, n)
#autosomes
tall3 <- tall2 %>% filter(XA == "A") %>% mutate(group = case_when((n==1) ~ "1",
                                                                  (n > 1 & n < 4) ~ "2-3",
                                                                     (n > 3 & n < 7) ~ "4-6",
                                                                     (n>6 & n<11) ~ "7-10",
                                                                     (n>10 & n<21) ~ "11-20",
                                                                      (n>20 & n<146) ~ "21-145"))
  
tall3$group <- factor(tall3$group, levels = c("1", "2-3", "4-6", "7-10", "11-20", "21-145"))

combined_per_gene_plot <- ggplot(tall3, aes(x = group, y = nn))+
  geom_bar(stat = "identity", fill = "coral2", alpha = 0.8)+
  ylab("")+
  coord_flip()+
  xlab("Exonic SNVs per gene")+
  theme_classic(base_size = 16)
  #theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0))
combined_per_gene_plot
ggsave("snvs_per_gene_A.svg", width = 6, height = 6, units = "in")

#X chromosome
tall3X <- tall2 %>% filter(XA == "X") %>% mutate(group = case_when((n==1) ~ "1",
                                                                   (n > 1 & n < 4) ~ "2-3",
                                                                   (n > 3 & n < 7) ~ "4-6",
                                                                   (n>6) ~ "7-16"))
                                                                   
                                                                   


tall3X$group <- factor(tall3X$group, levels = c("1", "2-3", "4-6", "7-16"))

combined_per_gene_plot <- ggplot(tall3X, aes(x = group, y = nn))+
  geom_bar(stat = "identity", fill = "turquoise3", alpha = 0.8)+
  ylab("")+
  coord_flip()+
  xlab("Exonic SNVs per gene")+
  theme_classic(base_size = 16)
#  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0))
combined_per_gene_plot
ggsave("snvs_per_gene_X.svg", width = 6, height = 6, units = "in")

#FIGURE 5
#Distribution of monoallelic/biallelic/silent SNV expression 
#Of all SNVs expressed bialleically in at least one sample, what proportion are monoallelic, etc in what samples?


#SNVs passing quality filters with a minimum SNV expression level of 0.2 normalised cpm (n = 433,144), 42.7% (n = 184,994)
# ..were expressed biallelically in at least one knockdown sample or the U6M2 control kd_dip6
kd_dip6 <- kd_dip5 %>% pivot_longer(cols = c(CAPRIN1, CKAP4, HNRNPK, NCL, NONO, SYNCRIP, U6M2, NO_VECTOR), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE)

kd_dip6 <- kd_dip6 %>% filter(id != "NO_VECTOR")
kd_dip6$id = factor(kd_dip6$id, levels = c("CAPRIN1", "CKAP4", "SYNCRIP", "HNRNPK", "NCL", "NONO", "U6M2"))

monocount_exon1 <- kd_dip6 %>% filter(exon_start != "NA") %>% filter(caprin1_cpm > 0.2 | ckap4_cpm > 0.2 | hnrnpk_cpm > 0.2 | ncl_cpm > 0.2 | nono_cpm > 0.2 | syncrip_cpm > 0.2 | u6m2_cpm > 0.2)
monocount_exon <- monocount_exon1 %>% filter(map>=0.5 & map <=0.98) %>% dplyr::count(id, XA, .drop = FALSE)
monocount2_exon <- monocount_exon1 %>% filter(map>0.98) %>% count(id, XA, .drop = FALSE)
monocount3_exon <- monocount_exon1 %>% filter(map==0) %>% count(id, XA, .drop = FALSE)

monocount3_exon <- monocount3_exon %>% add_row(id = "CAPRIN1", XA = "X", n = 0) %>% add_row(id = "CKAP4", XA = "X", n = 0) %>% add_row(id = "NCL", XA = "X", n = 0) %>% add_row(id = "NONO", XA = "X", n = 0) %>% add_row(id = "SYNCRIP", XA = "X", n = 0)

monocount4_exon <- left_join(monocount_exon, monocount2_exon, by = c("id", "XA"))
monocount4_exon <- left_join(monocount4_exon, monocount3_exon, by = c("id", "XA"))

colnames(monocount4_exon) <- c("id", "XA", "biallelic", "monoallelic", "silent") 
monocount5_exon <- pivot_longer(monocount4_exon, cols = c(biallelic, monoallelic, silent), names_to = "type")

mono_plot_A_exon <- ggplot(monocount5_exon %>% filter(XA == "A"), aes(x=id, y=value, fill=type))+
  geom_col() +
  coord_flip() +
  ggtitle("Autosomes")+
  ylab("") +
  xlab("") +
  theme_classic() +
  theme(title = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  theme(axis.text.x = element_text(size = 12))+
  theme(legend.text = element_text(size = 14))+
  theme(legend.title = element_text(size = 14))
mono_plot_A_exon
ggsave("mono_plot_A_exon.svg", width = 6.85, height = 4.01, units = "in")

mono_plot_X_exon <- ggplot(monocount5_exon %>% filter(XA == "X"), aes(x=id, y=value, fill=type))+
  geom_col() +
  coord_flip() +
  ggtitle("X chromosome")+
  ylab("") +
  xlab("") +
  theme_classic() +
  theme(title = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  theme(axis.text.x = element_text(size = 12))+
  theme(legend.text = element_text(size = 14))+
  theme(legend.title = element_text(size = 14))
mono_plot_X_exon
ggsave("mono_plot_X_exon.svg", width = 6.85, height = 4.01, units = "in")

#Chi squared test of independence -  X v A
#Pearson's chi-squared test with Yates' continuity correction
#FOR exonIC ONLY

mc <- pivot_wider(monocount5_exon, names_from = c(XA), values_from = c(value))
results_CS <- data.frame(kd = character(0), p.value = numeric(0))
id = unique(mc$id)

for (j in id){
  as <- mc %>% dplyr::filter(id == j) %>% select(-id)
  tst <- chisq.test(as[,2:3])
  p.value <- tst$p.value
  res <- c(j, p.value)
  results_CS[nrow(results_CS)+1,] <- res
}

results_CS <- as_tibble(results_CS)
results_CS$p.value <- as.numeric(results_CS$p.value)
results_CS <- results_CS %>% mutate("p.signif" = case_when(p.value <= 0.0001 ~ "****", p.value <= 0.001 ~ "***", p.value <= 0.01 ~ "**", p.value <= 0.05 ~ "*", p.value > 0.05 ~ "NS"))
results_CS$p.value <- format(results_CS$p.value, scientific = TRUE, digits = 3)
write.table(results_CS, file="Chi_sq_mono_exon.txt", quote=FALSE, sep="\t", row.names = FALSE)
write.table(monocount5_exon, file = "SNV_combined_counts_mono_exon.txt", quote=FALSE, sep="\t", row.names = FALSE)

#Distribution of monoallelic expression - INTRONIC ONLY 
#Of all SNVs expressed biallleically in at least one sample, what proportion are monoallelic in what samples?
kd_dip6 <- kd_dip6 %>% filter(id != "NA")
kd_dip6$id = factor(kd_dip6$id, levels = c("CAPRIN1", "CKAP4", "SYNCRIP", "HNRNPK", "NCL", "NONO", "U6M2"))

monocount_intron1 <- kd_dip6 %>% filter(is.na(exon_start) & gene_start != "NA") %>% filter(caprin1_cpm > 0.2 | ckap4_cpm > 0.2 | hnrnpk_cpm > 0.2 | ncl_cpm > 0.2 | nono_cpm > 0.2 | syncrip_cpm > 0.2 | u6m2_cpm > 0.2)
monocount_intron <- monocount_intron1 %>% filter(map>=0.5 & map <=0.98) %>% dplyr::count(id, XA)
monocount2_intron <- monocount_intron1 %>% filter(map>0.98) %>% count(id, XA)
monocount3_intron <- monocount_intron1 %>% filter(map==0) %>% count(id, XA)
monocount4_intron <- cbind(monocount_intron, monocount2_intron$n, monocount3_intron$n)
colnames(monocount4_intron) <- c("id", "XA", "biallelic", "monoallelic", "silent")
monocount5_intron <- pivot_longer(monocount4_intron, cols = c(biallelic, monoallelic, silent), names_to = "type")

mono_plot_A_intron <- ggplot(monocount5_intron %>% filter(XA == "A"), aes(x=id, y=value, fill=type))+
  geom_col() +
  coord_flip() +
  ggtitle("Autosomes")+
  ylab("") +
  xlab("") +
  theme_classic() +
  theme(title = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  theme(axis.text.x = element_text(size = 12))+
  theme(legend.text = element_text(size = 14))+
  theme(legend.title = element_text(size = 14))
mono_plot_A_intron
ggsave("mono_plot_A_intron.svg", width = 6.85, height = 4.01, units = "in")

mono_plot_X_intron <- ggplot(monocount5_intron %>% filter(XA == "X"), aes(x=id, y=value, fill=type))+
  geom_col() +
  coord_flip() +
  ggtitle("X chromosome")+
  ylab("") +
  xlab("") +
  theme_classic() +
  theme(title = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  theme(axis.text.x = element_text(size = 12))+
  theme(legend.text = element_text(size = 14))+
  theme(legend.title = element_text(size = 14))
mono_plot_X_intron
ggsave("mono_plot_X_intron.svg", width = 6.85, height = 4.01, units = "in")


#Chi squared test of independence -  X v A
#Pearson's chi-squared test with Yates' continuity correction
#FOR INTRONIC ONLY
mc <- pivot_wider(monocount5_intron, names_from = c(XA), values_from = c(value))
results_CS <- data.frame(kd = character(0), p.value = numeric(0))
id = unique(mc$id)

for (j in id){
  as <- mc %>% dplyr::filter(id == j) %>% select(-id)
  tst <- chisq.test(as[,2:3])
  p.value <- tst$p.value
  res <- c(j, p.value)
  results_CS[nrow(results_CS)+1,] <- res
}

results_CS <- as_tibble(results_CS)
results_CS$p.value <- as.numeric(results_CS$p.value)
results_CS <- results_CS %>% mutate("p.signif" = case_when(p.value <= 0.0001 ~ "****", p.value <= 0.001 ~ "***", p.value <= 0.01 ~ "**", p.value <= 0.05 ~ "*", p.value > 0.05 ~ "NS"))
results_CS$p.value <- format(results_CS$p.value, scientific = TRUE, digits = 3)
write.table(results_CS, file="Chi_sq_mono_intron.txt", quote=FALSE, sep="\t", row.names = FALSE)
write.table(monocount5_intron, file = "SNV_combined_counts_mono_intron.txt", quote=FALSE, sep="\t", row.names = FALSE)


##FIGURE 6 Pairwise correlation of SNV MAPs - by SNV location in gene
#GGPAIR Correlations
#run alternatiely for autosomal and X
setwd("/Users/kimmcintyre/Library/CloudStorage/OneDrive-UNSW/SNP_calling_bcf/new_combined_vcf_mMonDom1/")

kd_dip5ex <- kd_dip5 %>% filter(exon_start != "NA") %>% filter(XA == "A") %>% distinct(pos, gene, .keep_all = TRUE)#Autosomal
#kd_dip5ex <- kd_dip5 %>% filter(exon_start != "NA") %>% filter(XA == "A") %>% distinct(pos, gene, .keep_all = TRUE)#X

as <- kd_dip5ex %>% dplyr::count(gene) %>% filter(n > 1) %>% dplyr::select(-n)

kd_dip5ex3 <- left_join(as, kd_dip5ex) %>% group_by(gene) %>% arrange(pos) %>% mutate(numbering = row_number()) %>% 
  dplyr::select(-chromosome, -XA, -gene_start, -gene_end, -exon_start, -exon_end, -qual)

kd3_cap <- kd_dip5ex3 %>% filter(caprin1_cpm>0.2) %>% dplyr::select(gene, pos, CAPRIN1, CAPRIN1_genecpm, numbering) %>% dplyr::rename(map = 3, genecpm = 4) %>% filter(genecpm > 1) %>% filter(map >= 0.5)
kd3_cap2 <- kd3_cap %>% group_by(gene) %>% pivot_wider(names_from = numbering, values_from = c(map, pos)) %>% ungroup() %>% relocate(gene, .after = pos_16) %>% relocate(genecpm, .after = gene)
colnames(kd3_cap2) <- as.character(1:ncol(kd3_cap2))

kd3_ckap <- kd_dip5ex3 %>% filter(ckap4_cpm>0.2) %>% dplyr::select(gene, pos, CKAP4, CKAP4_genecpm, numbering) %>% dplyr::rename(map = 3, genecpm = 4) %>% filter(genecpm > 1) %>% filter(map >= 0.5)
kd3_ckap2 <- kd3_ckap %>% group_by(gene) %>% pivot_wider(names_from = numbering, values_from = c(map, pos)) %>% ungroup() %>% relocate(gene, .after = pos_16) %>% relocate(genecpm, .after = gene)
colnames(kd3_ckap2) <- as.character(1:ncol(kd3_ckap2))

kd3_hk <- kd_dip5ex3 %>% filter(hnrnpk_cpm>0.2) %>% dplyr::select(gene, pos, HNRNPK, HNRNPK_genecpm, numbering) %>% dplyr::rename(map = 3, genecpm = 4) %>% filter(genecpm > 1) %>% filter(map >= 0.5)
kd3_hk2 <- kd3_hk %>% group_by(gene) %>% pivot_wider(names_from = numbering, values_from = c(map, pos)) %>% ungroup() %>% relocate(gene, .after = pos_16) %>% relocate(genecpm, .after = gene)
colnames(kd3_hk2) <- as.character(1:ncol(kd3_hk2))

kd3_ncl <- kd_dip5ex3 %>% filter(ncl_cpm>0.2) %>% dplyr::select(gene, pos, NCL, NCL_genecpm, numbering) %>% dplyr::rename(map = 3, genecpm = 4) %>% filter(genecpm > 1) %>% filter(map >= 0.5)
kd3_ncl2 <- kd3_ncl %>% group_by(gene) %>% pivot_wider(names_from = numbering, values_from = c(map, pos)) %>% ungroup() %>% relocate(gene, .after = pos_16) %>% relocate(genecpm, .after = gene)
colnames(kd3_ncl2) <- as.character(1:ncol(kd3_ncl2))

kd3_nono <- kd_dip5ex3 %>% filter(nono_cpm>0.2) %>% dplyr::select(gene, pos, NONO, NONO_genecpm, numbering) %>% dplyr::rename(map = 3, genecpm = 4) %>% filter(genecpm > 1) %>% filter(map >= 0.5)
kd3_nono2 <- kd3_nono %>% group_by(gene) %>% pivot_wider(names_from = numbering, values_from = c(map, pos)) %>% ungroup() %>% relocate(gene, .after = pos_16) %>% relocate(genecpm, .after = gene)
colnames(kd3_nono2) <- as.character(1:ncol(kd3_nono2))

kd3_syncrip <- kd_dip5ex3 %>% filter(syncrip_cpm>0.2) %>% dplyr::select(gene, pos, SYNCRIP, SYNCRIP_genecpm, numbering) %>% dplyr::rename(map = 3, genecpm = 4) %>% filter(genecpm > 1) %>% filter(map >= 0.5)
kd3_syncrip2 <- kd3_syncrip %>% group_by(gene) %>% pivot_wider(names_from = numbering, values_from = c(map, pos)) %>% ungroup() %>% relocate(gene, .after = pos_16) %>% relocate(genecpm, .after = gene)
colnames(kd3_syncrip2) <- as.character(1:ncol(kd3_syncrip2))

kd3_u6 <- kd_dip5ex3 %>% filter(u6m2_cpm>0.2) %>% dplyr::select(gene, pos, U6M2, U6M2_genecpm, numbering) %>% dplyr::rename(map = 3, genecpm = 4) %>% filter(genecpm > 1) %>% filter(map >= 0.5)
kd3_u62 <- kd3_u6 %>% group_by(gene) %>% pivot_wider(names_from = numbering, values_from = c(map, pos)) %>% ungroup() %>% relocate(gene, .after = pos_16) %>% relocate(genecpm, .after = gene)
colnames(kd3_u62) <- as.character(1:ncol(kd3_u62))

kd3_nv <- kd_dip5ex3 %>% dplyr::select(gene, pos, NO_VECTOR, NO_VECTOR_genecpm, numbering) %>% dplyr::rename(map = 3, genecpm = 4) %>% filter(genecpm > 1) %>% filter(map >= 0.5)
kd3_nv2 <- kd3_nv %>% group_by(gene) %>% pivot_wider(names_from = numbering, values_from = c(map, pos)) %>% ungroup() %>% relocate(gene, .after = pos_16) %>% relocate(genecpm, .after = gene)
colnames(kd3_nv2) <- as.character(1:ncol(kd3_nv2))

#HEATPLOT MATRIX - Autosomal
setwd("/Users/.../snv_correlations/snv_correlations_autosomal/")

cap2_matrix <- ggcorr(kd3_cap2[,1:20], midpoint = 0.4, limits = c(0.05,0.8)) + theme_classic(base_size = 24) + labs(subtitle = "CAPRIN1")
cap2_matrix
ggsave("snvheat_cap_A_20_snv.svg", width = 10, height = 10, units = "in")

ckap2_matrix <- ggcorr(kd3_ckap2[,1:20], midpoint = 0.4, limits = c(0.05,0.8)) + theme_classic(base_size = 24) + labs(subtitle = "CKAP4")
ckap2_matrix
ggsave("snvheat_ckap_A_20_snv.svg", width = 10, height = 10, units = "in")

hk2_matrix <- ggcorr(kd3_hk2[,1:20], midpoint = 0.4, limits = c(0.05,0.8)) + theme_classic(base_size = 24) + labs(subtitle = "HNRNPK")
hk2_matrix
ggsave("snvheat_hk_A_20_snv.svg", width = 10, height = 10, units = "in")

ncl2_matrix <- ggcorr(kd3_ncl2[,1:20], midpoint = 0.4, limits = c(0.05,0.8)) + theme_classic(base_size = 24) + labs(subtitle = "NCL")
ncl2_matrix
ggsave("snvheat_ncl_A_20_snv.svg", width = 10, height = 10, units = "in")

nono2_matrix <- ggcorr(kd3_nono2[,1:20], midpoint = 0.4, limits = c(0.05,0.8)) + theme_classic(base_size = 24) + labs(subtitle = "NONO")
nono2_matrix
ggsave("snvheat_nono_A_20_snv.svg", width = 10, height = 10, units = "in")

syncrip2_matrix <- ggcorr(kd3_syncrip2[,1:20], midpoint = 0.4, limits = c(0.05,0.8)) + theme_classic(base_size = 24) + labs(subtitle = "SYNCRIP")
syncrip2_matrix
ggsave("snvheat_syncrip_A_20_snv.svg", width = 10, height = 10, units = "in")

u62_matrix <- ggcorr(kd3_u62[,1:20], midpoint = 0.4, limits = c(0.05,0.8)) + theme_classic(base_size = 24) + labs(subtitle = "U6M2")
u62_matrix
ggsave("snvheat_u6_A_20_snv.svg", width = 10, height = 10, units = "in")

nv2_matrix <- ggcorr(kd3_nv2[,1:30], midpoint = 0.4, limits = c(0,0.8)) + theme_classic(base_size = 24) + labs(subtitle = "NO VECTOR")
nv2_matrix
ggsave("snvheat_nv_A_30_snv.svg", width = 10, height = 10, units = "in")


p1 <- as.ggplot(cap2_matrix)
p2 <- as.ggplot(ckap2_matrix)
p3 <- as.ggplot(hk2_matrix)
p4 <- as.ggplot(ncl2_matrix)
p5 <- as.ggplot(nono2_matrix)
p6 <- as.ggplot(syncrip2_matrix)
p7 <- as.ggplot(u62_matrix)

# use patchwork to arrange them together
plots <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + plot_layout(ncol = 4) 
plots

ggsave(file = "snv_heatmaps_all_exon_02_autosomal.svg", plot = plots,  device = "svg", width = 15, height = 13, units = "in")

#HEATPLOT MATRIX - X chromosome
setwd("/Users/.../snv_correlations/snv_correlations_X/")

cap2_matrix <- ggcorr(kd3_cap2[,1:3], midpoint = 0.4, limits = c(0.05,0.8)) + theme_classic(base_size = 24) + labs(subtitle = "CAPRIN1")
cap2_matrix
ggsave("snvheat_cap_X_20_snv.svg", width = 10, height = 10, units = "in")

ckap2_matrix <- ggcorr(kd3_ckap2[,1:3], midpoint = 0.4, limits = c(0.05,0.8)) + theme_classic(base_size = 24) + labs(subtitle = "CKAP4")
ckap2_matrix
ggsave("snvheat_ckap_X_20_snv.svg", width = 10, height = 10, units = "in")

hk2_matrix <- ggcorr(kd3_hk2[,1:3], midpoint = 0.4, limits = c(0.05,0.8)) + theme_classic(base_size = 24) + labs(subtitle = "HNRNPK")
hk2_matrix
ggsave("snvheat_hk_X_20_snv.svg", width = 10, height = 10, units = "in")

ncl2_matrix <- ggcorr(kd3_ncl2[,1:3], midpoint = 0.4, limits = c(0.05,0.8)) + theme_classic(base_size = 24) + labs(subtitle = "NCL")
ncl2_matrix
ggsave("snvheat_ncl_X_20_snv.svg", width = 10, height = 10, units = "in")

nono2_matrix <- ggcorr(kd3_nono2[,1:3], midpoint = 0.4, limits = c(0.05,0.8)) + theme_classic(base_size = 24) + labs(subtitle = "NONO")
nono2_matrix
ggsave("snvheat_nono_X_20_snv.svg", width = 10, height = 10, units = "in")

syncrip2_matrix <- ggcorr(kd3_syncrip2[,1:3], midpoint = 0.4, limits = c(0.05,0.8)) + theme_classic(base_size = 24) + labs(subtitle = "SYNCRIP")
syncrip2_matrix
ggsave("snvheat_syncrip_X_20_snv.svg", width = 10, height = 10, units = "in")

u62_matrix <- ggcorr(kd3_u62[,1:3], midpoint = 0.4, limits = c(0.05,0.8)) + theme_classic(base_size = 24) + labs(subtitle = "U6M2")
u62_matrix
ggsave("snvheat_u6_X_20_snv.svg", width = 10, height = 10, units = "in")

p1 <- as.ggplot(cap2_matrix)
p2 <- as.ggplot(ckap2_matrix)
p3 <- as.ggplot(hk2_matrix)
p4 <- as.ggplot(ncl2_matrix)
p5 <- as.ggplot(nono2_matrix)
p6 <- as.ggplot(syncrip2_matrix)
p7 <- as.ggplot(u62_matrix)

# use patchwork to arrange them together
plots <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + plot_layout(ncol = 4) 
plots
ggsave(file = "snv_heatmaps_all_exon_02_X.svg", plot = plots,  device = "svg", width = 15, height = 13, units = "in")


#FIGURE 7 - SCATTERPLOTS
#Correlation with gene expression based on RNA-seq cpm counts 

setwd("/Users/.../new_combined_vcf_mMonDom1/")

#FIGURE 7A
kd_dip9_caprin <- kd_dip5 %>% filter(U6M2 >=0.5 | CAPRIN1 >=0.5 ) %>% filter(u6m2_cpm >0.5 | caprin1_cpm>0.5) %>% dplyr::select(-CKAP4, -HNRNPK, -NCL, -NONO, -SYNCRIP, -NO_VECTOR) %>% pivot_longer(cols = c(CAPRIN1, U6M2), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE) %>% dplyr::rename(kd = CAPRIN1_genecpm)  
kd_dip9_ckap <- kd_dip5 %>% filter(U6M2 >=0.5 | CKAP4 >=0.5 ) %>% filter(u6m2_cpm >0.5 | ckap4_cpm >0.5) %>% dplyr::select(-CAPRIN1, -HNRNPK, -NCL, -NONO, -SYNCRIP, -NO_VECTOR) %>% pivot_longer(cols = c(CKAP4, U6M2), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE) %>% dplyr::rename(kd = CKAP4_genecpm)
kd_dip9_hnrnpk <- kd_dip5 %>% filter(U6M2 >=0.5  | HNRNPK >=0.5 ) %>% filter(u6m2_cpm >0.5 | hnrnpk_cpm >0.5) %>% dplyr::select(-CKAP4, -CAPRIN1, -NCL, -NONO, -SYNCRIP, -NO_VECTOR) %>% pivot_longer(cols = c(HNRNPK, U6M2), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE) %>% dplyr::rename(kd = HNRNPK_genecpm)
kd_dip9_ncl <- kd_dip5 %>% filter(U6M2 >=0.5  | NCL >=0.5 ) %>% filter(u6m2_cpm >0.5 | ncl_cpm >0.5) %>% dplyr::select(-CKAP4, -HNRNPK, -CAPRIN1, -NONO, -SYNCRIP, -NO_VECTOR) %>% pivot_longer(cols = c(NCL, U6M2), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE) %>% dplyr::rename(kd = NCL_genecpm)
kd_dip9_nono <- kd_dip5 %>% filter(U6M2 >=0.5 | NONO >=0.5 ) %>% filter(u6m2_cpm >0.5 | nono_cpm >0.5) %>% dplyr::select(-CKAP4, -HNRNPK, -NCL, -CAPRIN1, -SYNCRIP, -NO_VECTOR) %>% pivot_longer(cols = c(NONO, U6M2), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE) %>% dplyr::rename(kd = NONO_genecpm)
kd_dip9_syncrip <- kd_dip5 %>% filter(U6M2 >=0.5  | SYNCRIP >=0.5 ) %>% filter(u6m2_cpm >0.5 | syncrip_cpm >0.5) %>% dplyr::select(-CKAP4, -HNRNPK, -NCL, -NONO, -CAPRIN1, -NO_VECTOR) %>% pivot_longer(cols = c(SYNCRIP, U6M2), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE) %>% dplyr::rename(kd = SYNCRIP_genecpm)

dflist <- list(kd_dip9_caprin, kd_dip9_ckap, kd_dip9_hnrnpk, kd_dip9_ncl, kd_dip9_nono, kd_dip9_syncrip)
num_files = length(dflist)
id <- list("CAPRIN1", "CKAP4", "HNRNPK", "NCL", "NONO", "SYNCRIP")
cpm_scatters = list()

for (j in 1:6){
  
  kd_dip9 <- dflist[[j]]
  kd <- id[[j]]
  cpm_scatters[[j]] <- ggplot(kd_dip9, aes(x = log2(U6M2_genecpm), y = log2(kd)))+
    geom_point(aes(color = XA), size = 0.3)+
    geom_abline(slope = 1, intercept = 0, color = "black") +
    theme_classic(base_size = 16) +
    ylab("")+
    xlab("")+
    theme(legend.position = "none")+
    labs(subtitle = paste0(kd))+
    ggpubr::stat_cor(data = kd_dip9 %>% filter(XA == "A"), aes(label = after_stat(paste0("A_",..r..))), label.x = 2, label.y = 13, method = "pearson", size = 5)+
    ggpubr::stat_cor(data = kd_dip9 %>% filter(XA == "X"), aes(label = after_stat(paste0("X_",..r..))), label.x = 2, label.y = 11, method = "pearson", size = 5)
  print(cpm_scatters[[j]])
}
p1 <- as.ggplot(cpm_scatters[[1]])
p2 <- as.ggplot(cpm_scatters[[2]])
p3 <- as.ggplot(cpm_scatters[[3]])
p4 <- as.ggplot(cpm_scatters[[4]])
p5 <- as.ggplot(cpm_scatters[[5]])
p6 <- as.ggplot(cpm_scatters[[6]])

# use patchwork to arrange them together
plots <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 1) 
plots
ggsave(file = "gene_cpm_log2_scatterplots.tiff", plot = plots,  device = "tiff", width = 4, height = 16, units = "in")


#FIGURE 7B,C - Correlation of SNV expression (individual SNVs) v gene in which SNV located (from RNA-seq)
#Run alternately for exonic and intronic SNVs
kd_dip9_caprin <- kd_dip5 %>% filter( CAPRIN1 >=0.5) %>% filter(caprin1_cpm>0.2) %>% filter(CAPRIN1_genecpm>1) %>% select(-CKAP4, -HNRNPK, -NCL, -NONO, -SYNCRIP) %>% pivot_longer(cols = c(CAPRIN1, U6M2), names_to = "id", values_to = "map") %>% dplyr::rename(kd_qr = caprin1_cpm) %>% dplyr::rename(rseq = CAPRIN1_genecpm) %>% distinct(chromosome, pos, id, .keep_all = TRUE) 
kd_dip9_ckap <- kd_dip5 %>% filter(CKAP4 >=0.5 ) %>% filter( ckap4_cpm >0.2) %>% filter(CKAP4_genecpm>1) %>% select(-CAPRIN1, -HNRNPK, -NCL, -NONO, -SYNCRIP) %>% pivot_longer(cols = c(CKAP4, U6M2), names_to = "id", values_to = "map") %>% dplyr::rename(kd_qr = ckap4_cpm) %>% dplyr::rename(rseq = CKAP4_genecpm) %>% distinct(chromosome, pos, id, .keep_all = TRUE)
kd_dip9_hnrnpk <- kd_dip5 %>% filter( HNRNPK >=0.5 ) %>% filter( hnrnpk_cpm >0.2) %>% filter(HNRNPK_genecpm>1) %>% select(-CKAP4, -CAPRIN1, -NCL, -NONO, -SYNCRIP) %>% pivot_longer(cols = c(HNRNPK, U6M2), names_to = "id", values_to = "map") %>% dplyr::rename(kd_qr = hnrnpk_cpm)  %>% dplyr::rename(rseq = HNRNPK_genecpm)  %>% distinct(chromosome, pos, id, .keep_all = TRUE)
kd_dip9_ncl <- kd_dip5 %>% filter( NCL >=0.5 ) %>% filter(ncl_cpm >0.2) %>% filter(NCL_genecpm>1) %>% select(-CKAP4, -HNRNPK, -CAPRIN1, -NONO, -SYNCRIP) %>% pivot_longer(cols = c(NCL, U6M2), names_to = "id", values_to = "map") %>% dplyr::rename(kd_qr = ncl_cpm) %>% dplyr::rename(rseq = NCL_genecpm)  %>% distinct(chromosome, pos, id, .keep_all = TRUE)
kd_dip9_nono <- kd_dip5 %>% filter( NONO >=0.5 ) %>% filter( nono_cpm >0.2) %>% filter(NONO_genecpm>1) %>% select(-CKAP4, -HNRNPK, -NCL, -CAPRIN1, -SYNCRIP) %>% pivot_longer(cols = c(NONO, U6M2), names_to = "id", values_to = "map") %>% dplyr::rename(kd_qr = nono_cpm) %>% dplyr::rename(rseq = NONO_genecpm)  %>% distinct(chromosome, pos, id, .keep_all = TRUE)
kd_dip9_syncrip <- kd_dip5 %>% filter( SYNCRIP >=0.5 ) %>% filter( syncrip_cpm >0.2) %>% filter(SYNCRIP_genecpm>1) %>% select(-CKAP4, -HNRNPK, -NCL, -NONO, -CAPRIN1) %>% pivot_longer(cols = c(SYNCRIP, U6M2), names_to = "id", values_to = "map") %>% dplyr::rename(kd_qr = syncrip_cpm) %>% dplyr::rename(rseq = SYNCRIP_genecpm)  %>% distinct(chromosome, pos, id, .keep_all = TRUE)
kd_dip9_u6 <- kd_dip5 %>% filter( U6M2 >=0.5 ) %>% filter( u6m2_cpm >0.2) %>% filter(U6M2_genecpm>1) %>% select(-CKAP4, -HNRNPK, -NCL, -NONO, -CAPRIN1) %>% pivot_longer(cols = c(SYNCRIP, U6M2), names_to = "id", values_to = "map") %>% dplyr::rename(kd_qr = u6m2_cpm) %>% dplyr::rename(rseq = U6M2_genecpm)  %>% distinct(chromosome, pos, id, .keep_all = TRUE)


dflist <- list(kd_dip9_caprin, kd_dip9_ckap, kd_dip9_hnrnpk, kd_dip9_ncl, kd_dip9_nono, kd_dip9_syncrip, kd_dip9_u6)
num_files = length(dflist)
id <- list("CAPRIN1", "CKAP4", "HNRNPK", "NCL", "NONO", "SYNCRIP", "U6M2")
snv_count_vcpm_scatters = list()

for (j in 1:7){
  kd_dip9 <- dflist[[j]]
  kd <- id[[j]]
  
  kd_dip9 <- kd_dip9 %>% filter(exon_start != "NA") #for exonic SNVs
  #kd_dip9 <- kd_dip9 %>% dplyr::filter(gene != "NA" & is.na(exon_start)) #for intronic SNVs
  
  kd_dip9 <- kd_dip9 %>% distinct(chromosome, gene, .keep_all = TRUE)
  
  snv_count_vcpm_scatters[[j]] <- ggplot(kd_dip9, aes(x = log2(rseq), y = log2(kd_qr)))+
    geom_point(aes(color = XA), size = 0.3)+
    theme_classic(base_size = 16) +
    ylab("")+
    xlab("")+
    theme(legend.position = "none")+
    labs(subtitle = paste0(kd))+
    ggpubr::stat_cor(data = kd_dip9 %>% filter(XA == "A"), aes(label = after_stat(paste0("A_",..r..))), label.x = 0, label.y = 6.4, method = "pearson", size = 5)+
    ggpubr::stat_cor(data = kd_dip9 %>% filter(XA == "X"), aes(label = after_stat(paste0("X_",..r..))), label.x = 0, label.y = 5, method = "pearson", size = 5)
  print(snv_count_vcpm_scatters[[j]])
}

p1 <- as.ggplot(snv_count_vcpm_scatters[[1]])
p2 <- as.ggplot(snv_count_vcpm_scatters[[2]])
p3 <- as.ggplot(snv_count_vcpm_scatters[[3]])
p4 <- as.ggplot(snv_count_vcpm_scatters[[4]])
p5 <- as.ggplot(snv_count_vcpm_scatters[[5]])
p6 <- as.ggplot(snv_count_vcpm_scatters[[6]])
p7 <- as.ggplot(snv_count_vcpm_scatters[[7]])
# use patchwork to arrange them together
plots <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + plot_layout(ncol = 1) 
plots
ggsave(file = "snv_count_scatterplots_total_bygene_indSNVs_exons.tiff", plot = plots,  device = "tiff", width = 4, height = 18, units = "in")
#ggsave(file = "snv_count_scatterplots_total_bygene_indSNVs_introns.tiff", plot = plots,  device = "tiff", width = 4, height = 18, units = "in")

#FIGURE 8 X chromosome SNV expression -by gene + associated tables
setwd("/Users/.../new_combined_vcf_mMonDom1/")

#Filter for exonic X SMVs
xgenes <- kd_dip5 %>% filter(XA == "X") %>% filter(exon_start != "NA") %>% distinct(gene, pos, .keep_all = TRUE) %>% group_by(gene) %>% mutate(count = n()) %>% ungroup()

xgenes3 <- xgenes %>% group_by(gene) %>% mutate(capmed = median(CAPRIN1), ckapmed = median(CKAP4), hkmed = median(HNRNPK), nclmed = median(NCL),
                                                nonomed = median(NONO), synmed = median(SYNCRIP), u6med = median(U6M2), novectormed = median(NO_VECTOR)) %>% ungroup %>% distinct(chromosome, gene, .keep_all = TRUE) 

xgenes3 <- xgenes3 %>% filter(CAPRIN1_genecpm > 1 | CKAP4_genecpm > 1 | HNRNPK_genecpm > 1 | NCL_genecpm > 1 | NONO_genecpm > 1 |
                                SYNCRIP_genecpm > 1 | NO_VECTOR_genecpm > 1)

xgenes4 <- xgenes3 %>% dplyr::select(pos, gene, gene_start, gene_end, capmed, ckapmed, hkmed, nclmed, nonomed, synmed, u6med, novectormed,
                                     CAPRIN1_genecpm, CKAP4_genecpm, HNRNPK_genecpm, NCL_genecpm, NONO_genecpm, SYNCRIP_genecpm, U6M2_genecpm, NO_VECTOR_genecpm, count)

write.table(xgenes4, file = "Xgenes_monoallelic.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

options(scipen = 999) #to suppress scientific notation

#PLOT in 3 panels based on X chromosome locations
cf1 <- ggplot(xgenes4 %>% filter(gene_start < 10000000))+
  geom_segment(aes(x = gene_start, xend = gene_start, y = 0.5, yend = 1), alpha = 0.3)+
  geom_point(aes(x = gene_start, y = capmed), colour = "red", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = ckapmed), colour = "orange", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = hkmed), colour = "green", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = nclmed), colour = "pink1", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = nonomed), colour = "purple", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = synmed), colour = "yellow3", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = u6med), colour = "blue", shape = 18, size = 3, alpha = 0.7)+
  geom_text(aes(x = gene_start, y = 1.05, label = gene), angle = 90, hjust = 0)+
  geom_text(aes(x = gene_start, y = 0.47, label = count), hjust = 0.5)+
  theme_classic(base_size = 18)+
  coord_cartesian(ylim = c(0.45, 1.5))+
  ylab("")+
  xlab("")
cf1

cf2 <- ggplot(xgenes4 %>% filter(gene_start > 10000000 & gene_start < 40000000))+
  geom_segment(aes(x = gene_start, xend = gene_start, y = 0.5, yend = 1), alpha = 0.3)+
  geom_point(aes(x = gene_start, y = capmed), colour = "red", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = ckapmed), colour = "orange", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = hkmed), colour = "green", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = nclmed), colour = "pink1", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = nonomed), colour = "purple", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = synmed), colour = "yellow3", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = u6med), colour = "blue", shape = 18, size = 3, alpha = 0.7)+
  geom_text(aes(x = gene_start, y = 1.05, label = gene), angle = 90, hjust = 0)+
  geom_text(aes(x = gene_start, y = 0.47, label = count), hjust = 0.5)+
  theme_classic(base_size = 16)+
  coord_cartesian(ylim = c(0.45, 1.5))+
  ylab("")+
  xlab("")
cf2

cf3 <- ggplot(xgenes4 %>% filter(gene_start > 40000000))+
  geom_segment(aes(x = gene_start, xend = gene_start, y = 0.5, yend = 1), alpha = 0.3)+
  geom_point(aes(x = gene_start, y = capmed), colour = "red", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = ckapmed), colour = "orange", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = hkmed), colour = "green", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = nclmed), colour = "pink1", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = nonomed), colour = "purple", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = synmed), colour = "yellow3", shape = 18, size = 3, alpha = 0.7)+
  geom_point(aes(x = gene_start, y = u6med), colour = "blue", shape = 18, size = 3, alpha = 0.7)+
  geom_text(aes(x = gene_start, y = 1.05, label = gene), angle = 90, hjust = 0)+
  geom_text(aes(x = gene_start, y = 0.47, label = count), hjust = 0.5)+
  theme_classic(base_size = 18)+
  coord_cartesian(ylim = c(0.45, 1.5))+
  ylab("")+
  xlab("")
cf3

p1 <- as.ggplot(cf1)
p2 <- as.ggplot(cf2)
p3 <- as.ggplot(cf3)

# use patchwork to arrange them together
plots <- p1 + p2 + p3+ plot_layout(ncol = 1) 
plots

ggsave(file = "Xgenes_allelic.svg", plot = plots, device = "svg", width = 16, height = 13, units = "in")

#Genes monoallelic in all samples
sd <- xgenes4 %>% filter(capmed > 0.95 & ckapmed > 0.95 & hkmed > 0.95 & nclmed > 0.95 &
                           nonomed > 0.95 & synmed > 0.95 & u6med > 0.95)


#TABLE 3-7
#Genes monoallelic in control, but escape XCI in knockdown
sd2 <- xgenes4 %>% filter(u6med > 0.98) %>% filter(capmed < 0.9 | ckapmed < 0.9 | hkmed < 0.9 | nclmed < 0.9 |
                                                     nonomed < 0.9 | synmed < 0.9)

#TABLE 3-6
#Genes monoallelic in knockdown, but escape XCI in U6M2
sd3 <- xgenes4 %>% filter(u6med < 0.9) %>% filter(capmed > 0.98 | ckapmed > 0.98 | hkmed > 0.98 | nclmed > 0.98 |
                                                    nonomed > 0.98 | synmed > 0.98)

#TABLE 3-5 Genes escaping XCI
xgenes_cap1 <- xgenes4 %>% dplyr::select(gene, capmed, CAPRIN1_genecpm) %>% dplyr::filter(CAPRIN1_genecpm >= 1)
xgenes_cap2 <- xgenes_cap1 %>% dplyr::filter(capmed<=0.9)
med_map_cap <- median(xgenes_cap2$capmed)#0.7017

xgenes_ck1 <- xgenes4 %>% dplyr::select(gene, ckapmed, CKAP4_genecpm) %>% dplyr::filter(CKAP4_genecpm >= 1)
xgenes_ck2 <- xgenes_ck1 %>% dplyr::filter(ckapmed<=0.9)
med_map_ck <- median(xgenes_ck2$ckapmed)#0.7083

xgenes_hk1 <- xgenes4 %>% dplyr::select(gene, hkmed, HNRNPK_genecpm) %>% dplyr::filter(HNRNPK_genecpm >= 1)
xgenes_hk2 <- xgenes_hk1 %>% dplyr::filter(hkmed<=0.9)
med_map_hk <- median(xgenes_hk2$hkmed)#0.65625

xgenes_ncl1 <- xgenes4 %>% dplyr::select(gene, nclmed, NCL_genecpm) %>% dplyr::filter(NCL_genecpm >= 1)
xgenes_ncl2 <- xgenes_ncl1 %>% dplyr::filter(nclmed<=0.9)
med_map_ncl <- median(xgenes_ncl2$nclmed)#0.68333

xgenes_nono1 <- xgenes4 %>% dplyr::select(gene, nonomed, NONO_genecpm) %>% dplyr::filter(NONO_genecpm >= 1)
xgenes_nono2 <- xgenes_nono1 %>% dplyr::filter(nonomed<=0.9)
med_map_nono <- median(xgenes_nono2$nonomed)#0.6857

xgenes_syn1 <- xgenes4 %>% dplyr::select(gene, synmed, SYNCRIP_genecpm) %>% dplyr::filter(SYNCRIP_genecpm >= 1)
xgenes_syn2 <- xgenes_syn1 %>% dplyr::filter(synmed<=0.9)
med_map_syn <- median(xgenes_syn2$synmed)#0.705

xgenes_u61 <- xgenes4 %>% dplyr::select(gene, u6med, U6M2_genecpm) %>% dplyr::filter(U6M2_genecpm >= 1)
xgenes_u62 <- xgenes_u61 %>% dplyr::filter(u6med<=0.9)
med_map_u6 <- median(xgenes_u62$u6med)#0.729

xgenes_nv1 <- xgenes4 %>% dplyr::select(gene, novectormed, NO_VECTOR_genecpm) %>% dplyr::filter(NO_VECTOR_genecpm >= 1)
xgenes_nv2 <- xgenes_nv1 %>% dplyr::filter(novectormed<=0.9)
med_map_nv <- median(xgenes_nv2$novectormed)#0.652


#FIGURE 9 SNV major allele proportion
#BOXPLOTS, histograms Figure 9A, B, C
#DEFINE SNV SETS for each knockdown BASED ON KNOCKDOWN + U6M2 control
#filter by each kd individually to retain only specific kd entries. Need to do each kd separately as each set of u6m2 snvs is unique to kd
#create list of separate dataframes
setwd("/Users/.../new_combined_vcf_mMonDom1/")
options(scipen = 0) #use scientfic notation

#filtering for exonic SNVs only 
kd_dip5a <- kd_dip5 %>% filter(exon_start != "NA")

#ALL - including monoallelic (filter for SNV MAP > 0.50)  **changed SNV expression to 0.5 cpm (cf 0.2 cpm) to remove lowly expresses SNVs
kd_dip9_caprin <- kd_dip5a %>% filter(CAPRIN1 >=0.5 | U6M2 >=0.5 ) %>% filter(u6m2_cpm>0.5 | caprin1_cpm>0.5) %>% dplyr::select(-CKAP4, -HNRNPK, -NCL, -NONO, -SYNCRIP, -NO_VECTOR) %>% pivot_longer(cols = c(CAPRIN1, U6M2), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE)
kd_dip9_ckap <- kd_dip5a %>% filter(CKAP4 >=0.5 | U6M2 >=0.5 ) %>% filter(u6m2_cpm>0.5 | ckap4_cpm>0.5) %>% dplyr::select(-CAPRIN1, -HNRNPK, -NCL, -NONO, -SYNCRIP, -NO_VECTOR) %>% pivot_longer(cols = c(CKAP4, U6M2), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE)
kd_dip9_hnrnpk <- kd_dip5a %>% filter(HNRNPK >=0.5 | U6M2 >=0.5 ) %>% filter(u6m2_cpm>0.5 | hnrnpk_cpm>0.5) %>% dplyr::select(-CKAP4, -CAPRIN1, -NCL, -NONO, -SYNCRIP, -NO_VECTOR) %>% pivot_longer(cols = c(HNRNPK, U6M2), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE)
kd_dip9_ncl <- kd_dip5a %>% filter(NCL >=0.5 | U6M2 >=0.5 ) %>% filter(u6m2_cpm>0.5 | ncl_cpm>0.5) %>% dplyr::select(-CKAP4, -HNRNPK, -CAPRIN1, -NONO, -SYNCRIP, -NO_VECTOR) %>% pivot_longer(cols = c(NCL, U6M2), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE)
kd_dip9_nono <- kd_dip5a %>% filter(NONO >=0.5 | U6M2 >=0.5 ) %>% filter(u6m2_cpm>0.5 | nono_cpm>0.5) %>% dplyr::select(-CKAP4, -HNRNPK, -NCL, -CAPRIN1, -SYNCRIP, -NO_VECTOR) %>% pivot_longer(cols = c(NONO, U6M2), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE)
kd_dip9_syncrip <- kd_dip5a %>% filter(SYNCRIP >=0.5 | U6M2 >=0.5 ) %>% filter(u6m2_cpm>0.5 | syncrip_cpm>0.5) %>% dplyr::select(-CKAP4, -HNRNPK, -NCL, -NONO, -CAPRIN1, -NO_VECTOR) %>% pivot_longer(cols = c(SYNCRIP, U6M2), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE)
kd_dip9_no_vector <- kd_dip5a %>% filter(NO_VECTOR >=0.5 | U6M2 >=0.5 ) %>% filter(u6m2_cpm>0.5 | no_vector_cpm>0.5) %>% dplyr::select(-CKAP4, -HNRNPK, -NCL, -NONO, -CAPRIN1, -SYNCRIP) %>% pivot_longer(cols = c(NO_VECTOR, U6M2), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE)

dflist <- list(kd_dip9_caprin, kd_dip9_ckap, kd_dip9_hnrnpk, kd_dip9_ncl, kd_dip9_nono, kd_dip9_syncrip)

num_files = length(dflist)
target <- list("CAPRIN1", "CKAP4", "HNRNPK", "NCL", "NONO", "SYNCRIP")

kd_map_boxplots = list()
kd_map_overlay_histo_A = list()
kd_map_overlay_histo_X = list()
kd_map_counts <- data.frame(kd = character(0), Akdc = numeric(0), Au6c = numeric(0), Xkdc = numeric(0), Xu6c = numeric(0), moodAp = numeric(0), moodXp = numeric(0), iqrAkd = numeric(0), iqrXkd = numeric(0), iqrAu6 = numeric(0), iqrXu6 = numeric(0))
results_WX <- data.frame(kd = character(0), p.value_a = numeric(0), p.value_x = numeric(0), p.value_kd = numeric(0), p.value_u6 = numeric(0))

for (j in 1:6){
  
  kd_dip9 <- dflist[[j]]
  
  kd <- target[[j]]
  
  setA <- kd_dip9 %>% filter(XA == "A")
  setX <- kd_dip9 %>% filter(XA == "X")
  #counts and medians  
  Akd <- setA %>% filter(id == kd) 
  Au6 <- setA %>% filter(id == "U6M2") 
  Xkd <- setX %>% filter(id == kd) 
  Xu6 <- setX %>% filter(id == "U6M2") 
 Akdc <- Akd %>% dplyr::count(id) 
 Au6c <- Au6 %>% dplyr::count(id) 
  Xkdc <- Xkd %>% dplyr::count(id) 
  Xu6c <- Xu6 %>% dplyr::count(id) 
  moodA <- mood.medtest(map ~ id, data = setA)
  moodX <- mood.medtest(map ~ id, data = setX)
  moodAp <- moodA$p.value
  moodXp <- moodX$p.value
  iqrAkd <- IQR(Akd$map)
  iqrXkd <- IQR(Xkd$map)
  iqrAu6 <- IQR(Au6$map)
  iqrXu6 <- IQR(Xu6$map)
  res1 <- c(kd, Akdc$n, Au6c$n, Xkdc$n, Xu6c$n, moodAp, moodXp, iqrAkd, iqrXkd, iqrAu6, iqrXu6)
  kd_map_counts[nrow(kd_map_counts)+1,] <- res1
  
  #overlay histograms of distribution of map in kd v control  
  kd_map_overlay_histo_A[[j]] <- ggplot(setA, aes(x = map, color = id))+
    geom_histogram(bins = 20, alpha = 0.10, position = "identity", aes(fill = id))+
    #geom_density(aes(x = map, y = ..density..), fill = "#F8766D", alpha = 0.5)+
    theme_classic(base_size = 16)+
    labs(subtitle = paste0(target[[j]]))+
   scale_y_continuous(limits = c(0, 12000))+ #for exonic snvs
    scale_color_manual(values = c("dodgerblue3", "orange"))+
    scale_fill_manual(values = c("dodgerblue3", "orange"))+
    xlab("")+
    ylab("")+
    theme(legend.title = element_blank())+
    theme(legend.position = "none")
  print(kd_map_overlay_histo_A[[j]])
  
  kd_map_overlay_histo_X[[j]] <- ggplot(setX, aes(x = map, color = id))+
    geom_histogram(bins = 20, alpha = 0.10, position = "identity", aes(fill = id))+
    theme_classic(base_size = 16)+
    labs(subtitle = paste0(target[[j]]))+
    scale_x_continuous(limits = c(0, 1.05))+
    scale_y_continuous(limits = c(0, 60))+ #for exonic snvs
    scale_color_manual(values = c("dodgerblue3", "orange"))+
    scale_fill_manual(values = c("dodgerblue3", "orange"))+
    xlab("")+
    ylab("")+
    theme(legend.title = element_blank())+
    theme(legend.position = "none")
  print(kd_map_overlay_histo_X[[j]])
  
  #boxplots
  kd_map_boxplots[[j]] <- ggplot(kd_dip9, aes(x = id, y = map))+
    geom_boxplot(notch = TRUE, outliers = FALSE, aes(fill = XA))+ 
    theme_classic(base_size = 20)+
    ylab("")+
    xlab("")+
    scale_y_continuous(breaks = c(0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1.0))+
    coord_cartesian(ylim = c(0.2, 1.2))+
    geom_text(data = Akd %>% dplyr::count(id), aes(y = 0.25, label = n), hjust = 1.5, size = 5)+
    geom_text(data = Xkd %>% dplyr::count(id), aes(y = 0.25, label = n), hjust = -0.5, size = 5)+
    geom_text(data = Au6 %>% dplyr::count(id), aes(y = 0.25, label = n), hjust = 1.5, size = 5)+
    geom_text(data = Xu6 %>% dplyr::count(id), aes(y = 0.25, label = n), hjust = -0.5, size = 5)+
    stat_summary(data = Akd, fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = 1.5, vjust = -1))+
    stat_summary(data = Xkd, fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = -1, vjust = -1))+
    stat_summary(data = Au6, fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = 1.5, vjust = -1))+
    stat_summary(data = Xu6, fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = -1, vjust = -1))+
    theme(legend.title = element_blank())+
    theme(legend.position = "none")
  print(kd_map_boxplots[[j]])
  
  #wilcoxon two-sample tests for X v A and for kd v control  
  WXa <- rstatix::wilcox_test(data = setA, map ~ id, paired = TRUE)
  WXx <- rstatix::wilcox_test(data = setX, map ~ id, paired = TRUE)
  WXkd <- rstatix::wilcox_test(data = kd_dip9 %>% filter(id == kd), map ~ XA, paired = FALSE)
  WXu6 <- rstatix::wilcox_test(data = kd_dip9 %>% filter(id == "U6M2"), map ~ XA, paired = FALSE)
  res2 <- c(kd, WXa$p, WXx$p, WXkd$p, WXu6$p)
  results_WX[nrow(results_WX)+1,] <- res2
  
}  

kd_map_counts <- as_tibble(kd_map_counts)
kd_map_counts$moodAp <- as.numeric(kd_map_counts$moodAp)
kd_map_counts$moodXp <- as.numeric(kd_map_counts$moodXp)
kd_map_counts$moodAp <- format(kd_map_counts$moodAp, scientific =  TRUE, digits = 3)
kd_map_counts$moodXp <- format(kd_map_counts$moodXp, scientific =  TRUE, digits = 3)
write.table(kd_map_counts, file="XA_kd_map_counts_allmono_exon_05.txt", quote=FALSE, sep="\t", row.names = FALSE)

results_WX <- as_tibble(results_WX)
results_WX$p.value_kd <- as.numeric(results_WX$p.value_kd)
results_WX$p.value_u6 <- as.numeric(results_WX$p.value_u6)
results_WX$p.value_a <- as.numeric(results_WX$p.value_a)
results_WX$p.value_x <- as.numeric(results_WX$p.value_x)
results_WX <- results_WX %>% mutate("p.signifa" = case_when(p.value_a <= 0.0001 ~ "****", p.value_a <= 0.001 ~ "***", p.value_a <= 0.01 ~ "**", p.value_a <= 0.05 ~ "*", p.value_a > 0.05 ~ "NS"))
results_WX <- results_WX %>% mutate("p.signifx" = case_when(p.value_x <= 0.0001 ~ "****", p.value_x <= 0.001 ~ "***", p.value_x <= 0.01 ~ "**", p.value_x <= 0.05 ~ "*", p.value_x > 0.05 ~ "NS"))
results_WX <- results_WX %>% mutate("p.signifkd" = case_when(p.value_kd <= 0.0001 ~ "****", p.value_kd <= 0.001 ~ "***", p.value_kd <= 0.01 ~ "**", p.value_kd <= 0.05 ~ "*", p.value_kd > 0.05 ~ "NS"))
results_WX <- results_WX %>% mutate("p.signifu6" = case_when(p.value_u6 <= 0.0001 ~ "****", p.value_u6 <= 0.001 ~ "***", p.value_u6 <= 0.01 ~ "**", p.value_u6 <= 0.05 ~ "*", p.value_u6 > 0.05 ~ "NS"))
write.table(results_WX, file="XA_kd_map_WX_allmono_exon_05.txt", quote=FALSE, sep="\t", row.names = FALSE)

p1 <- as.ggplot(kd_map_boxplots[[1]])
p2 <- as.ggplot(kd_map_boxplots[[2]])
p3 <- as.ggplot(kd_map_boxplots[[3]])
p4 <- as.ggplot(kd_map_boxplots[[4]])
p5 <- as.ggplot(kd_map_boxplots[[5]])
p6 <- as.ggplot(kd_map_boxplots[[6]])

# use patchwork to arrange them together
plots <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3) 
plots
ggsave(file = "XA_kd_map_boxplots_allmono_exon_05.svg", plot = plots,  device = "svg", width = 15, height = 12, units = "in")

pnude <- as.ggplot(kd_map_boxplots[[7]])
pnude
ggsave(file = "XA_kd_map_boxplots_nude_allmono_exon.svg", plot = pnude,  device = "svg", width = 5, height = 4, units = "in")

p7 <- as.ggplot(kd_map_overlay_histo_A[[1]])
p8 <- as.ggplot(kd_map_overlay_histo_A[[2]])
p9 <- as.ggplot(kd_map_overlay_histo_A[[3]])
p10 <- as.ggplot(kd_map_overlay_histo_A[[4]])
p11 <- as.ggplot(kd_map_overlay_histo_A[[5]])
p12 <- as.ggplot(kd_map_overlay_histo_A[[6]])
plots2 <- p7 + p8 + p9 + p10 + p11 + p12 + plot_layout(ncol = 6) 
plots2
ggsave(file = "XA_kd_map_overlay_histo_A_allmono_exon_05.svg", plot = plots2,  device = "svg", width = 18, height = 3.5, units = "in")

pnude2 <- as.ggplot(kd_map_overlay_histo_A[[7]])
pnude2
ggsave(file = "XA_kd_map_histos_A_nude_allmono_exon.svg", plot = pnude2,  device = "svg", width = 5, height = 4, units = "in")

p13 <- as.ggplot(kd_map_overlay_histo_X[[1]])
p14 <- as.ggplot(kd_map_overlay_histo_X[[2]])
p15 <- as.ggplot(kd_map_overlay_histo_X[[3]])
p16 <- as.ggplot(kd_map_overlay_histo_X[[4]])
p17 <- as.ggplot(kd_map_overlay_histo_X[[5]])
p18 <- as.ggplot(kd_map_overlay_histo_X[[6]])
plots3 <- p13 + p14 + p15 + p16 + p17 + p18 + plot_layout(ncol = 6) 
plots3
ggsave(file = "XA_kd_map_overlay_histo_X_allmono_exon_05.svg", plot = plots3,  device = "svg", width = 18, height = 3.5, units = "in")

pnude3 <- as.ggplot(kd_map_overlay_histo_X[[7]])
pnude3
ggsave(file = "XA_kd_map_histos_X_nude_allmono_exon.svg", plot = pnude3,  device = "svg", width = 6, height = 6, units = "in")


#FIGURE 9D, E
#BOOTSTRAPPING TO CLARIFY MEDIAN VALUES - BASED ON https://library.virginia.edu/data/articles/the-wilcoxon-rank-sum-test
#By default the wilcox.test() function will calculate exact p-values if the samples contains less than 50 finite values and there are no ties in the values.
#uses boot package
#The idea is to resample the data (with replacement) many times, say 1000 times, each time taking a difference in medians. 
#We then take the median of those 1000 differences to estimate the difference in medians. 
#We can then find a confidence interval based on our 1000 differences. An easy way is to use the 2.5th and 97.5th percentiles as the upper and lower bounds of a 95% confidence interval.
#d is argument for data
#i is argument to index the data
#function using med.diff will then resample data and return differnece in medians for resampled data
#bootstrap for autosomal - resample (with replacement) 1000x, taking difference in medians. 
#Then take median of those 1000 diff in medians to better estimate diff in median. Then calculate confidence interval
#use boot() function to resample data 1000 times, taking a difference in medians each time, and saving the results into an object called boot.out.
#The boot.out object is a list object. The element named "t" contains the 1000 differences in medians. 
#Taking the median of those values gives us a point estimate of the estimated difference in medians. 

options(scipen = 0)#use scientfic notation

getwd()
dflist <- list(kd_dip9_caprin, kd_dip9_ckap, kd_dip9_hnrnpk, kd_dip9_ncl, kd_dip9_nono, kd_dip9_syncrip)
#dflist <- list(kd_dip9_caprin, kd_dip9_ckap, kd_dip9_hnrnpk, kd_dip9_ncl, kd_dip9_nono, kd_dip9_syncrip, kd_dip9_no_vector)
num_files = length(dflist)
target <- list("CAPRIN1", "CKAP4", "HNRNPK", "NCL", "NONO", "SYNCRIP")
results_boot <- data.frame(kd = character(0), med_med = numeric(0), conf.int1 = numeric(0), conf.int2 = numeric(0))

for (j in 1:6){
  
  kd_dip9 <- dflist[[j]]
  
  kd <- target[[j]]
  
  setA <- kd_dip9 %>% filter(XA == "A")
  #setX <- kd_dip9 %>% filter(XA == "X")

med.diff <- function(d, i) {
  setA <- d[i,] 
  median(setA$map[setA$id == kd]) - 
    median(setA$map[setA$id == "U6M2"])
}

boot.out <- boot(data = setA, statistic = med.diff, R = 1000)
median(boot.out$t) 
med_med <- median(boot.out$t)#- this is median difference of medians
bc <- boot.ci(boot.out, type = "perc") #confidence intervals. based on 1000 bootstrap replicates
conf.int1 <- bc$percent[4]
conf.int2 <- bc$percent[5]
res <- c(kd, med_med, conf.int1, conf.int2)
results_boot[nrow(results_boot)+1,] <- res
}

results_boot$med_med <- as.numeric(results_boot$med_med)
results_boot$conf.int1 <- as.numeric(results_boot$conf.int1)
results_boot$conf.int2 <- as.numeric(results_boot$conf.int2)

write.table(results_boot, file="bootstrapping_median_of_medians_A_05.txt", quote=FALSE, sep="\t", row.names = FALSE)

#write.table(results_boot, file="bootstrapping_median_of_medians_X_05.txt", quote=FALSE, sep="\t", row.names = FALSE)

#plot of bootstrapped medians + confidence intervals
#options(scipen = 999) #don't use sci notation

boot_plot2 <- ggplot(results_boot) +
  geom_point(aes(x = kd, y = med_med), colour = "blue")+
  geom_segment(aes(x = kd, y = conf.int1, yend = conf.int2))+
  geom_hline(yintercept = 0, colour = "blue")+
  theme_classic(base_size = 16)+
  labs(title = "Bootstrapping of medians")+
 # labs(subtitle = "95% confidence interval for median of difference in medians")+
  xlab("")+
  ylab("difference in medians (knockdown - control)")
boot_plot2
ggsave(file = "XA_kd_map_exon_median_bootstrapping_A_05.svg", plot = boot_plot2,  device = "svg", width = 8, height = 4.5, units = "in")
#ggsave(file = "XA_kd_map_exon_median_bootstrapping_X_05.svg", plot = boot_plot2,  device = "svg", width = 8, height = 4.5, units = "in")


#FIGURE 10 + associated Supp figures - scatterplots of exonic SNV expression cf allelic balance
setwd("/Users/.../new_combined_vcf_mMonDom1/")
options(scipen = 0) #use scientfic notation

kd_dip10 <- kd_dip5 %>% filter(exon_start != "NA") #filter for exonic SVNs
#FIGURE 10 SCATTER PLOTS - INDIVIDUAL SNV CPM V INDIVIDUAL SNV MAP
#select for exonic SNVs
kd_caprin <- kd_dip10 %>% filter(CAPRIN1 >=0.5) %>% filter(caprin1_cpm >0.5) %>% mutate(cpm = caprin1_cpm) %>% mutate(kd = CAPRIN1) %>% distinct(pos, gene, .keep_all = TRUE)
kd_ckap4 <- kd_dip10 %>% filter(CKAP4 >=0.5) %>% filter(ckap4_cpm >0.5) %>% mutate(cpm = ckap4_cpm) %>% mutate(kd = CKAP4) %>% distinct(pos, gene, .keep_all = TRUE)
kd_hk <- kd_dip10 %>% filter(HNRNPK >=0.5) %>% filter(hnrnpk_cpm >0.5) %>% mutate(cpm = hnrnpk_cpm) %>% mutate(kd = HNRNPK) %>% distinct(pos, gene, .keep_all = TRUE)
kd_ncl <- kd_dip10 %>% filter(NCL >=0.5) %>% filter(ncl_cpm >0.5) %>% mutate(cpm = ncl_cpm)  %>% mutate(kd = NCL) %>% distinct(pos, gene, .keep_all = TRUE)
kd_nono <- kd_dip10 %>% filter(NONO >=0.5) %>% filter(nono_cpm >0.5) %>% mutate(cpm = nono_cpm)  %>% mutate(kd = NONO) %>% distinct(pos, gene, .keep_all = TRUE)
kd_syncrip <- kd_dip10 %>% filter(SYNCRIP >=0.5) %>% filter(syncrip_cpm >0.5) %>% mutate(cpm = syncrip_cpm)  %>% mutate(kd = SYNCRIP) %>% distinct(pos, gene, .keep_all = TRUE)
kd_u6m2 <- kd_dip10 %>% filter(U6M2 >=0.5) %>% filter(u6m2_cpm >0.5) %>% mutate(cpm = u6m2_cpm) %>% mutate(kd = U6M2) %>% distinct(pos, gene, .keep_all = TRUE)

dflist <- list(kd_caprin, kd_ckap4, kd_hk, kd_ncl, kd_nono, kd_syncrip, kd_u6m2)
num_files = length(dflist)
id <- list("CAPRIN1", "CKAP4", "HNRNPK", "NCL", "NONO", "SYNCRIP", "U6M2")
result_eqn <- data.frame(kd = character(0), lm_eqn_A = character(0), lm_eqn_x = character(0))
cpm_scatters = list()

for (j in 1:7){
  
  kd_dip10gene <- dflist[[j]]
  kd <- id[[j]]
  
  
  kd_dip10gene_A <- kd_dip10gene %>% filter(XA == "A")
  kd_dip10gene_X <- kd_dip10gene %>% filter(XA == "X")
  
  lm_eqn_A <- function(kd_dip10gene_A) {
    m <- lm(log2(cpm) ~ kd, kd_dip10gene_A);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(unname(coef(m)[1]), digits = 2),
                          b = format(unname(coef(m)[2]), digits = 2),
                          r2 = format(summary(m)$r.squared, digits = 2)))
    as.character(as.expression(eq));
  }
  
  lm_eqn_X <- function(kd_dip10gene_X) {
    m <- lm(log2(cpm) ~ kd, kd_dip10gene_X);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(unname(coef(m)[1]), digits = 2),
                          b = format(unname(coef(m)[2]), digits = 2),
                          r2 = format(summary(m)$r.squared, digits = 2)))
    as.character(as.expression(eq));
  }
  
  res <- c(kd, lm_eqn_A(kd_dip10gene_A), lm_eqn_A(kd_dip10gene_X))
  result_eqn[nrow(result_eqn)+1,] <- res
  
  cpm_scatters[[j]] <- ggplot(kd_dip10gene %>% filter(median(kd) >= 0.5), aes(x = kd, y = log2(cpm)))+
    #geom_point(size = 0.3)+
    geom_point(aes(color = XA), size = 0.3, alpha = 0.5)+
    geom_smooth(method = lm, aes(color = XA))+
    theme_classic(base_size = 16) +
    ylab("")+
    xlab("")+
    theme(legend.position = "none")+
    labs(subtitle = paste0(kd))+
    ggpubr::stat_cor(data = kd_dip10gene %>% filter(XA == "A"), label.x = 0.7, label.y = 4.5, method = "pearson", size = 4)+
    ggpubr::stat_cor(data = kd_dip10gene %>% filter(XA == "X"), label.x = 0.7, label.y = 4.1, method = "pearson", size = 4)
  print(cpm_scatters[[j]])
  
}

result_eqn <- as_tibble(result_eqn)

p1 <- as.ggplot(cpm_scatters[[1]])
p2 <- as.ggplot(cpm_scatters[[2]])
p3 <- as.ggplot(cpm_scatters[[3]])
p4 <- as.ggplot(cpm_scatters[[4]])
p5 <- as.ggplot(cpm_scatters[[5]])
p6 <- as.ggplot(cpm_scatters[[6]])
p7 <- as.ggplot(cpm_scatters[[7]])
# use patchwork to arrange them together
# geom_smooth default plots SE with 95% confidence interval
plots <- p1 + p2 + p3 + p4 + p5 + p6 + p7+ plot_layout(ncol = 3) 
plots
ggsave(file = "cpm_scatters_SNVind.tiff", plot = plots,  device = "tiff", width = 12, height = 12, units = "in")


#Supplementary Figure S3-1
#SCATTER PLOTS - GENE RNA-SEQ CPM V GENE MEDIAN SNV MAP
#select for exonic SNVs
# geom_smooth default plots SE with 95% confidence interval

kd_caprin <- kd_dip10 %>% filter(CAPRIN1 >=0.5) %>% filter(caprin1_cpm >0.5) %>% group_by(chromosome, gene) %>% mutate(CAPRIN1 = median(CAPRIN1)) %>% ungroup %>% mutate(kd = CAPRIN1) %>% mutate(cpm = CAPRIN1_genecpm) %>% distinct(gene, .keep_all = TRUE)
kd_ckap4 <- kd_dip10 %>% filter(CKAP4 >=0.5) %>% filter(ckap4_cpm >0.5) %>% group_by(chromosome, gene) %>% mutate(CKAP4 = median(CKAP4)) %>% ungroup %>% mutate(kd = CKAP4) %>% mutate(cpm = CKAP4_genecpm) %>% distinct(gene, .keep_all = TRUE)
kd_hk <- kd_dip10 %>% filter(HNRNPK >=0.5) %>% filter(hnrnpk_cpm >0.5) %>% group_by(chromosome, gene) %>% mutate(HNRNPK = median(HNRNPK)) %>% ungroup %>% mutate(kd = HNRNPK) %>% mutate(cpm = HNRNPK_genecpm) %>% distinct(gene, .keep_all = TRUE)
kd_ncl <- kd_dip10 %>% filter(NCL >=0.5) %>% filter(ncl_cpm >0.5) %>% group_by(chromosome, gene) %>% mutate(NCL = median(NCL)) %>% ungroup %>% mutate(kd = NCL) %>% mutate(cpm = NCL_genecpm) %>% distinct(gene, .keep_all = TRUE)
kd_nono <- kd_dip10 %>% filter(NONO >=0.5) %>% filter(nono_cpm >0.5) %>% group_by(chromosome, gene) %>% mutate(NONO = median(NONO)) %>% ungroup %>% mutate(kd = NONO) %>% mutate(cpm = NONO_genecpm) %>% distinct(gene, .keep_all = TRUE)
kd_syncrip <- kd_dip10 %>% filter(SYNCRIP >=0.5) %>% filter(syncrip_cpm >0.5) %>% group_by(chromosome, gene) %>% mutate(SYNCRIP = median(SYNCRIP)) %>% ungroup %>% mutate(kd = SYNCRIP) %>% mutate(cpm = SYNCRIP_genecpm) %>% distinct(gene, .keep_all = TRUE)
kd_u6m2 <- kd_dip10 %>% filter(U6M2 >=0.5) %>% filter(u6m2_cpm >0.5) %>% group_by(chromosome, gene) %>% mutate(U6M2 = median(U6M2)) %>% mutate(kd = U6M2) %>% ungroup %>% mutate(cpm = U6M2_genecpm) %>% distinct(gene, .keep_all = TRUE)

dflist <- list(kd_caprin, kd_ckap4, kd_hk, kd_ncl, kd_nono, kd_syncrip, kd_u6m2)
num_files = length(dflist)
id <- list("CAPRIN1", "CKAP4", "HNRNPK", "NCL", "NONO", "SYNCRIP", "U6M2")
result_eqn <- data.frame(kd = character(0), lm_eqn_A = character(0), lm_eqn_x = character(0))
cpm_scatters = list()

for (j in 1:7){
  
  kd_dip10gene <- dflist[[j]]
  kd <- id[[j]]
  
  kd_dip10gene_A <- kd_dip10gene %>% filter(XA == "A")
  kd_dip10gene_X <- kd_dip10gene %>% filter(XA == "X")
  
  lm_eqn_A <- function(kd_dip10gene_A) {
    m <- lm(log2(cpm) ~ kd, kd_dip10gene_A);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(unname(coef(m)[1]), digits = 2),
                          b = format(unname(coef(m)[2]), digits = 2),
                          r2 = format(summary(m)$r.squared, digits = 2)))
    as.character(as.expression(eq));
  }
  
  lm_eqn_X <- function(kd_dip10gene_X) {
    m <- lm(log2(cpm) ~ kd, kd_dip10gene_X);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(unname(coef(m)[1]), digits = 2),
                          b = format(unname(coef(m)[2]), digits = 2),
                          r2 = format(summary(m)$r.squared, digits = 2)))
    as.character(as.expression(eq));
  }
  
  
  res <- c(kd, lm_eqn_A(kd_dip10gene_A), lm_eqn_A(kd_dip10gene_X))
  result_eqn[nrow(result_eqn)+1,] <- res
  
  cpm_scatters[[j]] <- ggplot(kd_dip10gene, aes(x = kd, y = log2(cpm)))+
    #geom_point(size = 0.3)+
    geom_point(aes(color = XA), size = 0.3, alpha = 0.5)+
    geom_smooth(method = lm, aes(color = XA))+
    theme_classic(base_size = 16) +
    ylab("")+
    xlab("")+
    theme(legend.position = "none")+
    labs(subtitle = paste0(kd))+
    ggpubr::stat_cor(data = kd_dip10gene %>% filter(XA == "A"), label.x = 0.7, label.y = 12, method = "pearson", size = 4)+
    ggpubr::stat_cor(data = kd_dip10gene %>% filter(XA == "X"), label.x = 0.7, label.y = 11, method = "pearson", size = 4)
  
  print(cpm_scatters[[j]])
  
}

result_eqn <- as_tibble(result_eqn)

p1 <- as.ggplot(cpm_scatters[[1]])
p2 <- as.ggplot(cpm_scatters[[2]])
p3 <- as.ggplot(cpm_scatters[[3]])
p4 <- as.ggplot(cpm_scatters[[4]])
p5 <- as.ggplot(cpm_scatters[[5]])
p6 <- as.ggplot(cpm_scatters[[6]])
p7 <- as.ggplot(cpm_scatters[[7]])

# use patchwork to arrange them together

plots <- p1 + p2 + p3 + p4 + p5 + p6 + p7+ plot_layout(ncol = 3) 
plots
ggsave(file = "cpm_scatters_gene.tiff", plot = plots,  device = "tiff", width = 12, height = 12, units = "in")


#FIGURE 11 Gene expression (RNA-seq) based on allelic balance
#PLOT EXPRESSION (MPILEUP MAPPED READS) BASED ON SKEW/BALANCE IN SNV EXPRESSION BUT PLOT LOG2 EXPRESSION (RATHER THAN RATIO) TO CAPTURE SNVS WHERE EXPRESSION IN U6M2 = ZERO
#Figure 3-11 plots RNA-deq gene expression v gene median SNV MAP
#Supp Figure S3-2 plots individual SNV expression v individual SNV MAP

setwd("/Users/kimmcintyre/.../new_combined_vcf_mMonDom1/")

#read in normalised cpm counts ... caprin_cpm , etc - replace CAPRIN1_qr, etc (not normalised)

#USE this for RNA-seq cpm expression - Figure 3-11
kd_dip11 <- kd_dip5 %>% dplyr::select(-CAPRIN1_qr, -CKAP4_qr, -HNRNPK_qr, -NCL_qr, -NONO_qr, -SYNCRIP_qr, -U6M2_qr, -NO_VECTOR_qr, -NO_VECTOR_genecpm)
#rename gene cpm
colnames(kd_dip11)[18:24] <- c("CAPRIN1_qr", "CKAP4_qr", "HNRNPK_qr", "NCL_qr", "NONO_qr", "SYNCRIP_qr", "U6M2_qr")

kd_dip11 <- kd_dip11 %>% filter(exon_start != "NA")

kd_dip11g_caprin <- kd_dip11 %>% filter(U6M2 >=0.5 | CAPRIN1 >=0.5) %>% filter(u6m2_cpm >0.5 | caprin1_cpm >0.5) %>% 
  group_by(chromosome, gene) %>% mutate(CAPRIN1 = median(CAPRIN1)) %>% mutate(U6M2 = median(U6M2)) %>% ungroup %>%
  dplyr::select(-CKAP4, -HNRNPK, -NCL, -NONO, -SYNCRIP) %>% dplyr::select(-CKAP4_qr, -HNRNPK_qr, -NCL_qr, -NONO_qr, -SYNCRIP_qr) %>%
  pivot_longer(cols = c(CAPRIN1, U6M2), names_to = "id", values_to = "map") %>% pivot_longer(cols = c(CAPRIN1_qr, U6M2_qr), names_to = "id2", values_to = "qr")

kd_dip11g_ckap <- kd_dip11 %>% filter(U6M2 >=0.5 | CKAP4 >=0.5) %>% filter(u6m2_cpm >0.5 | ckap4_cpm >0.5) %>%
  group_by(chromosome, gene) %>% mutate(CKAP4 = median(CKAP4)) %>% mutate(U6M2 = median(U6M2)) %>% ungroup %>%
  dplyr::select(-CAPRIN1, -HNRNPK, -NCL, -NONO, -SYNCRIP) %>% dplyr::select(-CAPRIN1_qr, -HNRNPK_qr, -NCL_qr, -NONO_qr, -SYNCRIP_qr) %>%
  pivot_longer(cols = c(CKAP4, U6M2), names_to = "id", values_to = "map") %>% pivot_longer(cols = c(CKAP4_qr, U6M2_qr), names_to = "id2", values_to = "qr")

kd_dip11g_hnrnpk <- kd_dip11 %>% filter(U6M2 >=0.5 | HNRNPK >=0.5) %>% filter(u6m2_cpm >0.5 | hnrnpk_cpm >0.5) %>%
  group_by(chromosome, gene) %>% mutate(HNRNPK = median(HNRNPK)) %>% mutate(U6M2 = median(U6M2)) %>% ungroup %>%
  dplyr::select(-CKAP4, -CAPRIN1, -NCL, -NONO, -SYNCRIP) %>% dplyr::select(-CKAP4_qr, -CAPRIN1_qr, -NCL_qr, -NONO_qr, -SYNCRIP_qr) %>%
  pivot_longer(cols = c(HNRNPK, U6M2), names_to = "id", values_to = "map") %>% pivot_longer(cols = c(HNRNPK_qr, U6M2_qr), names_to = "id2", values_to = "qr") 

kd_dip11g_ncl <- kd_dip11 %>% filter(U6M2 >=0.5 | NCL >=0.5) %>% filter(u6m2_cpm >0.5 | ncl_cpm >0.5) %>%
  group_by(chromosome, gene) %>% mutate(NCL = median(NCL)) %>% mutate(U6M2 = median(U6M2)) %>% ungroup %>%
  dplyr::select(-CKAP4, -HNRNPK, -CAPRIN1, -NONO, -SYNCRIP) %>% dplyr::select(-CKAP4_qr, -HNRNPK_qr, -CAPRIN1_qr, -NONO_qr, -SYNCRIP_qr) %>%
  pivot_longer(cols = c(NCL, U6M2), names_to = "id", values_to = "map") %>% pivot_longer(cols = c(NCL_qr, U6M2_qr), names_to = "id2", values_to = "qr")

kd_dip11g_nono <- kd_dip11 %>% filter(U6M2 >=0.5 | NONO >=0.5) %>% filter(u6m2_cpm >0.5 | nono_cpm >0.5) %>% 
  group_by(chromosome, gene) %>% mutate(NONO = median(NONO)) %>% mutate(U6M2 = median(U6M2)) %>% ungroup %>%
  dplyr::select(-CKAP4, -HNRNPK, -NCL, -CAPRIN1, -SYNCRIP) %>% dplyr::select(-CKAP4_qr, -HNRNPK_qr, -NCL_qr, -CAPRIN1_qr, -SYNCRIP_qr) %>%
  pivot_longer(cols = c(NONO, U6M2), names_to = "id", values_to = "map") %>% pivot_longer(cols = c(NONO_qr, U6M2_qr), names_to = "id2", values_to = "qr") 

kd_dip11g_syncrip <- kd_dip11 %>% filter(U6M2 >=0.5 | SYNCRIP >=0.5) %>% filter(u6m2_cpm >0.5 | syncrip_cpm >0.5) %>% 
  group_by(chromosome, gene) %>% mutate(SYNCRIP = median(SYNCRIP)) %>% mutate(U6M2 = median(U6M2)) %>% ungroup %>%
  dplyr::select(-CKAP4, -HNRNPK, -NCL, -NONO, -CAPRIN1) %>% dplyr::select(-CKAP4_qr, -HNRNPK_qr, -NCL_qr, -NONO_qr, -CAPRIN1_qr) %>%
  pivot_longer(cols = c(SYNCRIP, U6M2), names_to = "id", values_to = "map") %>% pivot_longer(cols = c(SYNCRIP_qr, U6M2_qr), names_to = "id2", values_to = "qr")

#AUTOSOMAL
dflist <- list(kd_dip11g_caprin, kd_dip11g_ckap, kd_dip11g_hnrnpk, kd_dip11g_ncl, kd_dip11g_nono, kd_dip11g_syncrip) #for GENE MEDIAN
id <- list("CAPRIN1", "CKAP4", "HNRNPK", "NCL", "NONO", "SYNCRIP")
id2 <- list("CAPRIN1_qr", "CKAP4_qr", "HNRNPK_qr", "NCL_qr", "NONO_qr", "SYNCRIP_qr")
exp_plots = list()
exp_binned_histo = list()
category <- list("balanced", "inter", "skewed")
result_WX <- data.frame(kd = character(0), cat = character(0), medikd = numeric(0), mediu6 = numeric(0), p.value = numeric(0))
result_KW <- data.frame(kd = character(0), p.valuekd = numeric(0), p.valueu6 = numeric(0))
result_DT <- data.frame(kd = character(0), kd = character(0), group1 = character(0), group2 = character(0), n1 = numeric(0), n1 = numeric(0), statistic = numeric(0), p = numeric(0), p.adj = numeric(0), p.adj.signif = character(0))

#FOR AUTOSOMAL  
for (j in 1:6){
  kd_dip11s <- dflist[[j]]
  kd <- id[[j]]
  kd2 <- id2[[j]]
  
  kd_dip11_ex <- kd_dip11s %>% mutate(cat = (case_when(map >= 0.50 & map < 0.65 ~ "balanced", 
                                                       (map >= 0.65 & map <= 0.85) ~ "inter",
                                                       map > 0.85 ~ "skewed"))) %>% mutate(log2qr = log2(qr)) %>% ungroup()
  
  kd_dip11_ex$log2qr %>% as.numeric(kd_dip11_ex$log2qr)
  
  #FOR gene medians
  specA <- kd_dip11_ex %>% filter(cat != "NA") %>% filter(XA == "A") %>% filter(exon_start != "NA") %>% filter(gene != "NA") %>% filter(log2qr != "-Inf") %>% distinct(chromosome, gene, id, id2, .keep_all = TRUE)
  
  for (k in category){   
    data <- specA %>% filter(id != "U6M2") %>% dplyr::filter(cat == k)
    WX <- rstatix::pairwise_wilcox_test(data, log2qr ~ id2, paired = FALSE)
    datakd <- data %>% filter(id2 != "U6M2_qr")
    datau6 <- data %>% filter(id2 == "U6M2_qr")
    medi_kd <- median(datakd$log2qr)
    medi_u6 <- median(datau6$log2qr)
    res <- c(kd, k, medi_kd, medi_u6, WX$p)
    result_WX[nrow(result_WX)+1,] <- res
  }
  
  exp_plots[[j]] = ggplot(specA %>% filter(id != "U6M2"), aes(x = cat, y = log2qr))+
    geom_boxplot(aes(fill = id2), notch = TRUE, outliers = FALSE, alpha = 0.8)+
    theme_classic(base_size = 16)+
    ylab("")+
    xlab("")+
    labs(subtitle = paste0(kd))+
    geom_text(data = specA %>% filter(cat == "balanced") %>% filter(id != "U6M2") %>% filter(id2 == kd2), stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=1, size=3)+
    geom_text(data = specA %>% filter(cat == "balanced") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), stat = "count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=-0.5, size=3)+
    geom_text(data = specA %>% filter(cat == "inter") %>% filter(id != "U6M2") %>% filter(id2 == kd2), stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=1, size=3)+
    geom_text(data = specA %>% filter(cat == "inter") %>% filter(id != "U6M2")  %>% filter(id2 == "U6M2_qr"), stat = "count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=-0.5, size=3)+
    geom_text(data = specA %>% filter(cat == "skewed") %>% filter(id != "U6M2") %>% filter(id2 == kd2), stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=1, size=3)+
    geom_text(data = specA %>% filter(cat == "skewed") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), stat = "count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=-0.5, size=3)+
    stat_summary(data = specA %>% filter(cat == "balanced") %>% filter(id != "U6M2") %>% filter(id2 == kd2), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = 1, vjust = -0.5))+
    stat_summary(data = specA %>% filter(cat == "balanced") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = -0.5, vjust = -0.5))+
    stat_summary(data = specA %>% filter(cat == "inter") %>% filter(id != "U6M2") %>% filter(id2 == kd2), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = 1, vjust = -0.5))+
    stat_summary(data = specA %>% filter(cat == "inter") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = -0.5, vjust = -0.5))+
    stat_summary(data = specA %>% filter(cat == "skewed") %>% filter(id != "U6M2") %>% filter(id2 == kd2), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = 1, vjust = -0.5))+
    stat_summary(data = specA %>% filter(cat == "skewed") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = -0.5, vjust = -0.5))+
    theme(legend.position="none")
  print(exp_plots[[j]])
  
  exp_binned_histo[[j]] = ggplot()+
    geom_density(data = specA %>% filter(id != "U6M2"), aes(x = map, y = ..density..), fill = "#F8766D", alpha = 0.5)+
    geom_density(data = specA %>% filter(id == "U6M2"), aes(x = map, y = ..density..), fill = "#00BFC4", alpha = 0.5)+
    theme_classic(base_size = 16)+
    scale_x_continuous(limits = c(0.5, 1))+
    xlab("")+
    ylab("")+ #gene number
    theme(legend.position = "none")+
    labs(subtitle = paste0(kd))
  print(exp_binned_histo[[j]])
  
  KWkd <- rstatix::kruskal_test(data = specA %>% filter(id != "U6M2") %>% filter(id2 != "U6M2_qr"), log2qr ~ cat)
  KWu6 <- rstatix::kruskal_test(data = specA %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), log2qr ~ cat)
  KWkd$p <- format(KWkd$p, scientific = TRUE, digits = 3)
  KWu6$p <- format(KWu6$p, scientific = TRUE, digits = 3)
  res <- c(kd, KWkd$p, KWu6$p)
  result_KW[nrow(result_KW)+1,] <- res
  
  DTkd <- rstatix::dunn_test(data = specA %>% filter(id != "U6M2") %>% filter(id2 != "U6M2_qr"), log2qr ~ cat, p.adjust.method = "holm")
  DTu6 <- rstatix::dunn_test(data = specA %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), log2qr ~ cat, p.adjust.method = "holm")
  res1 <- c(kd, kd2, DTkd[1,2:9])
  res2 <- c(kd, kd2, DTkd[2,2:9])
  res3 <- c(kd, kd2, DTkd[3,2:9])
  res4 <- c(kd, "U6M2_qr", DTu6[1,2:9])
  res5 <- c(kd, "U6M2_qr", DTu6[2,2:9])
  res6 <- c(kd, "U6M2_qr", DTu6[3,2:9])
  result_DT[nrow(result_DT)+1,] <- res1
  result_DT[nrow(result_DT)+1,] <- res2
  result_DT[nrow(result_DT)+1,] <- res3
  result_DT[nrow(result_DT)+1,] <- res4
  result_DT[nrow(result_DT)+1,] <- res5
  result_DT[nrow(result_DT)+1,] <- res6
  
}

p1 <- as.ggplot(exp_plots[[1]])
p2 <- as.ggplot(exp_plots[[2]])
p3 <- as.ggplot(exp_plots[[3]])
p4 <- as.ggplot(exp_plots[[4]])
p5 <- as.ggplot(exp_plots[[5]])
p6 <- as.ggplot(exp_plots[[6]])

plots <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3) 
plots
ggsave(file = "exp_binned_plots_exon_A_gene_rnaseq_05.svg", plot = plots,  device = "svg", width = 15, height = 12, units = "in")

p7 <- as.ggplot(exp_binned_histo[[1]])
p8 <- as.ggplot(exp_binned_histo[[2]])
p9 <- as.ggplot(exp_binned_histo[[3]])
p10 <- as.ggplot(exp_binned_histo[[4]])
p11 <- as.ggplot(exp_binned_histo[[5]])
p12 <- as.ggplot(exp_binned_histo[[6]])

plots <- p7 + p8 + p9 + p10 + p11 + p12 + plot_layout(ncol = 3)
plots
ggsave(file = "exp_binned_histos_exon_A_gene_rnaseq_05.svg", plot = plots,  device = "svg", width = 13, height = 7, units = "in")

result_WX <- as_tibble(result_WX)
result_WX$p.value <- as.numeric(result_WX$p.value)
result_WX$medikd <- as.numeric(result_WX$medikd)
result_WX$mediu6 <- as.numeric(result_WX$mediu6)
result_WX <- result_WX %>% mutate("p.signif" = case_when(p.value <= 0.0001 ~ "****", p.value <= 0.001 ~ "***", p.value <= 0.01 ~ "**", p.value <= 0.05 ~ "*", p.value > 0.05 ~ "NS"))
result_WX$p.value <- format(result_WX$p.value, scientific = TRUE, digits = 3)
result_WX$medikd <- format(result_WX$medikd, scientific = TRUE, digits = 5)
result_WX$mediu6 <- format(result_WX$mediu6, scientific = TRUE, digits = 5)
write.table(result_WX, file="Wilcox_test_exp_binned_exon_A_gene_rnaseq_05.txt", quote=FALSE, sep="\t", row.names = FALSE)

result_KW <- as_tibble(result_KW)
result_KW$p.valuekd <- as.numeric(result_KW$p.valuekd)
result_KW$p.valueu6 <- as.numeric(result_KW$p.valueu6)
result_KW$p.valuekd <- format(result_KW$p.valuekd, scientific = TRUE, digits = 3)
result_KW$p.valueu6 <- format(result_KW$p.valueu6, scientific = TRUE, digits = 3)
write.table(result_KW, file="KW_test_exp_binned_exon_A_gene_rnaseq_05.txt", quote=FALSE, sep="\t", row.names = FALSE)

result_DT <- as_tibble(result_DT)
write.table(result_DT, file="DT_test_exp_binned_exon_A_gene_rnaseq_05.txt", quote=FALSE, sep="\t", row.names = FALSE)


#FOR X chromosome
dflist <- list(kd_dip11g_caprin, kd_dip11g_ckap, kd_dip11g_hnrnpk, kd_dip11g_ncl, kd_dip11g_nono, kd_dip11g_syncrip) #for GENE MEDIAN
id <- list("CAPRIN1", "CKAP4", "HNRNPK", "NCL", "NONO", "SYNCRIP")
id2 <- list("CAPRIN1_qr", "CKAP4_qr", "HNRNPK_qr", "NCL_qr", "NONO_qr", "SYNCRIP_qr")
exp_plots = list()
exp_binned_histo = list()
category <- list("balanced", "skewed")
result_WX <- data.frame(kd = character(0), cat = character(0), medikd = numeric(0), mediu6 = numeric(0), p.value = numeric(0))
result_KW <- data.frame(kd = character(0), p.valuekd = numeric(0), p.valueu6 = numeric(0))
result_DT <- data.frame(kd = character(0), kd = character(0), group1 = character(0), group2 = character(0), n1 = numeric(0), n1 = numeric(0), statistic = numeric(0), p = numeric(0), p.adj = numeric(0), p.adj.signif = character(0))

for (j in 1:6){
  kd_dip11s <- dflist[[j]]
  kd <- id[[j]]
  kd2 <- id2[[j]]
  
  kd_dip11_ex <- kd_dip11s %>% mutate(cat = (case_when(map >= 0.50 & map <= 0.80 ~ "balanced", map > 0.80 ~ "skewed"))) %>% mutate(log2qr = log2(qr)) %>% ungroup()
  
  kd_dip11_ex$log2qr %>% as.numeric(kd_dip11_ex$log2qr)
  
  #FOR gene medians
  specA <- kd_dip11_ex %>% filter(cat != "NA") %>% filter(XA == "X") %>% filter(exon_start != "NA") %>% filter(gene != "NA") %>% filter(log2qr != "-Inf") %>% distinct(chromosome, gene, id, id2, .keep_all = TRUE)
  
  for (k in category){   
    data <- specA %>% filter(id != "U6M2") %>% dplyr::filter(cat == k)
    WX <- rstatix::pairwise_wilcox_test(data, log2qr ~ id2, paired = FALSE)
    datakd <- data %>% filter(id2 != "U6M2_qr")
    datau6 <- data %>% filter(id2 == "U6M2_qr")
    medi_kd <- median(datakd$log2qr)
    medi_u6 <- median(datau6$log2qr)
    res <- c(kd, k, medi_kd, medi_u6, WX$p)
    result_WX[nrow(result_WX)+1,] <- res
  }
  
  exp_plots[[j]] = ggplot(specA %>% filter(id != "U6M2"), aes(x = cat, y = log2qr))+
    geom_boxplot(aes(fill = id2), notch = TRUE, outliers = FALSE, alpha = 0.8)+
    theme_classic(base_size = 16)+
    ylab("")+
    xlab("")+
    labs(subtitle = paste0(kd))+
    geom_text(data = specA %>% filter(cat == "balanced") %>% filter(id != "U6M2") %>% filter(id2 == kd2), stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=1, size=3)+
    geom_text(data = specA %>% filter(cat == "balanced") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), stat = "count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=-0.5, size=3)+
    geom_text(data = specA %>% filter(cat == "skewed") %>% filter(id != "U6M2") %>% filter(id2 == kd2), stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=1, size=3)+
    geom_text(data = specA %>% filter(cat == "skewed") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), stat = "count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=-0.5, size=3)+
    stat_summary(data = specA %>% filter(cat == "balanced") %>% filter(id != "U6M2") %>% filter(id2 == kd2), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = 1, vjust = -0.5))+
    stat_summary(data = specA %>% filter(cat == "balanced") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = -0.5, vjust = -0.5))+
    stat_summary(data = specA %>% filter(cat == "skewed") %>% filter(id != "U6M2") %>% filter(id2 == kd2), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = 1, vjust = -0.5))+
    stat_summary(data = specA %>% filter(cat == "skewed") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = -0.5, vjust = -0.5))+
    theme(legend.position="none")
  print(exp_plots[[j]])
  
  exp_binned_histo[[j]] = ggplot()+
    geom_density(data = specA %>% filter(id != "U6M2"), aes(x = map, y = ..density..), fill = "#F8766D", alpha = 0.5)+
    geom_density(data = specA %>% filter(id == "U6M2"), aes(x = map, y = ..density..), fill = "#00BFC4", alpha = 0.5)+
    theme_classic(base_size = 16)+
    scale_x_continuous(limits = c(0.5, 1))+
    xlab("")+
    ylab("")+ #gene number
    theme(legend.position = "none")+
    labs(subtitle = paste0(kd))
  print(exp_binned_histo[[j]])
  
  KWkd <- rstatix::kruskal_test(data = specA %>% filter(id != "U6M2") %>% filter(id2 != "U6M2_qr"), log2qr ~ cat)
  KWu6 <- rstatix::kruskal_test(data = specA %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), log2qr ~ cat)
  KWkd$p <- format(KWkd$p, scientific = TRUE, digits = 3)
  KWu6$p <- format(KWu6$p, scientific = TRUE, digits = 3)
  res <- c(kd, KWkd$p, KWu6$p)
  result_KW[nrow(result_KW)+1,] <- res
  
  DTkd <- rstatix::dunn_test(data = specA %>% filter(id != "U6M2") %>% filter(id2 != "U6M2_qr"), log2qr ~ cat, p.adjust.method = "holm")
  DTu6 <- rstatix::dunn_test(data = specA %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), log2qr ~ cat, p.adjust.method = "holm")
  res1 <- c(kd, kd2, DTkd[1,2:9])
  res2 <- c(kd, kd2, DTkd[2,2:9])
  res3 <- c(kd, "U6M2_qr", DTu6[1,2:9])
  res4 <- c(kd, "U6M2_qr", DTu6[2,2:9])
  result_DT[nrow(result_DT)+1,] <- res1
  result_DT[nrow(result_DT)+1,] <- res2
  result_DT[nrow(result_DT)+1,] <- res3
  result_DT[nrow(result_DT)+1,] <- res4
  
  
}

p1 <- as.ggplot(exp_plots[[1]])
p2 <- as.ggplot(exp_plots[[2]])
p3 <- as.ggplot(exp_plots[[3]])
p4 <- as.ggplot(exp_plots[[4]])
p5 <- as.ggplot(exp_plots[[5]])
p6 <- as.ggplot(exp_plots[[6]])

plots <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3) #+ plot_annotation(title = "Autosomal specific SNVs \nGene expression by allelic balance")
plots
ggsave(file = "exp_binned_plots_exon_X_gene_rnaseq_05.svg", plot = plots,  device = "svg", width = 10, height = 8, units = "in")

p7 <- as.ggplot(exp_binned_histo[[1]])
p8 <- as.ggplot(exp_binned_histo[[2]])
p9 <- as.ggplot(exp_binned_histo[[3]])
p10 <- as.ggplot(exp_binned_histo[[4]])
p11 <- as.ggplot(exp_binned_histo[[5]])
p12 <- as.ggplot(exp_binned_histo[[6]])

plots <- p7 + p8 + p9 + p10 + p11 + p12 + plot_layout(ncol = 3)# + plot_annotation(title = "Autosomal specific SNVs \nGene expression by allelic balance")
plots
ggsave(file = "exp_binned_histos_exon_X_gene_rnaseq_05.svg", plot = plots,  device = "svg", width = 10, height = 7, units = "in")

result_WX <- as_tibble(result_WX)
result_WX$p.value <- as.numeric(result_WX$p.value)
result_WX$medikd <- as.numeric(result_WX$medikd)
result_WX$mediu6 <- as.numeric(result_WX$mediu6)
result_WX <- result_WX %>% mutate("p.signif" = case_when(p.value <= 0.0001 ~ "****", p.value <= 0.001 ~ "***", p.value <= 0.01 ~ "**", p.value <= 0.05 ~ "*", p.value > 0.05 ~ "NS"))
result_WX$p.value <- format(result_WX$p.value, scientific = TRUE, digits = 3)
result_WX$medikd <- format(result_WX$medikd, scientific = TRUE, digits = 5)
result_WX$mediu6 <- format(result_WX$mediu6, scientific = TRUE, digits = 5)
write.table(result_WX, file="Wilcox_test_exp_binned_exon_X_gene_rnaseq_05.txt", quote=FALSE, sep="\t", row.names = FALSE)

result_KW <- as_tibble(result_KW)
result_KW$p.valuekd <- as.numeric(result_KW$p.valuekd)
result_KW$p.valueu6 <- as.numeric(result_KW$p.valueu6)
result_KW$p.valuekd <- format(result_KW$p.valuekd, scientific = TRUE, digits = 3)
result_KW$p.valueu6 <- format(result_KW$p.valueu6, scientific = TRUE, digits = 3)
write.table(result_KW, file="KW_test_exp_binned_intron_X_snvcpm_05.txt", quote=FALSE, sep="\t", row.names = FALSE)

result_DT <- as_tibble(result_DT)
write.table(result_DT, file="DT_test_exp_binned_exon_X_gene_rnaseq_05.txt", quote=FALSE, sep="\t", row.names = FALSE)


#SUPP FIGURE S3-2
#USE this for individual SNV analysis 
kd_dip11 <- kd_dip5 %>% dplyr::select(-CAPRIN1_genecpm, -CKAP4_genecpm, -HNRNPK_genecpm, -NCL_genecpm, -NONO_genecpm, -SYNCRIP_genecpm, -U6M2_genecpm, -NO_VECTOR_genecpm, -CAPRIN1_qr, -CKAP4_qr, -HNRNPK_qr, -NCL_qr, -NONO_qr, -SYNCRIP_qr, -U6M2_qr, -NO_VECTOR_qr, -no_vector_cpm)
colnames(kd_dip11)[18:24] <- c("CAPRIN1_qr", "CKAP4_qr", "HNRNPK_qr", "NCL_qr", "NONO_qr", "SYNCRIP_qr", "U6M2_qr")

kd_dip11_caprin <- kd_dip11 %>% filter(U6M2 >=0.5 | CAPRIN1 >=0.5) %>% filter(U6M2_qr >0.5 | CAPRIN1_qr >0.5) %>% 
  dplyr::select(-CKAP4, -HNRNPK, -NCL, -NONO, -SYNCRIP) %>% dplyr::select(-CKAP4_qr, -HNRNPK_qr, -NCL_qr, -NONO_qr, -SYNCRIP_qr) %>%
  pivot_longer(cols = c(CAPRIN1, U6M2), names_to = "id", values_to = "map") %>% pivot_longer(cols = c(CAPRIN1_qr, U6M2_qr), names_to = "id2", values_to = "qr")

kd_dip11_ckap <- kd_dip11 %>% filter(U6M2 >=0.5 | CKAP4 >=0.5) %>% filter(U6M2_qr >0.5 | CKAP4_qr >0.5) %>% 
  dplyr::select(-CAPRIN1, -HNRNPK, -NCL, -NONO, -SYNCRIP) %>% dplyr::select(-CAPRIN1_qr, -HNRNPK_qr, -NCL_qr, -NONO_qr, -SYNCRIP_qr) %>%
  pivot_longer(cols = c(CKAP4, U6M2), names_to = "id", values_to = "map") %>% pivot_longer(cols = c(CKAP4_qr, U6M2_qr), names_to = "id2", values_to = "qr")

kd_dip11_hnrnpk <- kd_dip11 %>% filter(U6M2 >=0.5 | HNRNPK >=0.5) %>% filter(U6M2_qr >0.5 | HNRNPK_qr >0.5) %>% 
  dplyr::select(-CKAP4, -CAPRIN1, -NCL, -NONO, -SYNCRIP) %>% dplyr::select(-CKAP4_qr, -CAPRIN1_qr, -NCL_qr, -NONO_qr, -SYNCRIP_qr) %>%
  pivot_longer(cols = c(HNRNPK, U6M2), names_to = "id", values_to = "map") %>% pivot_longer(cols = c(HNRNPK_qr, U6M2_qr), names_to = "id2", values_to = "qr") 

kd_dip11_ncl <- kd_dip11 %>% filter(U6M2 >=0.5 | NCL >=0.5) %>% filter(U6M2_qr >0.5 | NCL_qr >0.5) %>% 
  dplyr::select(-CKAP4, -HNRNPK, -CAPRIN1, -NONO, -SYNCRIP) %>% dplyr::select(-CKAP4_qr, -HNRNPK_qr, -CAPRIN1_qr, -NONO_qr, -SYNCRIP_qr) %>%
  pivot_longer(cols = c(NCL, U6M2), names_to = "id", values_to = "map") %>% pivot_longer(cols = c(NCL_qr, U6M2_qr), names_to = "id2", values_to = "qr")

kd_dip11_nono <- kd_dip11 %>% filter(U6M2 >=0.5 | NONO >=0.5) %>% filter(U6M2_qr >0.5 | NONO_qr >0.5) %>% 
  dplyr::select(-CKAP4, -HNRNPK, -NCL, -CAPRIN1, -SYNCRIP) %>% dplyr::select(-CKAP4_qr, -HNRNPK_qr, -NCL_qr, -CAPRIN1_qr, -SYNCRIP_qr) %>%
  pivot_longer(cols = c(NONO, U6M2), names_to = "id", values_to = "map") %>% pivot_longer(cols = c(NONO_qr, U6M2_qr), names_to = "id2", values_to = "qr") 

kd_dip11_syncrip <- kd_dip11 %>% filter(U6M2 >=0.5 | SYNCRIP >=0.5) %>% filter(U6M2_qr >0.5 | SYNCRIP_qr >0.5) %>% 
  dplyr::select(-CKAP4, -HNRNPK, -NCL, -NONO, -CAPRIN1) %>% dplyr::select(-CKAP4_qr, -HNRNPK_qr, -NCL_qr, -NONO_qr, -CAPRIN1_qr) %>%
  pivot_longer(cols = c(SYNCRIP, U6M2), names_to = "id", values_to = "map") %>% pivot_longer(cols = c(SYNCRIP_qr, U6M2_qr), names_to = "id2", values_to = "qr")


#AUTOSOMAL
dflist <- list(kd_dip11_caprin, kd_dip11_ckap, kd_dip11_hnrnpk, kd_dip11_ncl, kd_dip11_nono, kd_dip11_syncrip) #for SNVs
id <- list("CAPRIN1", "CKAP4", "HNRNPK", "NCL", "NONO", "SYNCRIP")
id2 <- list("CAPRIN1_qr", "CKAP4_qr", "HNRNPK_qr", "NCL_qr", "NONO_qr", "SYNCRIP_qr")
exp_plots = list()
exp_binned_histo = list()
category <- list("balanced", "inter", "skewed")
result_WX <- data.frame(kd = character(0), cat = character(0), medikd = numeric(0), mediu6 = numeric(0), p.value = numeric(0))
result_KW <- data.frame(kd = character(0), p.valuekd = numeric(0), p.valueu6 = numeric(0))
result_DT <- data.frame(kd = character(0), kd = character(0), group1 = character(0), group2 = character(0), n1 = numeric(0), n1 = numeric(0), statistic = numeric(0), p = numeric(0), p.adj = numeric(0), p.adj.signif = character(0))

#FOR AUTOSOMAL  
for (j in 1:6){
  
  kd_dip11s <- dflist[[j]]
  kd <- id[[j]]
  kd2 <- id2[[j]]

  kd_dip11_ex <- kd_dip11s %>% mutate(cat = (case_when(map >= 0.50 & map < 0.65 ~ "balanced", 
                                                      (map >= 0.65 & map <= 0.85) ~ "inter",
                                                      map > 0.85 ~ "skewed"))) %>% mutate(log2qr = log2(qr)) %>% ungroup()
  
  kd_dip11_ex$log2qr %>% as.numeric(kd_dip11_ex$log2qr)
  
  #FOR individual SNVs - exons
specA <- kd_dip11_ex %>% filter(cat != "NA") %>% filter(XA == "A") %>% filter(exon_start != "NA") %>% filter(log2qr != "-Inf") %>% distinct(chromosome, pos, id, id2, .keep_all = TRUE)
  
  for (k in category){   
    data <- specA %>% filter(id != "U6M2") %>% dplyr::filter(cat == k)
    WX <- rstatix::pairwise_wilcox_test(data, log2qr ~ id2, paired = FALSE)
    datakd <- data %>% filter(id2 != "U6M2_qr")
    datau6 <- data %>% filter(id2 == "U6M2_qr")
    medi_kd <- median(datakd$log2qr)
    medi_u6 <- median(datau6$log2qr)
    res <- c(kd, k, medi_kd, medi_u6, WX$p)
    result_WX[nrow(result_WX)+1,] <- res
  }
  
  exp_plots[[j]] = ggplot(specA %>% filter(id != "U6M2"), aes(x = cat, y = log2qr))+
    geom_boxplot(aes(fill = id2), notch = TRUE, outliers = FALSE, alpha = 0.8)+
    theme_classic(base_size = 16)+
    ylab("")+
    xlab("")+
    labs(subtitle = paste0(kd))+
    geom_text(data = specA %>% filter(cat == "balanced") %>% filter(id != "U6M2") %>% filter(id2 == kd2), stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=1, size=3)+
    geom_text(data = specA %>% filter(cat == "balanced") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), stat = "count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=-0.5, size=3)+
    geom_text(data = specA %>% filter(cat == "inter") %>% filter(id != "U6M2") %>% filter(id2 == kd2), stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=1, size=3)+
    geom_text(data = specA %>% filter(cat == "inter") %>% filter(id != "U6M2")  %>% filter(id2 == "U6M2_qr"), stat = "count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=-0.5, size=3)+
    geom_text(data = specA %>% filter(cat == "skewed") %>% filter(id != "U6M2") %>% filter(id2 == kd2), stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=1, size=3)+
    geom_text(data = specA %>% filter(cat == "skewed") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), stat = "count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=-0.5, size=3)+
    stat_summary(data = specA %>% filter(cat == "balanced") %>% filter(id != "U6M2") %>% filter(id2 == kd2), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = 1, vjust = -0.5))+
    stat_summary(data = specA %>% filter(cat == "balanced") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = -0.5, vjust = -0.5))+
    stat_summary(data = specA %>% filter(cat == "inter") %>% filter(id != "U6M2") %>% filter(id2 == kd2), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = 1, vjust = -0.5))+
    stat_summary(data = specA %>% filter(cat == "inter") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = -0.5, vjust = -0.5))+
    stat_summary(data = specA %>% filter(cat == "skewed") %>% filter(id != "U6M2") %>% filter(id2 == kd2), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = 1, vjust = -0.5))+
    stat_summary(data = specA %>% filter(cat == "skewed") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = -0.5, vjust = -0.5))+
    theme(legend.position="none")
  print(exp_plots[[j]])
  
  exp_binned_histo[[j]] = ggplot()+
    geom_density(data = specA %>% filter(id != "U6M2"), aes(x = map, y = ..density..), fill = "#F8766D", alpha = 0.5)+
    geom_density(data = specA %>% filter(id == "U6M2"), aes(x = map, y = ..density..), fill = "#00BFC4", alpha = 0.5)+
    theme_classic(base_size = 16)+
    scale_x_continuous(limits = c(0.5, 1))+
    xlab("")+
    ylab("")+ #gene number
    theme(legend.position = "none")+
    labs(subtitle = paste0(kd))
  print(exp_binned_histo[[j]])
  
  KWkd <- rstatix::kruskal_test(data = specA %>% filter(id != "U6M2") %>% filter(id2 != "U6M2_qr"), log2qr ~ cat)
  KWu6 <- rstatix::kruskal_test(data = specA %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), log2qr ~ cat)
  KWkd$p <- format(KWkd$p, scientific = TRUE, digits = 3)
  KWu6$p <- format(KWu6$p, scientific = TRUE, digits = 3)
  res <- c(kd, KWkd$p, KWu6$p)
  result_KW[nrow(result_KW)+1,] <- res
  
  DTkd <- rstatix::dunn_test(data = specA %>% filter(id != "U6M2") %>% filter(id2 != "U6M2_qr"), log2qr ~ cat, p.adjust.method = "holm")
  DTu6 <- rstatix::dunn_test(data = specA %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), log2qr ~ cat, p.adjust.method = "holm")
  res1 <- c(kd, kd2, DTkd[1,2:9])
  res2 <- c(kd, kd2, DTkd[2,2:9])
  res3 <- c(kd, kd2, DTkd[3,2:9])
  res4 <- c(kd, "U6M2_qr", DTu6[1,2:9])
  res5 <- c(kd, "U6M2_qr", DTu6[2,2:9])
  res6 <- c(kd, "U6M2_qr", DTu6[3,2:9])
  result_DT[nrow(result_DT)+1,] <- res1
  result_DT[nrow(result_DT)+1,] <- res2
  result_DT[nrow(result_DT)+1,] <- res3
  result_DT[nrow(result_DT)+1,] <- res4
  result_DT[nrow(result_DT)+1,] <- res5
  result_DT[nrow(result_DT)+1,] <- res6
}

p1 <- as.ggplot(exp_plots[[1]])
p2 <- as.ggplot(exp_plots[[2]])
p3 <- as.ggplot(exp_plots[[3]])
p4 <- as.ggplot(exp_plots[[4]])
p5 <- as.ggplot(exp_plots[[5]])
p6 <- as.ggplot(exp_plots[[6]])

plots <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3) #+ plot_annotation(title = "Autosomal specific SNVs \nGene expression by allelic balance")
plots
ggsave(file = "exp_binned_plots_exon_A_SNV_05.svg", plot = plots,  device = "svg", width = 15, height = 12, units = "in")

p7 <- as.ggplot(exp_binned_histo[[1]])
p8 <- as.ggplot(exp_binned_histo[[2]])
p9 <- as.ggplot(exp_binned_histo[[3]])
p10 <- as.ggplot(exp_binned_histo[[4]])
p11 <- as.ggplot(exp_binned_histo[[5]])
p12 <- as.ggplot(exp_binned_histo[[6]])

plots <- p7 + p8 + p9 + p10 + p11 + p12 + plot_layout(ncol = 3)# + plot_annotation(title = "Autosomal specific SNVs \nGene expression by allelic balance")
plots
ggsave(file = "exp_binned_histos_exon_A_SNV_05.svg", plot = plots,  device = "svg", width = 13, height = 7, units = "in")

result_WX <- as_tibble(result_WX)
result_WX$p.value <- as.numeric(result_WX$p.value)
result_WX$medikd <- as.numeric(result_WX$medikd)
result_WX$mediu6 <- as.numeric(result_WX$mediu6)
result_WX <- result_WX %>% mutate("p.signif" = case_when(p.value <= 0.0001 ~ "****", p.value <= 0.001 ~ "***", p.value <= 0.01 ~ "**", p.value <= 0.05 ~ "*", p.value > 0.05 ~ "NS"))
result_WX$p.value <- format(result_WX$p.value, scientific = TRUE, digits = 3)
result_WX$medikd <- format(result_WX$medikd, scientific = TRUE, digits = 5)
result_WX$mediu6 <- format(result_WX$mediu6, scientific = TRUE, digits = 5)
write.table(result_WX, file="Wilcox_test_exp_binned_exon_A_SNV_05.txt", quote=FALSE, sep="\t", row.names = FALSE)

result_KW <- as_tibble(result_KW)
result_KW$p.valuekd <- as.numeric(result_KW$p.valuekd)
result_KW$p.valueu6 <- as.numeric(result_KW$p.valueu6)
result_KW$p.valuekd <- format(result_KW$p.valuekd, scientific = TRUE, digits = 3)
result_KW$p.valueu6 <- format(result_KW$p.valueu6, scientific = TRUE, digits = 3)
write.table(result_KW, file="KW_test_exp_binned_exon_A_SNV_05.txt", quote=FALSE, sep="\t", row.names = FALSE)

result_DT <- as_tibble(result_DT)
write.table(result_DT, file="DT_test_exp_binned_exon_A_SNV_05.txt", quote=FALSE, sep="\t", row.names = FALSE)


#FOR X chromosome
dflist <- list(kd_dip11_caprin, kd_dip11_ckap, kd_dip11_hnrnpk, kd_dip11_ncl, kd_dip11_nono, kd_dip11_syncrip) #for SNVs
id <- list("CAPRIN1", "CKAP4", "HNRNPK", "NCL", "NONO", "SYNCRIP")
id2 <- list("CAPRIN1_qr", "CKAP4_qr", "HNRNPK_qr", "NCL_qr", "NONO_qr", "SYNCRIP_qr")
exp_plots = list()
exp_binned_histo = list()
category <- list("balanced", "skewed")
result_WX <- data.frame(kd = character(0), cat = character(0), medikd = numeric(0), mediu6 = numeric(0), p.value = numeric(0))
result_KW <- data.frame(kd = character(0), p.valuekd = numeric(0), p.valueu6 = numeric(0))
result_DT <- data.frame(kd = character(0), kd = character(0), group1 = character(0), group2 = character(0), n1 = numeric(0), n1 = numeric(0), statistic = numeric(0), p = numeric(0), p.adj = numeric(0), p.adj.signif = character(0))

for (j in 1:6){
  
  kd_dip11s <- dflist[[j]]
  kd <- id[[j]]
  kd2 <- id2[[j]]
  
  kd_dip11_ex <- kd_dip11s %>% mutate(cat = (case_when(map >= 0.50 & map <= 0.80 ~ "balanced", map > 0.80 ~ "skewed"))) %>% mutate(log2qr = log2(qr)) %>% ungroup()
  
  kd_dip11_ex$log2qr %>% as.numeric(kd_dip11_ex$log2qr)
  
  #FOR individual SNVs - exons
 specA <- kd_dip11_ex %>% filter(cat != "NA") %>% filter(XA == "X") %>% filter(exon_start != "NA") %>% filter(log2qr != "-Inf") %>% distinct(chromosome, pos, id, id2, .keep_all = TRUE)

  for (k in category){   
    data <- specA %>% filter(id != "U6M2") %>% dplyr::filter(cat == k)
    WX <- rstatix::pairwise_wilcox_test(data, log2qr ~ id2, paired = FALSE)
    datakd <- data %>% filter(id2 != "U6M2_qr")
    datau6 <- data %>% filter(id2 == "U6M2_qr")
    medi_kd <- median(datakd$log2qr)
    medi_u6 <- median(datau6$log2qr)
    res <- c(kd, k, medi_kd, medi_u6, WX$p)
    result_WX[nrow(result_WX)+1,] <- res
  }
  
  exp_plots[[j]] = ggplot(specA %>% filter(id != "U6M2"), aes(x = cat, y = log2qr))+
    geom_boxplot(aes(fill = id2), notch = TRUE, outliers = FALSE, alpha = 0.8)+
    theme_classic(base_size = 16)+
    ylab("")+
    xlab("")+
    labs(subtitle = paste0(kd))+
    geom_text(data = specA %>% filter(cat == "balanced") %>% filter(id != "U6M2") %>% filter(id2 == kd2), stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=1, size=3)+
    geom_text(data = specA %>% filter(cat == "balanced") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), stat = "count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=-0.5, size=3)+
    geom_text(data = specA %>% filter(cat == "skewed") %>% filter(id != "U6M2") %>% filter(id2 == kd2), stat="count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=1, size=3)+
    geom_text(data = specA %>% filter(cat == "skewed") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), stat = "count", aes(label=..count..), y=-Inf, vjust=-0.5, hjust=-0.5, size=3)+
    stat_summary(data = specA %>% filter(cat == "balanced") %>% filter(id != "U6M2") %>% filter(id2 == kd2), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = 1, vjust = -0.5))+
    stat_summary(data = specA %>% filter(cat == "balanced") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = -0.5, vjust = -0.5))+
    stat_summary(data = specA %>% filter(cat == "skewed") %>% filter(id != "U6M2") %>% filter(id2 == kd2), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = 1, vjust = -0.5))+
    stat_summary(data = specA %>% filter(cat == "skewed") %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = -0.5, vjust = -0.5))+
    theme(legend.position="none")
  print(exp_plots[[j]])
  
  exp_binned_histo[[j]] = ggplot()+
    geom_density(data = specA %>% filter(id != "U6M2"), aes(x = map, y = ..density..), fill = "#F8766D", alpha = 0.5)+
    geom_density(data = specA %>% filter(id == "U6M2"), aes(x = map, y = ..density..), fill = "#00BFC4", alpha = 0.5)+
    theme_classic(base_size = 16)+
    scale_x_continuous(limits = c(0.5, 1))+
    xlab("")+
    ylab("")+ #gene number
    theme(legend.position = "none")+
    labs(subtitle = paste0(kd))
  print(exp_binned_histo[[j]])
  
  KWkd <- rstatix::kruskal_test(data = specA %>% filter(id != "U6M2") %>% filter(id2 != "U6M2_qr"), log2qr ~ cat)
  KWu6 <- rstatix::kruskal_test(data = specA %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), log2qr ~ cat)
  KWkd$p <- format(KWkd$p, scientific = TRUE, digits = 3)
  KWu6$p <- format(KWu6$p, scientific = TRUE, digits = 3)
  res <- c(kd, KWkd$p, KWu6$p)
  result_KW[nrow(result_KW)+1,] <- res
  
  DTkd <- rstatix::dunn_test(data = specA %>% filter(id != "U6M2") %>% filter(id2 != "U6M2_qr"), log2qr ~ cat, p.adjust.method = "holm")
  DTu6 <- rstatix::dunn_test(data = specA %>% filter(id != "U6M2") %>% filter(id2 == "U6M2_qr"), log2qr ~ cat, p.adjust.method = "holm")
  res1 <- c(kd, kd2, DTkd[1,2:9])
  res2 <- c(kd, kd2, DTkd[2,2:9])
  res3 <- c(kd, "U6M2_qr", DTu6[1,2:9])
  res4 <- c(kd, "U6M2_qr", DTu6[2,2:9])
  result_DT[nrow(result_DT)+1,] <- res1
  result_DT[nrow(result_DT)+1,] <- res2
  result_DT[nrow(result_DT)+1,] <- res3
  result_DT[nrow(result_DT)+1,] <- res4
}

p1 <- as.ggplot(exp_plots[[1]])
p2 <- as.ggplot(exp_plots[[2]])
p3 <- as.ggplot(exp_plots[[3]])
p4 <- as.ggplot(exp_plots[[4]])
p5 <- as.ggplot(exp_plots[[5]])
p6 <- as.ggplot(exp_plots[[6]])

plots <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3)
plots
ggsave(file = "exp_binned_plots_exon_X_SNV_05.svg", plot = plots,  device = "svg", width = 10, height = 8, units = "in")

p7 <- as.ggplot(exp_binned_histo[[1]])
p8 <- as.ggplot(exp_binned_histo[[2]])
p9 <- as.ggplot(exp_binned_histo[[3]])
p10 <- as.ggplot(exp_binned_histo[[4]])
p11 <- as.ggplot(exp_binned_histo[[5]])
p12 <- as.ggplot(exp_binned_histo[[6]])

plots <- p7 + p8 + p9 + p10 + p11 + p12 + plot_layout(ncol = 3)# + plot_annotation(title = "Autosomal specific SNVs \nGene expression by allelic balance")
plots
ggsave(file = "exp_binned_histos_exon_X_SNV_05.svg", plot = plots,  device = "svg", width = 10, height = 7, units = "in")

result_WX <- as_tibble(result_WX)
result_WX$p.value <- as.numeric(result_WX$p.value)
result_WX$medikd <- as.numeric(result_WX$medikd)
result_WX$mediu6 <- as.numeric(result_WX$mediu6)
result_WX <- result_WX %>% mutate("p.signif" = case_when(p.value <= 0.0001 ~ "****", p.value <= 0.001 ~ "***", p.value <= 0.01 ~ "**", p.value <= 0.05 ~ "*", p.value > 0.05 ~ "NS"))
result_WX$p.value <- format(result_WX$p.value, scientific = TRUE, digits = 3)
result_WX$medikd <- format(result_WX$medikd, scientific = TRUE, digits = 5)
result_WX$mediu6 <- format(result_WX$mediu6, scientific = TRUE, digits = 5)
write.table(result_WX, file="Wilcox_test_exp_binned_exon_X_SNV_05.txt", quote=FALSE, sep="\t", row.names = FALSE)

result_KW <- as_tibble(result_KW)
result_KW$p.valuekd <- as.numeric(result_KW$p.valuekd)
result_KW$p.valueu6 <- as.numeric(result_KW$p.valueu6)
result_KW$p.valuekd <- format(result_KW$p.valuekd, scientific = TRUE, digits = 3)
result_KW$p.valueu6 <- format(result_KW$p.valueu6, scientific = TRUE, digits = 3)
write.table(result_KW, file="KW_test_exp_binned_exon_X_SNV_05.txt", quote=FALSE, sep="\t", row.names = FALSE)

result_DT <- as_tibble(result_DT)
write.table(result_DT, file="DT_test_exp_binned_exon_X_SNV_05.txt", quote=FALSE, sep="\t", row.names = FALSE)


#FIGURE 3-12 
#AUTOSOMAL UPREGULATED RNA-SEQ KNOCKDOWNS - SKEW CHANGE 
setwd("/Users/.../MDO_diploid_knockdowns/MDO_diploid_kd_normalised_counts_fc_mMonDom1/")

up_caprin <- read.table("up_A_caprin_genes.txt",header=TRUE, sep="\t")
up_ckap <- read.table("up_A_ckap_genes.txt",header=TRUE, sep="\t")
up_hnrnpk <- read.table("up_A_hnrnpk_genes.txt",header=TRUE, sep="\t")
up_syncrip <- read.table("up_A_syncrip_genes.txt",header=TRUE, sep="\t")

kd_dip5a <- kd_dip5 %>% filter(exon_start != "NA")
up_caprina <- left_join(up_caprin, kd_dip5a)
count_caprin <- up_caprina %>% dplyr::count(chromosome)

up_ckapa <- left_join(up_ckap, kd_dip5a)
count_ckap <- up_ckapa %>% dplyr::count(chromosome)

up_hnrnpka <- left_join(up_hnrnpk, kd_dip5a)
count_hnrnpk <- up_hnrnpka %>% dplyr::count(chromosome)

up_syncripa <- left_join(up_syncrip, kd_dip5a)
count_syncrip <- up_syncripa %>% dplyr::count(chromosome)

#INDIVIDUAL SNVS
kd_dip9_caprin <- up_caprina %>% filter(chromosome != "NA") %>% filter(CAPRIN1 >=0.5 | U6M2 >=0.5 ) %>% filter(u6m2_cpm>0.5 | caprin1_cpm>0.5) %>% dplyr::select(-CKAP4, -HNRNPK, -NCL, -NONO, -SYNCRIP, -NO_VECTOR) %>% pivot_longer(cols = c(CAPRIN1, U6M2), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE)
kd_dip9_ckap <- up_ckapa %>% filter(chromosome != "NA") %>% filter(CKAP4 >=0.5 | U6M2 >=0.5 ) %>% filter(u6m2_cpm>0.5 | ckap4_cpm>0.5) %>% dplyr::select(-CAPRIN1, -HNRNPK, -NCL, -NONO, -SYNCRIP, -NO_VECTOR) %>% pivot_longer(cols = c(CKAP4, U6M2), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE)
kd_dip9_hnrnpk <- up_hnrnpka %>% filter(chromosome != "NA") %>% filter(HNRNPK >=0.5 | U6M2 >=0.5 ) %>% filter(u6m2_cpm>0.5 | hnrnpk_cpm>0.5) %>% dplyr::select(-CKAP4, -CAPRIN1, -NCL, -NONO, -SYNCRIP, -NO_VECTOR) %>% pivot_longer(cols = c(HNRNPK, U6M2), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE)
kd_dip9_syncrip <- up_syncripa %>% filter(chromosome != "NA") %>% filter(SYNCRIP >=0.5 | U6M2 >=0.5 ) %>% filter(u6m2_cpm>0.5 | syncrip_cpm>0.5) %>% dplyr::select(-CKAP4, -HNRNPK, -NCL, -NONO, -CAPRIN1, -NO_VECTOR) %>% pivot_longer(cols = c(SYNCRIP, U6M2), names_to = "id", values_to = "map") %>% distinct(chromosome, pos, id, .keep_all = TRUE)

dflist <- list(kd_dip9_caprin, kd_dip9_ckap, kd_dip9_hnrnpk, kd_dip9_nono, kd_dip9_syncrip)

num_files = length(dflist)
target <- list("CAPRIN1", "CKAP4", "HNRNPK", "SYNCRIP")
kd_map_boxplots = list()
kd_map_overlay_histo_A = list()
kd_map_overlay_histo_X = list()
kd_map_counts <- data.frame(kd = character(0), Akdc = numeric(0), Au6c = numeric(0), moodAp = numeric(0))
results_WX <- data.frame(kd = character(0), p.value_a = numeric(0))

for (j in 1:5){
  
  kd_dip9 <- dflist[[j]]
  
  kd <- target[[j]]
  
  setA <- kd_dip9 
  
  #counts and medians  
  Akd <- setA %>% filter(id == kd) 
  Au6 <- setA %>% filter(id == "U6M2") 
  Akdc <- Akd %>% dplyr::count(id) 
  Au6c <- Au6 %>% dplyr::count(id) 
  moodA <- mood.medtest(map ~ id, data = setA)
  moodAp <- moodA$p.value
  res1 <- c(kd, Akdc$n, Au6c$n, moodAp)
  kd_map_counts[nrow(kd_map_counts)+1,] <- res1
  
  #overlay histograms of distribution of map in kd v control  
  kd_map_overlay_histo_A[[j]] <- ggplot(setA, aes(x = map, color = id))+
    geom_histogram(bins = 20, alpha = 0.10, position = "identity", aes(fill = id))+
    theme_classic(base_size = 16)+
    labs(subtitle = paste0(target[[j]]))+
    scale_y_continuous(limits = c(0, 180))+
    scale_x_continuous(limits = c(0, 1))+
    scale_color_manual(values = c("dodgerblue3", "orange"))+
    scale_fill_manual(values = c("dodgerblue3", "orange"))+
    xlab("")+
    ylab("")+
    theme(legend.title = element_blank())+
    theme(legend.position = "none")
  print(kd_map_overlay_histo_A[[j]])
  
  #boxplots
  kd_map_boxplots[[j]] <- ggplot(kd_dip9, aes(x = id, y = map, fill = id))+
    geom_boxplot(notch = TRUE, outliers = FALSE)+ 
    theme_classic(base_size = 20)+
    ylab("")+
    xlab("")+
    scale_y_continuous(breaks = c(0.4,0.5, 0.6, 0.7, 0.8, 0.9, 1.0))+
    coord_cartesian(ylim = c(0.2, 1.2))+
    geom_text(data = Akd %>% dplyr::count(id), aes(y = 0.25, label = n), hjust = 1, size = 5)+
    geom_text(data = Au6 %>% dplyr::count(id), aes(y = 0.25, label = n), hjust = 1, size = 5)+
    stat_summary(data = Akd, fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = 1, vjust = -1))+
    stat_summary(data = Au6, fun = "median", geom = "text", aes(label = round(after_stat(y),3), hjust = 1, vjust = -1))+
    theme(legend.title = element_blank())+
    theme(legend.position = "none")
  print(kd_map_boxplots[[j]])
  
  #wilcoxon two-sample tests for X v A and for kd v control  
  WXa <- rstatix::wilcox_test(data = setA, map ~ id, paired = TRUE)
  res2 <- c(kd, WXa$p)
  results_WX[nrow(results_WX)+1,] <- res2
  
}  

kd_map_counts <- as_tibble(kd_map_counts)
kd_map_counts$moodAp <- as.numeric(kd_map_counts$moodAp)
kd_map_counts$moodAp <- format(kd_map_counts$moodAp, scientific =  TRUE, digits = 3)
write.table(kd_map_counts, file="A_up_SNV_counts_all_05.txt", quote=FALSE, sep="\t", row.names = FALSE)

results_WX <- as_tibble(results_WX)
results_WX$p.value_a <- as.numeric(results_WX$p.value_a)
write.table(results_WX, file="A_up_gene_med_WX_05.txt", quote=FALSE, sep="\t", row.names = FALSE)

p1 <- as.ggplot(kd_map_boxplots[[1]])
p2 <- as.ggplot(kd_map_boxplots[[2]])
p3 <- as.ggplot(kd_map_boxplots[[3]])
p4 <- as.ggplot(kd_map_boxplots[[4]])

# use patchwork to arrange them together
#select only HNRNPK and CAPRIN1 for Figure 3-12
plots <- p1 + p3 + plot_layout(ncol = 2) 
plots
ggsave(file = "A_up_SNV_counts_boxplots_05.svg", plot = plots,  device = "svg", width = 6, height = 5, units = "in")

p7 <- as.ggplot(kd_map_overlay_histo_A[[1]])
p8 <- as.ggplot(kd_map_overlay_histo_A[[2]])
p9 <- as.ggplot(kd_map_overlay_histo_A[[3]])
p10 <- as.ggplot(kd_map_overlay_histo_A[[4]])

#select only HNRNPK and CAPRIN1 for Figure 3-12
plots2 <- p7 + p9 + plot_layout(ncol = 2) 
plots2
ggsave(file = "A_up_SNV_counts_histo_05.svg", plot = plots2,  device = "svg", width = 6, height = 4, units = "in")


