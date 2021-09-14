library(readr)
library(tidyverse)
library(ggplot2)
library(splitstackshape)
library(ggpubr)


source_data = "/Users/danielagaio/Gaio/github/metapigs_base/source_data/" # git 
middle_dir = "/Users/danielagaio/Gaio/github/metapigs_base/middle_dir/" # git 
stats_dir = "/Users/danielagaio/Gaio/github/metapigs_base/middle_dir/stats/" # git 
out_dir_git = "/Users/danielagaio/Gaio/github/metapigs_base/out/" # git 
out_dir = "/Users/danielagaio/Desktop/metapigs_base/phylosift/out/" # local 


readcounts_samples <- read.csv("~/Gaio/github/metapigs_base/source_data/readcounts_samples.csv", header=FALSE)
NROW(readcounts_samples)

r <- cSplit(readcounts_samples, "V1","_")
r$DNA_plate <- paste0(r$V1_1,"_",r$V1_2)
r <- cSplit(r, "V1_3",".")

r <- r %>%
  dplyr::select(DNA_plate,V1_3_1, V3)
colnames(r) <- c("DNA_plate","DNA_well","read_count")
head(r)

PD_1M <- read_csv("Gaio/github/metapigs_base/middle_dir/fpddat_clean")
PD_100K <- read_csv("Gaio/github/metapigs_base/middle_dir/old/fpddat_clean")

NROW(PD_1M)
NROW(PD_100K)

both <- inner_join(PD_100K,PD_1M, by = c("DNA_plate", "DNA_well")) %>%
  dplyr::select(DNA_plate,DNA_well,unrooted_pd.x,unrooted_pd.y,bwpd.x, bwpd.y)
colnames(both) <- gsub("x","PD_100k",colnames(both))
colnames(both) <- gsub("y","PD_1M",colnames(both))
both$DNA_plate <- gsub("P","plate_",both$DNA_plate)


both %>%
  dplyr::summarise(unroo_PD_100K_median=median(unrooted_pd.PD_100k),
                   unroo_PD_100K_sd=sd(unrooted_pd.PD_100k),
                unroo_PD_1M_median=median(unrooted_pd.PD_1M),
                unroo_PD_1M_sd=sd(unrooted_pd.PD_1M),
                bwpd_PD_100K_median=median(bwpd.PD_100k),
                bwpd_PD_100K_sd=sd(bwpd.PD_100k),
                bwpd_PD_1M_median=median(bwpd.PD_1M),
                bwpd_PD_1M_sd=sd(bwpd.PD_1M))
  
cor(x = both$unrooted_pd.PD_100k,y=both$unrooted_pd.PD_1M, method = "pearson")
cor(x = both$bwpd.PD_100k,y=both$bwpd.PD_1M, method = "pearson")

# merge read count info to PD info: 
df <- inner_join(both,r, by=c("DNA_plate","DNA_well"))
df <- na.omit(df)
NROW(df)

# median read count per sample, percentage corresponding to analysis (100k vs 1M)
df %>%
  dplyr::summarise(median=median(read_count),
                   sd=sd(read_count)) %>%
  dplyr::mutate(PD_100K=(100*100000)/median,
                PD_1M=(100*1000000)/median)







