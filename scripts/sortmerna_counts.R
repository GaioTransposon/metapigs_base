
library(readr)
library(data.table)
library(splitstackshape)
library(readxl)
library(dplyr)
library(tidyverse)



source_data = "/Users/12705859/metapigs_base/source_data/" # git 
middle_dir = "/Users/12705859/metapigs_base/middle_dir/" # git 
out_dir = "/Users/12705859/metapigs_base/out/" # git
out_dir_local = "/Users/12705859/Desktop/metapigs_base/reads_counts_distribution/" # local 


######################################################################


# load read counts file 

readcounts <- read_csv(paste0(middle_dir,"readcounts.csv"), 
                               col_names = TRUE)
head(readcounts)


######################################################################


# load log of reads extracted with sortmerna:

sortmerna_log <- read_csv(paste0(middle_dir,"sortmerna_log"), # from HPC: /shared/homes/12705859/sortmerna_aligned
                          col_names = FALSE)
head(sortmerna_log)


######################################################################


sink(file = paste0(out_dir,"sortmerna_stats.txt"), 
     append = FALSE, type = c("output"))
paste0("########### SORTMERNA STATS ###########")
sink()

# parse sortmerna log 
sortmerna_log$X2 <- rep(c("file","passed","failed"),1,nrow(sortmerna_log))

sortmerna_log <- as.data.frame(sortmerna_log)

grouped <- sortmerna_log %>% 
  group_by(group = as.integer(gl(n(), 3, n())))

grouped2 <- grouped %>%
  pivot_wider(names_from = X2, values_from=X1)

grouped2 <- cSplit(grouped2, "passed", " ")
grouped2 <- cSplit(grouped2, "failed", " ")
grouped2$file <- gsub("Reads file: /shared/homes/s1/pig_microbiome/sortmerna_16S/","",grouped2$file)
grouped2$file <- gsub("(_S).*","",grouped2$file)

grouped2$passed_8 <- as.numeric(gsub("[^0-9.]","",grouped2$passed_8))
grouped2$failed_8 <- as.numeric(gsub("[^0-9.]","",grouped2$failed_8))

grouped2 <- grouped2 %>%
  dplyr::select(file, passed_7)

colnames(grouped2) <- c("DNA_plate_well","passed_count")

sortme <- grouped2

merged <- merge(readcounts,sortme) %>%
  dplyr::mutate(perc_16S = (passed_count/counts)*100) 

head(merged)
NROW(merged)

sink(file = paste0(out_dir,"sortmerna_stats.txt"), 
     append = TRUE, type = c("output"))
paste0("##################################")
paste0("sum of reads that passed sortmerna filtering (e-value<=0.0001)) : ",
       sum(merged$passed_count))
paste0("##################################")
paste0("16S reads in samples : ")
summary(merged$perc_16S)
paste0("##################################")
paste0("16S reads in samples (% of tot. reads) : ")
summary(merged$perc_16S)
paste0("##################################")
paste0("16S reads in samples - 95 percent confidence interval : ")
t.test(merged$perc_16S)
paste0("##################################")
sink()


