
# run on HPC
# R: /shared/homes/12705859/miniconda3/envs/R36/bin/R

library(readr)
library(data.table)
library(splitstackshape)
library(readxl)
library(dplyr)
library(tidyverse)



HPC_dir = "/shared/homes/12705859/sortmerna_aligned/all_aligned_fq/filtered_fq/classified_out/" # HPC 


######################################################################

# load RDP classifications (file is 7.73Gb so this is run directly on HPC)
classifications_all <- read.delim(paste0(HPC_dir, "classifications_all.tsv"), header=FALSE)
head(classifications_all)

######################################################################


# load metadata 
mdat <- read_excel(paste0(HPC_dir,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)


######################################################################


# parse metadata 

mdat$`*collection_date` <- as.character(mdat$`*collection_date`)
mdat$Cohort <- gsub("Sows","Sows",mdat$Cohort)
mdat$Cohort <- gsub("D-scour","D-Scour", mdat$Cohort)

colnames(mdat)[colnames(mdat) == '*collection_date'] <- 'collection_date'

mdat <- mdat %>%
  dplyr::select(DNA_plate,DNA_well,Cohort,collection_date,isolation_source)

mdat$DNA_plate_well <- paste0(mdat$DNA_plate,"_",mdat$DNA_well)
head(mdat)


######################################################################

classifications_all$V1<-gsub("_S.*","",classifications_all$V1)

colnames(classifications_all)[colnames(classifications_all) == 'V1'] <- 'DNA_plate_well'

df <- merge(classifications_all,mdat)

head(df)


# function to get numbers: 
get_me_stats <- function(DF) {
  out <- DF %>% 
    dplyr::summarise(n=n()) %>%
    dplyr::mutate(perc=n/sum(n)*100) %>%
    dplyr::arrange(desc(perc))
  return(out)
}

# assign "pig" vs "mom" tag: 
df <- df %>%
  dplyr::mutate(is_mother=Cohort=="Mothers") 
head(df)



sink(file = paste0(HPC_dir,"RDP_analysis.txt"), 
     append = FALSE, type = c("output"))
paste0("########### RDP analysis ###########")
paste0("##################################")
paste0("Confidence of phylum classification (scale 0-1)")
summary(df$V11) 
paste0("Distribution by domain (remember sortmerna reference used was bacterial 16S)")
df %>% 
  group_by(V6) %>% #phylum
  get_me_stats()
paste0("##################################")
paste0("Most common phyla (%)")
df %>% 
  group_by(V9) %>% #phylum
  get_me_stats()
paste0("##################################")
paste0("Most common classes (%)")
df %>% 
  group_by(V12) %>% #class
  get_me_stats()
paste0("##################################")
paste0("Most common families (%)")
df %>% 
  group_by(V18) %>% #family
  get_me_stats()
paste0("##################################")
paste0("Most common genera (%)")
df %>% 
  group_by(V21) %>% #genus
  get_me_stats()
paste0("##################################")
paste0("Most common phyla in mothers (%)")
df %>% 
  dplyr::filter(is_mother==TRUE) %>%
  group_by(V9) %>% #phylum
  get_me_stats()
paste0("##################################")
paste0("Most common phyla in piglets (%)")
df %>% 
  dplyr::filter(is_mother==FALSE) %>%
  group_by(V9) %>% #phylum
  get_me_stats()
paste0("##################################")
sink()


