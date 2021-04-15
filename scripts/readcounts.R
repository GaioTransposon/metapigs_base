
library(readr)
library(splitstackshape)
library(readxl)
library(dplyr)
library(tidyverse)
library(data.table)


source_data = "/Users/12705859/metapigs_base/source_data/" # git 
middle_dir = "/Users/12705859/metapigs_base/middle_dir/" # git 
out_dir = "/Users/12705859/metapigs_base/out/" # git
out_dir_local = "/Users/12705859/Desktop/metapigs_base/reads_counts_distribution/" # local 


sink(file = paste0(out_dir,"read_counts_stats.txt"), 
     append = FALSE, type = c("output"))
paste0("########### READ COUNT STATS ###########")
sink()

###########################################################################################

# raw read counts (only R1)
readcounts_samples <- read_csv(paste0(source_data,"readcounts_samples.csv"),
                               col_names = FALSE)
head(readcounts_samples)


######################################################################


# clean paired reads 
lib <- read_delim(paste0(middle_dir,"clean_paired_lib_sizes_final.tsv"),
                  "\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)

sink(file = paste0(out_dir,"read_counts_stats.txt"), 
     append = TRUE, type = c("output"))
paste0("##################################") 
paste0("Summary of read counts from clean paired libraries ",
       summary(as.data.frame(lib$X2)))
paste0("##################################")
paste0("Sum of read counts from clean paired libraries : ",
       sum(as.data.frame(lib$X2)))
paste0("##################################")
paste0("number of samples included (only piglet samples) : ", NROW(lib))     
paste0("##################################") 
sink()

######################################################################


# load metadata 
mdat <- read_excel(paste0(source_data,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)


######################################################################

# split first column to contain plate and well
reads <- cSplit(readcounts_samples, "X1", ".")

reads$X2 <- NULL
reads$X1_2 <- NULL
reads$X1_3 <- NULL
reads$X1_4 <- NULL

# stats read counts: 
summary(reads$X3*2)
sum(reads$X3*2)

lowreads_samples <- hist(head(reads$X3)[
  readcounts_samples$X3<=40000], 
  breaks=20, 
  xlim=c(0,40000),
  xlab = "read counts",
  ylab = "Frequency",
  main= NULL)

histogram <- hist(readcounts_samples$X3, breaks = 200,
                  main = "Read count distribution across samples",
                  xlab = "read counts",
                  ylab = "Frequency (sample counts)")

pdf(paste0(out_dir_local, "readcounts_distribution.pdf"), onefile = TRUE)
histogram
dev.off()

pdf(paste0(out_dir_local, "readcounts_distribution_all_and_low.pdf"),onefile = TRUE)
par(mfrow=c(1,1))
hist(readcounts_samples$X3[
  readcounts_samples$X3<=165291790], 
     breaks = seq(from=1, to=165291790, by=1000000), #165291790/1 000 000 -> 200 bins! 
     cex.axis=1,cex.lab=1.5, cex.main=2,
     main = "Read count distribution across samples",
     xlab = "read counts",
     ylab = "Frequency (N samples)")
par(fig=c(0.45, 0.99, 0.2, 1), new = T) 
hist(readcounts_samples$X3[
  readcounts_samples$X3<=30000], 
  breaks = seq(from=0, to=30000, by=1500), #30000/1500 -> 20 bins
  xlim=c(0,30000),
  xlab = NULL,
  ylab = NULL,
  main= NULL)
dev.off()



######################################################################

# parse metadata to merge to read count info 

mdat$`*collection_date` <- as.character(mdat$`*collection_date`)
mdat$Cohort <- gsub("Sows","Sows",mdat$Cohort)
mdat$Cohort <- gsub("D-scour","D-Scour", mdat$Cohort)

colnames(mdat)[colnames(mdat) == '*collection_date'] <- 'collection_date'

mdat <- mdat %>%
  dplyr::select(DNA_plate,DNA_well,Cohort,collection_date,isolation_source)

mdat$DNA_plate_well <- paste0(mdat$DNA_plate,"_",mdat$DNA_well)
head(mdat)


# merge read counts to metadata

colnames(reads) <- c("counts","DNA_plate_well")

df <- merge(mdat, reads)
head(df)


# function to get stats per cohort: 
get_me_stats <- function(DF) {
  out <- DF %>% 
    dplyr::summarise(min=min(counts*2),
                     median=median(counts*2),
                     mean=mean(counts*2),
                     sd=sd(counts*2),
                     max=max(counts*2),
                     samples=n(),
                     sum=sum(counts*2))
  return(out)
}

z <- df %>% 
  group_by(Cohort) %>% 
  get_me_stats()

sink(file = paste0(out_dir,"read_counts_stats.txt"), 
     append = TRUE, type = c("output"))
paste0("##################################") 
paste0("Summary of read counts (R1*2) ",
       summary(as.data.frame(df$counts*2)))
paste0("##################################")
paste0("Standard deviation of read counts : ",
       sd(df$counts*2))
paste0("##################################")
paste0("Sum of read counts (R1*2): ",
       sum(as.data.frame(df$counts*2)))
paste0("##################################")
paste0("number of samples included (these are all - true - samples): ", NROW(df))    
paste0("##################################")
paste0("Included cohorts : ",
       unique(df$Cohort))
paste0("##################################")
paste0("Reads counts summary per cohort (R1*2) : ")
as.data.frame(z)
paste0("##################################")   
sink()


fwrite(x=df, file=paste0(middle_dir,"readcounts.csv"), sep = ",")

