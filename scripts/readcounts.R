
library(readr)
library(splitstackshape)
library(readxl)
library(dplyr)


source_data = "/Users/12705859/metapigs_base/source_data/" # git 
middle_dir = "/Users/12705859/metapigs_base/middle_dir/" # git 
out_dir = "/Users/12705859/Desktop/metapigs_base/reads_counts_distribution/" # local 

###########################################################################################

# raw read counts (only R1)
readcounts_samples <- read_csv(paste0(source_data,"readcounts_samples.csv"),
                               col_names = FALSE)
tail(readcounts_samples)

######################################################################

# clean paired reads 
lib <- read_delim(paste0(middle_dir,"clean_paired_lib_sizes_final.tsv"),
                  "\t", escape_double = FALSE, col_names = FALSE,trim_ws = TRUE)

hist(lib$X2, breaks=100)

summary(lib$X2)
head(lib)
sum(lib$X2)
NROW(lib)
# this includes all piglet samples 
tail(lib)


######################################################################


# load metadata 
mdat <- read_excel(paste0(source_data,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)

mdat$`*collection_date` <- as.character(mdat$`*collection_date`)
mdat$Cohort <- gsub("Sows","Sows",mdat$Cohort)
mdat$Cohort <- gsub("D-scour","D-Scour", mdat$Cohort)


colnames(mdat)

colnames(mdat)[colnames(mdat) == '*collection_date'] <- 'collection_date'

mdat <- mdat %>%
  dplyr::select(DNA_plate,DNA_well,Cohort,collection_date,isolation_source)

mdat$DNA_plate_well <- paste0(mdat$DNA_plate,"_",mdat$DNA_well)
head(mdat)

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
                  ylab = "Frequency")

pdf(paste0(out_dir, "readcounts_distribution.pdf"), onefile = TRUE)
histogram
dev.off()

pdf("/Users/12705859/Desktop/metapigs_base/readcounts_distribution_all_and_low.pdf", onefile = TRUE)
par(mfrow=c(1,1))
hist(readcounts_samples$X3, breaks = 200, cex.axis=1,cex.lab=1.5, cex.main=2,
     main = "Read count distribution across samples",
     xlab = "read counts",
     ylab = "Frequency")
par(fig=c(0.45, 0.99, 0.2, 1), new = T) 
hist(readcounts_samples$X3[
  readcounts_samples$X3<=30000], 
  breaks=20, 
  xlim=c(0,30000),
  xlab = NULL,
  ylab = NULL,
  main= NULL)
dev.off()


######################################################################

# merge read counts to metadata

head(mdat)

colnames(reads) <- c("counts","DNA_plate_well")

df <- merge(mdat, reads)
head(df)
NROW(df)

summary(df$counts*2) # median=32,688,566  mean=35,118,723
sd(df$counts*2) # sd=20333765
sum(df$counts*2) # sum=31,993,156,356


unique(df$Cohort)

get_me_stats <- function(DF) {
  out <- DF %>% 
    dplyr::summarise(min=min(counts*2),
                  median=median(counts*2),
                  mean=mean(counts*2),
                  sd=sd(counts*2),
                  max=max(counts*2),
                  n=n(),
                  sum=sum(counts*2))
  return(out)
}

# neg controls
neg <- df %>% dplyr::filter(Cohort=="NegativeControl") #%>% dplyr::summarise(sum=sum(counts*2))

# pos controls
pos <- df %>% dplyr::filter(Cohort=="MockCommunity"|
                       Cohort=="PosControl_D-Scour"|
                       Cohort=="PosControl_ColiGuard")
# mothers alone
moms <- df %>% dplyr::filter(Cohort=="Mothers")  

# piggies alone
piggies <- df %>% dplyr::filter(!Cohort=="NegativeControl" & !Cohort=="MockCommunity"
                                & !Cohort=="PosControl_D-Scour" & !Cohort=="PosControl_ColiGuard"
                                & !Cohort=="Mothers")


get_me_stats(df)
get_me_stats(neg)
get_me_stats(pos)
get_me_stats(moms)
get_me_stats(piggies)



######################################################################

# log of reads extracted with sortmerna:

sortmerna_log <- read_csv(paste0(source_data,"sortmerna_log"), 
                          col_names = FALSE)

head(sortmerna_log)

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
  dplyr::select(file, passed_7, passed_8, failed_7, failed_8)

colnames(grouped2) <- c("DNA_plate_well","passed_count","passed_perc","failed_count","failed_perc")


sortmerna_filtering_stats <- grouped2

summary(sortmerna_filtering_stats)


sortme <- sortmerna_filtering_stats 

merged <- merge(df,sortme) %>%
  dplyr::mutate(perc_16S = (passed_count/counts)*100) 

head(merged)
NROW(merged)

hist(merged$perc_16S)
summary(merged$perc_16S)
sum(merged$counts*2) # tot umber of reads
sum(merged$passed_count) # tot number of 16S genes that passed E-value threshold 
sum(merged$failed_count) # tot number of 16S genes that failed E-value threshold 
head(merged)

# 95% confidence interval 
t.test(merged$perc_16S)

# 
# cat all_plates.tsv | cut -f 23 > perc_conf.tsv
# 
# awk '$1>0.1{c++} END{print c+0}' perc_conf.tsv
# 46885952/60586929*100-100
# 
# awk '$1>0.3{c++} END{print c+0}' perc_conf.tsv
# 32959839/60586929*100
# 
# awk '$1>0.9{c++} END{print c+0}' perc_conf.tsv
# 10852863/60586929*100
# 
# awk '$1>0.99{c++} END{print c+0}' perc_conf.tsv
# 4183468/60586929*100

