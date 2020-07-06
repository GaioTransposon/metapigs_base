
setwd("/Users/12705859/Desktop/metapigs_base/reads_counts_distribution/")

library(readr)
library(splitstackshape)

readcounts_samples <- read_csv("readcounts_samples.csv", 
                               col_names = FALSE)


# split first column to contain plate and well
reads <- cSplit(readcounts_samples, "X1", ".")

reads$X2 <- NULL
reads$X1_2 <- NULL
reads$X1_3 <- NULL
reads$X1_4 <- NULL

head(reads)

summary(reads)
sum(reads$X3)

# 15999031919
# nearly 16 billion reads! 

lowreads_samples <- hist(head(reads$X3)[
  readcounts_samples$X3<=40000], 
  breaks=20, 
  xlim=c(0,40000),
  xlab = "read counts",
  ylab = "Frequency",
  main= NULL)
lowreads_samples


pdf("/Users/12705859/Desktop/metapigs_base/readcounts_distribution.pdf", onefile = TRUE)
histogram <- hist(readcounts_samples$X3, breaks = 200,
     main = "Read counts distribution across samples",
     xlab = "read counts",
     ylab = "Frequency")
histogram
dev.off()

pdf("/Users/12705859/Desktop/metapigs_base/readcounts_distribution_all_and_low.pdf", onefile = TRUE)
par(mfrow=c(1,1))
hist(readcounts_samples$X3, breaks = 200,
     main = "Read counts distribution across samples",
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

# 
# 
# a <- subset(readcounts_samples, X3 <= 5000)
# a #10 samples with less than 5000 reads
# 
# a <- subset(readcounts_samples, X3 <= 10000)
# a #25 samples with less than 10000 reads
# 
# a <- subset(readcounts_samples, X3 <= 20000)
# a #43 samples with less than 20000 reads
# 
# a <- subset(readcounts_samples, X3 <= 50000)
# a #53 samples with less than 50000 reads
# 
# 
# # are these negative controls? 
# a <- subset(readcounts_samples, X3 <= 10000)
# a #25 samples with less than 10000 reads
# 
# 
# 
# 
# 
# 
# 
# # import NCBI table
# library(readxl)
# metadata <- read_excel("/Users/12705859/metapigs/source_data/Metagenome.environmental_20190308_2.xlsx",
#                        sheet = 1, 
#                        col_names = TRUE,
#                        col_types = NULL,
#                        skip = 12)
# 
# 
# library(splitstackshape)
# # split sample column to date and pig
# metadata <- cSplit(metadata, "*sample_name", "/")
# 
# # rename columns
# names(metadata)[names(metadata) == "*sample_name_1"] <- "date"
# names(metadata)[names(metadata) == "*sample_name_2"] <- "pig"
# 
# # keep only necessary columns 
# metadata <- metadata[,c(21,23,24,31,32)]
# 
# metadata$X1_1 <- paste0(metadata$DNA_plate,"_",metadata$DNA_well)
# 
# head(reads)
# head(metadata)
# 
# 
# df <- merge(metadata, reads, by = "X1_1")
# head(df)
# 
# 
# # subset rows of negative controls
# negs <- subset(df, grepl(paste("neg.control"), df$pig))
# negs
# 
# 
# summary(negs$X3)




