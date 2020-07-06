

# time trend with sortmerna data
# type of sample (library size) normalization: median sequencing depth
# figures not used for paper (metapigs_base)

library(readr)
library(dplyr)
library(splitstackshape)
library(phyloseq)
library(readxl)


setwd("~/Desktop/metapigs_base/sortmerna/")
basedir = "~/Desktop/metapigs_dry/"


# this is the long longed for silva tax dictionary from the sortmerna people
# who parsed the XML file from silva 

# columns are ordered: 
# reference sequence ID, 
# the second being the accession number and 
# the final column is the taxonomy.

silva_dict <- read_delim("silva_ids_acc_tax/silva-bac-16s-id90_accession_taxonomy.txt", 
                                                    "\t", escape_double = FALSE, col_names = FALSE, 
                                                    trim_ws = TRUE)

sortmeall_clean <- read_table2("~/Desktop/metapigs_base/sortmerna/sortmeall_clean.tsv", col_names = FALSE)

# get rid of the evalue now. not necessary anymore now we already filtered
sortmeall_clean$X3 <- NULL

# get a count (1 S16S x sample = 1 count)
NROW(sortmeall_clean)
sortmeall_clean_counts <- sortmeall_clean %>%
  group_by(X1,X2) %>%
  tally()
NROW(sortmeall_clean_counts)

######################################################################

# upload metadata from the gdtb metadata (already has everything) 

gtdb_df0 <- read_csv("gtdb_df0.txt", col_types = cols(pig=col_character()))

gtdb_df0 <- gtdb_df0 %>%
  dplyr::select(pig,date,cohort,weight,pen,breed,birth_day,nurse_mother,mother) %>%
  distinct()

head(gtdb_df0)


# upload metadata for DNA plate and DNA well info 

mdat <- read_excel(paste0(basedir,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)

mdat$pig = as.character(paste0(mdat$isolation_source))

######################################################################

# silva parsing 
head(silva_dict)
colnames(silva_dict) <- c("S16S","accession","full_tax")
all_silva_tags <- unique(silva_dict$S16S)
NROW(all_silva_tags)

# sortme parsing 
head(sortmeall_clean_counts)
colnames(sortmeall_clean_counts) <- c("sample","S16S","count")
NROW(sortmeall_clean_counts)

# let's check if all the tags I have in sortmerna match with the dictionary 
mines <- unique(sortmeall_clean_counts$S16S)
NROW(mines)

NROW(sortmeall_clean_counts)
dfs_joined <- left_join(sortmeall_clean_counts,silva_dict)
NROW(dfs_joined)
dfs_joined$sample

dfs_joined <- cSplit(dfs_joined,"sample","_")

dfs_joined$DNA_plate <- paste0(dfs_joined$sample_1,"_",dfs_joined$sample_2)
dfs_joined$DNA_well <- dfs_joined$sample_3

dfs_joined <- dfs_joined %>%
  dplyr::select(DNA_plate,DNA_well, S16S, count, accession, full_tax) %>%
  dplyr::mutate(sample=paste0(DNA_plate,"_",DNA_well))
NROW(dfs_joined)

mdat <- mdat %>% mutate(sample=paste0(DNA_plate,"_",DNA_well))
plates_wells_tokeep <- mdat$sample

dfs_joined_clean <- subset(dfs_joined, (sample %in% plates_wells_tokeep))
NROW(unique(dfs_joined_clean$sample))
NROW(dfs_joined_clean)


z <- inner_join(mdat,dfs_joined_clean) %>%
  dplyr::select(DNA_plate,DNA_well,pig,S16S,count,accession,full_tax,`*collection_date`,Cohort)
NROW(z)
head(z)


######################################################################
# dates and cohort abbrev conversion 



# tM <- "2017-01-30"
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-01-30", 
  replacement = "tM", 
  fixed = TRUE)
# t0 <- "2017-01-31" "2017-02-01" 
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-01-31", 
  replacement = "t0", 
  fixed = TRUE)
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-02-01", 
  replacement = "t0", 
  fixed = TRUE)

# t1 <- "2017-02-03" 
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-02-03", 
  replacement = "t1", 
  fixed = TRUE)

# t2 <- "2017-02-06" "2017-02-07" "2017-02-08"
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-02-06", 
  replacement = "t2", 
  fixed = TRUE)
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-02-07", 
  replacement = "t2", 
  fixed = TRUE)
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-02-08", 
  replacement = "t2", 
  fixed = TRUE)

# t3 <- "2017-02-10" 
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-02-10", 
  replacement = "t3", 
  fixed = TRUE)

# t4 <- "2017-02-14"
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-02-14", 
  replacement = "t4", 
  fixed = TRUE)

# t4 <- "2017-02-16" "2017-02-17" 
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-02-16", 
  replacement = "t5", 
  fixed = TRUE)
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-02-17", 
  replacement = "t5", 
  fixed = TRUE)

# t6 <- "2017-02-21" 
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-02-21", 
  replacement = "t6", 
  fixed = TRUE)

# t7 <- "2017-02-24" 
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-02-24", 
  replacement = "t7", 
  fixed = TRUE)

# t8 <- "2017-02-28" 
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-02-28", 
  replacement = "t8", 
  fixed = TRUE)

# t9 <- "2017-03-03" 
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-03-03", 
  replacement = "t9", 
  fixed = TRUE)

# t10 <- "2017-03-06" "2017-03-07" "2017-03-08" "2017-03-09" "2017-03-10"
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-03-06", 
  replacement = "t10", 
  fixed = TRUE)
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-03-07", 
  replacement = "t10", 
  fixed = TRUE)
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-03-08", 
  replacement = "t10", 
  fixed = TRUE)
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-03-09", 
  replacement = "t10", 
  fixed = TRUE)
z[8] <- lapply(
  z[8], 
  gsub, 
  pattern = "2017-03-10", 
  replacement = "t10", 
  fixed = TRUE)


# need to substitute cohort names that contain symbols: 
z[9] <- lapply(
  z[9], 
  gsub, 
  pattern = "Neomycin+D-scour", 
  replacement = "NeoD", 
  fixed = TRUE)
z[9] <- lapply(
  z[9], 
  gsub, 
  pattern = "Neomycin+ColiGuard", 
  replacement = "NeoC", 
  fixed = TRUE)
z[9] <- lapply(
  z[9], 
  gsub, 
  pattern = "D-scour", 
  replacement = "DScour", 
  fixed = TRUE)


z[9] <- lapply(
  z[9], 
  gsub, 
  pattern = "2017-08-14", 
  replacement = "no-t-pos", 
  fixed = TRUE)
z[9] <- lapply(
  z[9], 
  gsub, 
  pattern = "2018-01-24", 
  replacement = "no-t-pos", 
  fixed = TRUE)

colnames(z) <- c("DNA_plate","DNA_well","pig","S16S","count","accession","full_tax","date","cohort","sample")

######################################################################

# otu_table 

z$sample <- paste0(z$date,"_",z$pig)
otu_table <- z %>%
  dplyr::select(sample,S16S,count)

NROW(otu_table)
# mean of count from same species and sample
otu_table <- otu_table %>%
  group_by(sample,S16S) %>%
  dplyr::summarize(count = mean(count)) 
NROW(otu_table)

otu_table <- as.data.frame(otu_table)

head(otu_table)
unique(otu_table$sample)
#which(duplicated(otu_table))


# find minium non zero value 
the_minimum <- min(otu_table[,3][which(otu_table[,3]>0)])

otu_table <- otu_table %>% 
  pivot_wider(names_from=sample, values_from=count)


# add pseudocount 
otu_table[is.na(otu_table)] <- the_minimum/10


otu_table <- as.data.frame(otu_table)
rownames(otu_table) <- otu_table$S16S
otu_table$S16S <- NULL

otu_table <- as.matrix(otu_table)


######################################################################

# taxa_table 

taxa_table <- z %>%
  dplyr::select(S16S,accession,full_tax) %>%
  group_by(S16S,accession,full_tax) %>%
  distinct()

taxa_table <- as.data.frame(taxa_table)
taxa_table <- cSplit(taxa_table,"full_tax",";")
taxa_table <- cSplit(taxa_table,"full_tax_6"," ")


taxa_table <- taxa_table %>% mutate(Species = coalesce(full_tax_6_2,full_tax_6_3,full_tax_6_4)) %>%
  dplyr::select(S16S,accession,
                full_tax_1,
                full_tax_2,
                full_tax_3,
                full_tax_4,
                full_tax_5,
                full_tax_6_1,
                Species)

colnames(taxa_table) <- c("S16S","accession","Kingdom","Phylum","Class","Order",
                          "Family","Genus","Species")

unique(taxa_table$Species)
unique(taxa_table$Genus)


taxa_table <- as.data.frame(taxa_table)

rownames(taxa_table) <- taxa_table$S16S
taxa_table <- as.matrix(taxa_table)

######################################################################

# metadata 

sample_df <- z %>%
  dplyr::select(sample,pig,date,cohort) %>%
  group_by(sample) %>%
  distinct()

sample_df <- as.data.frame(sample_df)
rownames(sample_df) <- sample_df[,1]


######################################################################


# create phyloseq object

gOTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(taxa_table)
samples = sample_data(sample_df)


############################################################################################################

# PLOT

######################


# ORDINATION 

# NORMALIZATION BY MEDIAN SEQUENCING DEPTH
carbom <- phyloseq(gOTU,TAX,samples)
# SUBSETTING phyloseq obejct
carbom <- subset_samples(carbom, (date %in% c("t0","t1","t2","t3","t4","t5","t6","t7","t8","t9")))


# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)
# keep only very abundant OTUs
# taking gOTUs that represent at least 3% of the sample and present in at least 40 samples 
carbom_abund <- filter_taxa(carbom, function(x) sum(x > total*0.02) > 20, TRUE)

carbom_abund.ord <- ordinate(carbom_abund, "PCoA", "bray")

sm_ordination_plot <- plot_ordination(carbom_abund, carbom_abund.ord, type="samples", color="date") + 
  geom_point(size=1) +
  theme_bw()+
  theme(axis.title = element_text(size=9),
        axis.text = element_text(size=7))+
  guides(colour = guide_legend(nrow = 1))+
  theme(legend.position="top")

pdf("sortmerna_phylo_ordination.pdf")
sm_ordination_plot
sm_ordination_plot +
  facet_wrap(~cohort)
dev.off()


######################################################################


carbom <- phyloseq(gOTU,TAX,samples)
# SUBSETTING phyloseq obejct
carbom <- subset_samples(carbom, (date %in% c("t0","t1","t2","t3","t4","t5","t6","t7","t8","t9")))


# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(carbom))
standf = function(x, t=total) round(t * (x / sum(x)))
carbom = transform_sample_counts(carbom, standf)
sample_variables(carbom)


# keep only very abundant OTUs
# taking gOTUs that represent at least 3% of the sample and present in at least 40 samples 
carbom_abund <- filter_taxa(carbom, function(x) sum(x > total*0.02) > 20, TRUE)

# HEATMAP time - cohorts

pdf("sortmerna_phylo_heatmap.pdf")
plot_heatmap(carbom_abund, 
             taxa.label = "Genus", 
             sample.order = "date") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "Microbe Species Diversity")
plot_heatmap(carbom_abund, 
             taxa.label = "Genus", 
             sample.order = "date") +
  facet_wrap(~ cohort, switch = "x", scales = "free_x")+
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(label = "Microbe Species Diversity")
dev.off()


######################

