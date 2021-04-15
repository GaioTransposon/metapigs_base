
library(readr)
library(tidyverse)
library(dplyr)
library(robCompositions)
library(microbiome)
library(phyloseq)
library(ggplot2)
library(ape)
library(splitstackshape)
library(magrittr)
library(scales)
library(ggpubr)
library(pheatmap)

library(EnvStats)
library(readxl)
library(data.table)
library(FSA)
library(openxlsx)
library(purrr)
library(broom)

library(cowplot)
library(ggsci)
library(pairwiseCI)

library(readr)
library(dplyr)
library(readxl)
library(splitstackshape)
library(data.table)
library(robCompositions)
library(tidyr)
library(tidyverse)
library(ggbiplot)
library(magrittr)
library(ggpubr)
library(grDevices)
library(colorRamps)
library(EnvStats)
library(corrplot)
library(grid)
library(cowplot)
library(factoextra)
library(broom)
library(openxlsx)
library(ggplot2)


source_data = "/Users/12705859/metapigs_base/source_data/" # git 
middle_dir = "/Users/12705859/metapigs_base/middle_dir/" # git 
out_dir = "/Users/12705859/Desktop/metapigs_base/sortmerna/" # local 
out_dir_git = "/Users/12705859/metapigs_base/out/" # git 
###########################################################################################


# upload sortmerna output (from the original output we retained columns:plate_well, 16Sgene ID, e-value)
# and we filtered out all the 16S rRNA genes below e-30 threshold 
so <- read_table2(paste0(out_dir,"sortmeall_evaluefiltered.tsv"), col_names = FALSE)
so$X1 <- gsub("_S","", so$X1)
so <- so[,1:2]

#number of unique 16S rRNA genes in dataset (filtered with e-value < e-30)
NROW(unique(so$X2))

NROW(so)

so_work <- so

head(so_work)

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

######################################################################

# upload 16S names 

# 1. unzips the file (otherwise too large for github) and places it in out_dir
zipF<- paste0(source_data,"silva-bac-16s-id90_accession_taxonomy.txt.zip")
unzip(zipF,exdir=out_dir)

# 2. now it reads it from out_dir
silva <- read_table2(paste0(out_dir,"silva-bac-16s-id90_accession_taxonomy.txt"), 
                     col_names = FALSE)


######################################################################



# further reducing size: AIM: 

so_work$count <- as.numeric(paste0(1))

NROW(so_work)
so_done <- setDT(so_work)[,.(A = sum(count)), by = 'X1,X2']
NROW(so_done)
# much an improvement! 

#############

# once size of sortme file is reduced, we can parse it: 
so_done_temp <- cSplit(so_done,"X1","_")
so_done_temp$DNA_plate <- paste0(so_done_temp$X1_1,"_",so_done_temp$X1_2)
so_done_temp$DNA_well <- paste0(so_done_temp$X1_3)
head(so_done_temp)

so_done_temp <- so_done_temp %>%
  dplyr::select(DNA_plate,DNA_well,X2,A)

colnames(so_done_temp)[colnames(so_done_temp) == 'X2'] <- 'rRNA16S'
colnames(so_done_temp)[colnames(so_done_temp) == 'A'] <- 'count'
head(so_done_temp)

###############################

# time to merge to metadata!

NROW(so_done_temp)
NROW(mdat)

df <- left_join(so_done_temp,mdat)

NROW(df)
head(df)


# tM <- "2017-01-30"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-01-30",
  replacement = "tM",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-01-31",
  replacement = "t0",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-01",
  replacement = "t0",
  fixed = TRUE)

# t1 <- "2017-02-03"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-03",
  replacement = "t1",
  fixed = TRUE)

# t2 <- "2017-02-06" "2017-02-07" "2017-02-08"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-06",
  replacement = "t2",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-07",
  replacement = "t2",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-08",
  replacement = "t2",
  fixed = TRUE)

# t3 <- "2017-02-10"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-10",
  replacement = "t3",
  fixed = TRUE)

# t4 <- "2017-02-14"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-14",
  replacement = "t4",
  fixed = TRUE)

# t4 <- "2017-02-16" "2017-02-17"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-16",
  replacement = "t5",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-17",
  replacement = "t5",
  fixed = TRUE)

# t6 <- "2017-02-21"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-21",
  replacement = "t6",
  fixed = TRUE)

# t7 <- "2017-02-24"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-24",
  replacement = "t7",
  fixed = TRUE)

# t8 <- "2017-02-28"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-02-28",
  replacement = "t8",
  fixed = TRUE)

# t9 <- "2017-03-03"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-03-03",
  replacement = "t9",
  fixed = TRUE)

# t10 <- "2017-03-06" "2017-03-07" "2017-03-08" "2017-03-09" "2017-03-10"
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-03-06",
  replacement = "t10",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-03-07",
  replacement = "t10",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-03-08",
  replacement = "t10",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-03-09",
  replacement = "t10",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-03-10",
  replacement = "t10",
  fixed = TRUE)

df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2017-08-14",
  replacement = "no-t-pos",
  fixed = TRUE)
df[,6] <- lapply(
  df[,6],
  gsub,
  pattern = "2018-01-24",
  replacement = "no-t-pos",
  fixed = TRUE)


# no-t-neg for negative control
df <- df %>%
  dplyr::mutate(collection_date = if_else(is.na(collection_date), "no-t-neg", collection_date))

# group mothers vs piglets 
df <- df %>%
  dplyr::mutate(is_mom = if_else(collection_date=="tM", "mothers", "piglets"))


unique(df$collection_date)

df$sample <- paste0(df$collection_date,"_",df$isolation_source)

sample_df <- df %>% 
  dplyr::select(DNA_plate,DNA_well,Cohort,collection_date,isolation_source,sample,is_mom) %>%
  dplyr::select(sample,everything()) %>%
  dplyr::filter(!str_detect(collection_date, "^no-t")) %>%
  ungroup() %>%
  group_by(sample) %>%
  slice(1)
sample_df <- as.data.frame(sample_df)
rownames(sample_df) <- sample_df[,1]


# reorder dates 
sample_df$collection_date  = factor(sample_df$collection_date, levels=c("t0",
                                  "t1", 
                                  "t2",
                                  "t3",
                                  "t4",
                                  "t5",
                                  "t6",
                                  "t7",
                                  "t8",
                                  "t9",
                                  "t10",
                                  "tM",
                                  "no-t"))
# reorder cohorts 
sample_df$Cohort  = factor(sample_df$Cohort, levels=c("Control",
                                      "D-Scour", 
                                      "ColiGuard",
                                      "Neomycin",
                                      "Neomycin+D-Scour",
                                      "Neomycin+ColiGuard",
                                      "Mothers"))



###########

# parsing silva and joining info to my df: 

silva <- cSplit(silva,"X3",";")
colnames(silva) <- c("rRNA16S","access","King","Phy","Class","Order","Fam","Species")
silva_edit <- silva

# replace non-taxa with NA
silva_edit[ silva_edit == "uncultured" ] <- NA
silva_edit[ silva_edit == "Family" ] <- NA
silva_edit[ silva_edit == "Order" ] <- NA
silva_edit[ silva_edit == "Subgroup" ] <- NA
silva_edit[ silva_edit == "Incertae" ] <- NA
#
# extract first non-NA taxonomic definition --> to  rightmost column 
silva_edit_rightmost <- as.matrix(silva_edit)[cbind(1:nrow(silva_edit), 
                                                    max.col(!is.na(silva_edit), "last"))] 
silva_new <- cbind(silva_edit,silva_edit_rightmost)

# joining 
new <- inner_join(df,silva_new) 

NROW(new)

###########


# sum all the norm values that fall within same pig,date,taxa_2
df1 <- new %>%
  dplyr::filter(!str_detect(collection_date, "^no-t")) %>%
  dplyr::select(sample,rRNA16S, count, access,King,Phy,Class,Order,Fam,Species,silva_edit_rightmost) %>%
  group_by(sample,rRNA16S, access,King,Phy,Class,Order,Fam,Species,silva_edit_rightmost) %>%
  dplyr::summarise(last_count = sum(count))
head(df1)


tax <- df1 %>%
  ungroup() %>%
  dplyr::select(rRNA16S, access,King,Phy,Class,Order,Fam,Species,silva_edit_rightmost) %>%
  group_by(rRNA16S, access,King,Phy,Class,Order,Fam,Species,silva_edit_rightmost) %>%
  slice(1)
tax <- as.data.frame(tax)
rownames(tax) <- tax$rRNA16S
tax <- as.matrix(tax)


gotu <- df1 %>%
  ungroup() %>%
  dplyr::select(sample,rRNA16S,last_count) %>%
  pivot_wider(names_from=sample, values_from=last_count, values_fill = list(last_count = 0))

gotu <- as.data.frame(gotu)
rownames(gotu) <- gotu$rRNA16S
gotu$rRNA16S <- NULL




gOTU = otu_table(gotu, taxa_are_rows = TRUE)
NCOL(gOTU)
NROW(gOTU)

TAX = tax_table(tax)
NROW(TAX)
NCOL(TAX)

samples = sample_data(sample_df)
NROW(samples)
NCOL(samples)
tail(samples)

which(rownames(sample_df) == colnames(gotu))

carbom <- phyloseq(gOTU,TAX,samples)

summary(colSums(otu_table(gOTU)))
hist(colSums(otu_table(gOTU)), breaks=1000)


############################################################################################################
############################################################################################################
############################################################################################################


carbom <- phyloseq(gOTU,TAX,samples)
carbom
otu_table(carbom)
summary(colSums(otu_table(carbom)))
hist(colSums(otu_table(carbom)))
# removal of samples with low count (any sample with a count lower than 1k)
r <- which(colSums(otu_table(carbom))<10000) 
#r2 <- which(colSums(otu_table(carbom))>80000) 
#r <- cbind(r,r2)
to_remove <- rownames(as.data.frame(r))
carbom_clean <- prune_samples(!(sample_names(carbom) %in% to_remove), carbom)
hist(colSums(otu_table(carbom_clean)))




# BETA div: 

######################


# ORDINATION 

c <- carbom_clean

# SUBSETTING phyloseq obejct
c <- subset_samples(c, (collection_date %in% c("t0","t2","t4","t6","t8","t10","tM")))
# Normalize number of reads in each sample using median sequencing depth.
total = median(sample_sums(c))
standf = function(x, t=total) round(t * (x / sum(x)))
c = transform_sample_counts(c, standf)
sample_variables(c)

carbom_abund <- c 

carbom_abund.ord <- ordinate(carbom_abund, "PCoA")


so_ordination_plot <- plot_ordination(carbom_abund, carbom_abund.ord, 
                                      type="samples", 
                                      color="collection_date") + 
  geom_point(size=1) +
  theme_bw()+
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12), 
        legend.text = element_text(size=12),
        legend.title = element_text(size=14),
        strip.text = element_text(size=14))+
  guides(colour = guide_legend(nrow = 1))+
  theme(legend.position="top")

pdf(paste0(out_dir,"so_phylo_ordination.pdf"))
so_ordination_plot 
so_ordination_plot +
  facet_wrap(~Cohort)
dev.off()

mdat %>%
  group_by(Cohort) %>%
  dplyr::select(isolation_source) %>%
  distinct() %>%
  tally()

#######################################
# collection date variance calc: 

x <- so_ordination_plot$data %>%
  dplyr::select(collection_date,Axis.1,Axis.2)
x <- as.data.frame(x) 
x$collection_date <- as.character(x$collection_date)
head(x)

x %>%
  group_by(collection_date) %>%
  tally()

with(x, tapply(Axis.1, collection_date, var))
with(x, tapply(Axis.2, collection_date, var))


# Compute the analysis of variance
res.aov <- aov(Axis.1 ~ collection_date, data = x)
summary(res.aov)

pairwise.t.test(x$Axis.1, x$collection_date,
                p.adjust.method = "bonferroni")
pairwise.t.test(x$Axis.2, x$collection_date,
                p.adjust.method = "bonferroni")

#######################################


#######################################
# moms vs piglets 

x <- so_ordination_plot$data %>%
  dplyr::select(is_mom,Axis.1,Axis.2)
x <- as.data.frame(x) 

x %>%
  group_by(is_mom) %>%
  tally()

with(x, tapply(Axis.1, is_mom, var))
with(x, tapply(Axis.2, is_mom, var))


# Compute the analysis of variance
res.aov <- aov(Axis.1 ~ is_mom, data = x)
summary(res.aov)

pairwise.t.test(x$Axis.1, x$is_mom,
                p.adjust.method = "bonferroni")
pairwise.t.test(x$Axis.2, x$is_mom,
                p.adjust.method = "bonferroni")


#######################################


# ALPHA div: 

######################


carbom <- carbom_clean


# Prepare function for rarefaction: 

myrarefy_fun <- function(your_phyloseq_obj) {
  
  # removal of samples with low count (any sample with a count lower than 10k)
  r <- which(colSums(otu_table(your_phyloseq_obj))<10000) 
  to_remove <- rownames(as.data.frame(r))
  
  carbom_noFailSamples <- prune_samples(!(sample_names(your_phyloseq_obj) %in% to_remove), your_phyloseq_obj)
  
  # RAREFY
  carbom_rarefied = rarefy_even_depth(carbom_noFailSamples, 
                                      sample.size = min(sample_sums(carbom_noFailSamples)), 
                                      rngseed = 42)
  return(carbom_rarefied)
  
}


# DIVERSITY 

# NORMALIZATION BY RAREFACTION
#carbom <- phyloseq(gOTU,TAX,samples)
# SUBSETTING phyloseq obejct
carbom <- subset_samples(carbom, (collection_date %in% c("t0","t1","t2","t3","t4","t5", "t6","t7","t8","t9","t10","tM")))

# cut out samples with extremely low counts and RAREFY: 
carbom_rarefied <- myrarefy_fun(carbom)

c <- carbom_rarefied

so_diversity_samples <- plot_richness(c, measures=c("Chao1","Shannon", "ACE", "Shannon", 
                                                    "Simpson", "InvSimpson", "Fisher"), 
                                      color="date", x="date") +
  guides(colour = guide_legend(nrow = 1))+
  theme(legend.position="top")

######### plotting above results in a different way: 
# focus on three measures;
# whisker plots instead 

comparison <- data.frame(comparison=c("t2-t0", "t4-t2", "t6-t4", "t8-t6", "t10-t8","tM-t10"))

Chao1_data <- so_diversity_samples$data %>%
  filter(collection_date=="t0"|collection_date=="t2"|collection_date=="t4"|collection_date=="t6"|collection_date=="t8"|collection_date=="t10"|collection_date=="tM") %>%
  filter(variable=="Chao1")
Shannon_data <- so_diversity_samples$data %>%
  dplyr::filter(collection_date=="t0"|collection_date=="t2"|collection_date=="t4"|collection_date=="t6"|collection_date=="t8"|collection_date=="t10"|collection_date=="tM") %>%
  dplyr::filter(variable=="Shannon")
Simpson_data <- so_diversity_samples$data %>%
  dplyr::filter(collection_date=="t0"|collection_date=="t2"|collection_date=="t4"|collection_date=="t6"|collection_date=="t8"|collection_date=="t10"|collection_date=="tM") %>%
  dplyr::filter(variable=="Simpson")


estimates_CI_ttest_fun <- function(xxx) {
  
  # estimates and CIs: 
  apc <- pairwiseCI(value ~ collection_date, data=xxx,
                    method="Param.diff")
  
  # decompose (simply because object won't simply turn into dataframe)
  s <- summary(apc)
  s_estimate <- round(s$estimate,2)
  s_CI <- as.data.frame(s$conf.int)
  s_CI$comparison <- rownames(s_CI)
  s <- cbind(s_CI,s_estimate)
  # and obtain standard error
  s$se <- round(s$estimate-s$lower,2)

  # t-test with Bonferroni adjust: 
  apcTest <- pairwiseTest(value ~ collection_date, data=xxx,
                          method="t.test")
  t <- summary(apcTest, p.adjust.method = "bonferroni")
  
  st <- inner_join(s,t)
  st$variable <- unique(xxx$variable)
  colnames(st) <- c("lower","upper","comparison","estimate","se","p.val.adj","p.val.raw",
                    "group1","group2","variable")
  return(st)
}

Chao1_data_stats <- estimates_CI_ttest_fun(Chao1_data)
Shannon_data_stats <- estimates_CI_ttest_fun(Shannon_data)
Simpson_data_stats <- estimates_CI_ttest_fun(Simpson_data)



# plotting: 
Chao1_data_stats <- Chao1_data_stats %>%
  inner_join(., comparison, by="comparison") %>%
  dplyr::mutate(y.position=c(1500,1550,1600,1650,1700,1750))

Chao1_plot <- Chao1_data %>%
  ggplot(., aes(x=collection_date, y=value, color=collection_date)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=1)+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "time point",
       y = "Chao1") +
  theme_minimal() +
  stat_pvalue_manual(Chao1_data_stats, size=2, label = "p.val.adj")



Shannon_data_stats <- Shannon_data_stats %>%
  inner_join(., comparison, by="comparison") %>%
  dplyr::mutate(y.position=c(6,6.1,6.2,6.3,6.4,6.5))

Shannon_plot <- Shannon_data %>%
  ggplot(., aes(x=collection_date, y=value, color=collection_date)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=1)+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "time point",
       y = "Shannon") +
  theme_minimal() +
  stat_pvalue_manual(Shannon_data_stats, size=2, label = "p.val.adj")


Simpson_data_stats <- Simpson_data_stats %>%
  inner_join(., comparison, by="comparison") %>%
  dplyr::mutate(y.position=c(1.002,1.004,1.006,1.008,1.010,1.012))

Simpson_plot <- Simpson_data %>%
  ggplot(., aes(x=collection_date, y=value, color=collection_date)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size=1)+
  theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "time point",
       y = "Simpson") +
  theme_minimal() +
  stat_pvalue_manual(Simpson_data_stats, size = 2, label = "p.val.adj")

tosave <- ggarrange(Chao1_plot, 
                    Shannon_plot,
                    Simpson_plot,
                    ncol = 3, 
                    nrow=1,
                    labels=c("A","B","C"), 
                    common.legend = TRUE)

ggsave(filename = paste0(out_dir,"so_phylo_diversity_boxplot.pdf"), plot = tosave)


# save stats : 

so_diversity <- rbind(Shannon_data_stats,
                      Simpson_data_stats,
                      Chao1_data_stats)


also_interesting_to_add <- rbind(estimates_CI_ttest_fun(Chao1_data),
                            estimates_CI_ttest_fun(Shannon_data),
                            estimates_CI_ttest_fun(Simpson_data)) %>%
  dplyr::filter(comparison=="tM-t0"|comparison=="tM-t10"|comparison=="t10-t0") 

so_diversity <- rbind(so_diversity %>% dplyr::select(-y.position), also_interesting_to_add)

# save 
fwrite(x=so_diversity, file = paste0(out_dir_git,"sortmerna_diversity_giga.csv"), sep = ",")


######################################################################
######################################################################

z <- so_diversity %>%
  dplyr::filter(p.val.adj<=0.05)
View(z)






