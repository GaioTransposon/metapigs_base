


library(ggplot2)
library(wordspace)
library(pheatmap)
library(taxize)
library(readxl)
library(dplyr)
library(readr)

source_data = "/Users/12705859/metapigs_base/source_data/" # git 
middle_dir = "/Users/12705859/metapigs_base/middle_dir/" # git 
out_dir = "/Users/12705859/Desktop/metapigs_base/phylosift/guppy/" # local 

###########################################################################################

# visualizing guppy fat output 

# load metadata 
mdat <- read_excel(paste0(source_data,"Metagenome.environmental_20190308_2.xlsx"),
                   col_types = c("text", "numeric", "numeric", "text", "text",
                                 "text", "date", "text","text", "text", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "text", "text","text", "text", "text", "text",
                                 "text","text", "text", "text", "text", "text","text", "text"),
                   skip = 12)

# formatting metadata column names 
mdat$`*collection_date` <- as.character(mdat$`*collection_date`)
colnames(mdat)[colnames(mdat) == '*collection_date'] <- 'collection_date'
colnames(mdat)[colnames(mdat) == '*sample_name'] <- 'sample_name'

mdat$Cohort <- gsub("D-scour","D-Scour", mdat$Cohort)

mdat <- mdat %>%
  dplyr::select(isolation_source,collection_date,Cohort,DNA_plate,DNA_well,PigPen,sample_name)
head(mdat)

# load guppy fat output df

guppyfat_simplified <- read_csv(paste0(middle_dir,"guppyfat_simplified"),col_names = TRUE)


# merge guppy fat output df with metadata 

df <- right_join(guppyfat_simplified,mdat)

df <- df %>% 
  select(name,branch_length,branch_width,isolation_source,collection_date,Cohort,sample_name)


##########################################################################################
##########################################################################################

# POSITIVE CONTROLS 

# number of unique taxa identified in whole dataset from (guppy fat) reads placements onto tree
x <- df %>%
  filter(Cohort=="MockCommunity"|Cohort=="PosControl_D-Scour"|Cohort=="PosControl_ColiGuard")%>%
  group_by(sample_name) %>%
  top_n(n = 15, wt = branch_width) %>% # sensible number to not overcrowd heatmap
  select(name,sample_name)

x <- x %>%
  group_by(sample_name,name) %>%
  tally() %>%
  pivot_wider(names_from = sample_name, values_from = n) %>%
  mutate_all(~replace(., is.na(.), 0))

x <- as.data.frame(x)

# to matrix conversion
rownames(x) <- x[,1]
x[,1] <- NULL

##norm by rows
x_m <- as.matrix(x)
#x_m <- as.matrix(x_m[rowSums(x_m)>1,])
x_m <- normalize.rows(x_m, method = "euclidean", 
                        tol = 1e-6, inplace = TRUE)

guppyfat_pos_controls_plot <- pheatmap(x_m, display_numbers = T,
         #main="Mock community contamination",
         cluster_rows = F, cluster_cols = F, fontsize_number = 3,
         fontsize_row = 3)

pdf(paste0(out_dir,"guppyfat_positive_controls.pdf"))
guppyfat_pos_controls_plot
dev.off()


##########################################################################################
##########################################################################################

# PIGGIES TIME  

# number of unique taxa identified in whole dataset from (guppy fat) reads placements onto tree
x <- df %>%
  dplyr::filter(collection_date=="2017-01-31"|collection_date=="2017-02-07"|collection_date=="2017-02-14"|
           collection_date=="2017-02-21"|collection_date=="2017-02-28"|collection_date=="2017-03-03") %>%
  dplyr::filter(Cohort=="Control"|Cohort=="D-Scour"|Cohort=="ColiGuard"|
           Cohort=="Neomycin"|Cohort=="Neomycin+D-Scour"|Cohort=="Neomycin+ColiGuard") 

# re-order
x$Cohort <- factor(x$Cohort, 
                       levels=c("Control", 
                                "D-Scour", 
                                "ColiGuard",
                                "Neomycin",
                                "Neomycin+D-Scour",
                                "Neomycin+ColiGuard"))

# define the number top n of species to visualize (top n by branch width per time points)
mytops <- 20

###############
head(x)


# PIGGIES TIME  - ALL COHORTS
# PIGGIES TIME  - ALL COHORTS

# 1. sum the branch widths that fall within the same taxon name and collection date 
x1 <- x %>%
  group_by(name,collection_date) %>%
  dplyr::summarise(Sum_branch_width = sum(branch_width)) %>%
  dplyr::filter(!name==".")

# 2. normalize by collection date
x2 <- x1 %>%
  group_by(collection_date) %>%
  dplyr::mutate(Sum_branch_width=Sum_branch_width/sum(Sum_branch_width))

# 3. take top 30 taxa highest in abundance per collection date 
x3 <- x2 %>%
  group_by(collection_date) %>%
  top_n(30) 

# 4. filter out these matches (not taxa)
x4 <- x3 %>%
  dplyr::filter(!grepl("METAGENOME",name)) 

# 5. long to wide format; NA to zeros 
x5 <- x4  %>%
  pivot_wider(names_from = collection_date, values_from = Sum_branch_width) %>%
  mutate_all(~replace(., is.na(.), 0))

# set order of columns by date
x5 <- x5 %>%
  dplyr::select(name,`2017-01-31`,`2017-02-07`,`2017-02-14`,
         `2017-02-21`,`2017-02-28`,`2017-03-03`)

# to matrix conversion
x5 <- as.data.frame(x5)

# dates to timepoint conversion: 
colnames(x5) <- c("name","t0","t2","t4","t6","t8","t10")

rownames(x5) <- x5[,1]
x5[,1] <- NULL
x_m <- as.matrix(x5)

# optional: remove rows with few counts
x_m<- x_m[which(rowSums(x_m) > 0.01),]

guppyfat_time_plot <- pheatmap(x_m, display_numbers = T,angle_col = 0,
         cluster_rows = T, cluster_cols = F, fontsize_number = 6,
         fontsize_row = 6)

pdf(paste0(out_dir,"guppyfat_time.pdf"))
guppyfat_time_plot
dev.off()


# PIGGIES TIME  - BY COHORT
pdf(paste0(out_dir,"guppyfat_time_cohorts.pdf"))
cohorts <- unique(x$Cohort)
for (coho in cohorts) {
  # 1. sum the branch widths that fall within the same taxon name and collection date 
  x1 <- x %>%
    filter(Cohort==coho) %>%
    group_by(name,collection_date) %>%
    dplyr::summarise(Sum_branch_width = sum(branch_width)) %>%
    filter(!name==".")
  
  # 2. normalize by collection date
  x2 <- x1 %>%
    group_by(collection_date) %>%
    mutate(Sum_branch_width=Sum_branch_width/sum(Sum_branch_width))
  
  # 3. take top 10 taxa highest in abundance per collection date 
  x3 <- x2 %>%
    group_by(collection_date) %>%
    top_n(10) 
  
  # 4. filter out these matches (not taxa)
  x4 <- x3 %>%
    filter(!grepl("METAGENOME",name)) 
  
  # 5. long to wide format; NA to zeros 
  x5 <- x4  %>%
    pivot_wider(names_from = collection_date, values_from = Sum_branch_width) %>%
    mutate_all(~replace(., is.na(.), 0))
  
  # set order of columns by date
  x5 <- x5 %>%
    select(name,`2017-01-31`,`2017-02-07`,`2017-02-14`,
           `2017-02-21`,`2017-02-28`,`2017-03-03`)
  
  # to matrix conversion
  x5 <- as.data.frame(x5)
  rownames(x5) <- x5[,1]
  x5[,1] <- NULL
  x_m <- as.matrix(x5)
  
  # optional: remove rows with few counts
  x_m<- x_m[which(rowSums(x_m) > 0.02),]
  
  print(pheatmap(x_m, display_numbers = T,
           cluster_rows = F, cluster_cols = F, fontsize_number = 5,
           fontsize_row = 5, main = paste0(coho)))
}
dev.off()




##########################################################################################
##########################################################################################

