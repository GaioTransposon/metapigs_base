
###########################################################################################

# construct an empty dataframe to build on 
complete.df <- data.frame(
  var_explained = character(),
  branch_length = character(),
  branch_width = character(),
  blue = character(),
  green = character(),
  red = character(),
  taxa = character(),
  component = character(),
  PC_position = character(),
  file = character()
)

my.files = list.files(guppyout_dir,pattern=".txt.xml.txt")
my.files <- my.files[49:50]
# extract pca ones
#my.files <- my.files[-17] # removing problematic file (won't parse - it's Ja31 anyway) # not an issue anymore 20210812
#my.files <- my.files[-26] # removing problematic file (won't parse - it's Ja31 anyway) # not an issue anymore 20210812

for (textfile in my.files) {
  
  # read in file 
  my.df <- read_csv(file.path(guppyout_dir,textfile), col_names = FALSE)
  
  # extract file name 
  myfilename <- basename(textfile)
  
  startOFnewPC <- grep("<phylogeny rooted", my.df$X1)
  
  z <- split(my.df, cumsum(1:nrow(my.df) %in% startOFnewPC))
  
  PC1 <- as.data.frame(z$`1`)
  PC1$component = "PC1"
  PC1$var_explained <- as.character(z$`1`[2,])
  PC2 <- as.data.frame(z$`2`)
  PC2$component = "PC2"
  PC2$var_explained <- as.character(z$`2`[2,])
  PC3 <- as.data.frame(z$`3`)
  PC3$component = "PC3"
  PC3$var_explained <- as.character(z$`3`[2,])
  PC4 <- as.data.frame(z$`4`)
  PC4$component = "PC4"
  PC4$var_explained <- as.character(z$`4`[2,])
  PC5 <- as.data.frame(z$`5`)
  PC5$component = "PC5"
  PC5$var_explained <- as.character(z$`5`[2,])
  
  my.df <- rbind(PC1,PC2,PC3,PC4,PC5)
  
  # start the grepping of useful tree info 
  mylist <- grep("blue", my.df$X1)
  myend <- lapply(mylist, print)
  
  df <- data.frame(matrix(unlist(myend), nrow=length(myend), byrow=T))
  colnames(df) <- "G"
  head(df)
  mysel <- df %>%
    dplyr::mutate(FF=G-1) %>% # FF is green
    dplyr::mutate(E=FF-1) %>% # E is red
    dplyr::mutate(D=E-1) %>% # D is "color"
    dplyr::mutate(C=D-1) %>% # C is width
    dplyr::mutate(B=C-1) %>% # B is length
    dplyr::mutate(A=B-1) %>% # A is taxa
    dplyr::mutate(placeholder="placeholder")%>%
    pivot_longer(cols=A:G)
  
  mysel <- as.data.frame(mysel)
  row.names(mysel) <- mysel$value
  
  almostthere <- my.df[match(rownames(mysel), rownames(my.df)),] #my.df[match(rownames(mysel), rownames(my.df), nomatch=0),]
  almostthere <- cSplit(almostthere, "X1","</")
  almostthere$X1_2 <- NULL
  
  myprecious <- almostthere %>% 
    #dplyr::filter(X1_1 != '') %>%    # drop empty rows
    dplyr::mutate(key = rep(c('taxa', 'branch_length', 'branch_width', 'color','red','green','blue'), 
                            n() / 7), 
                  id = cumsum(key == 'taxa')) %>% 
    spread(key, X1_1) %>%
    dplyr::select(var_explained,branch_length,branch_width,blue,green,red,taxa,component)
  
  # remove all the symbols derived from the xml format
  myprecious <- as.data.frame(lapply(myprecious, function(y) gsub("<[^>]+>", "", y)))
  
  # round the variation explained by the PC, down to two digits 
  myprecious$var_explained <-round(as.numeric(as.character(myprecious[,1])) * 100,2)
  
  myprecious$PC_position <- paste0(myprecious$blue,"_",myprecious$green,"_",myprecious$red)
  # only two unique combos that make up for green or red (unique(myprecious$PC_position)
  
  # green standing for higher up in PC, red the opposite
  myprecious$PC_position <- gsub("165_194_102", "up",myprecious$PC_position)
  myprecious$PC_position <- gsub("98_141_252", "down",myprecious$PC_position)
  
  myprecious$file <- myfilename

  complete.df <- rbind(
    complete.df, 
    myprecious
  )
  
}


complete <- complete.df

#unique(jplace_df$file)
unique(complete$file)

# clean file name & parse file name 
complete$file <- gsub('_sel.txt.xml.txt', '', complete$file)
complete$file <- gsub('pca_piggies_group_A', 'groupA', complete$file)
complete$file <- gsub('pca_piggies_group_B', 'groupB', complete$file)
complete$file <- gsub('pca_piggies_CTRLNEO', 'groupC', complete$file)
complete$file <- gsub('pca_piggies_NEONEOD', 'groupD', complete$file)
complete$file <- gsub('pca_piggies_NEONEOC', 'groupE', complete$file)
complete$file <- gsub('pca_piggies_CTRLDs', 'groupF', complete$file)
complete$file <- gsub('pca_piggies_CTRLC', 'groupG', complete$file)
complete$file <- gsub('pca_piggies', 'all', complete$file)
complete$file <- gsub('pca_pos_controls', 'pos_tNONE', complete$file)
complete$file <- gsub('^all$', 'all_tALL', complete$file)


unique(complete$file)

complete <- cSplit(complete, "file","_")

colnames(complete)[colnames(complete) == 'file_1'] <- 'sample_type'
colnames(complete)[colnames(complete) == 'file_2'] <- 'guppied_date'

NROW(complete)
unique(complete$sample_type)
unique(complete$guppied_date)

###############

# simplify taxa


complete$taxa_simple <- complete$taxa %<>%
  gsub('\\{|\\}', '', .) %>% # removes curly brackets
  gsub('^_', '', .) %>% # removes the first _
  gsub('[0-9]+', '', .) %>% # removes digits
  gsub('_$', '', .) %>%  # removes the last _
  gsub('.*__', '', .)  # removes everything up to __ (keeping only the most specific)

sort(unique(complete$taxa_simple))


# store unique taxa, per PC, per down/up; keeping largest branch width and shortest branch length:


head(complete)

# some tests
a <- complete %>%
  group_by(var_explained,PC_position,taxa) %>%
  top_n(desc(branch_length), n = 1) %>%
  group_by(var_explained,PC_position,taxa) %>%
  top_n(branch_width, n = 1)
  
complete <- a


###############


colnames(complete)

# Store the minimum necessary info - no branch width included as we need this xml data 
# for the axes description only (as we don't need this info to this purpose, we can do "distinct")
simplified <- complete %>%
  select(sample_type, guppied_date,var_explained,
         component,PC_position,taxa_simple) %>%
  distinct()

unique(simplified$guppied_date)

# save both complete and simplified dataframes 
fwrite(x = complete, file = paste0(middle_dir,"guppy_xml_complete.df"))
fwrite(x = simplified, file = paste0(middle_dir,"guppy_xml_simplified.df"))



# functions to collect taxa associated with PCs:

#############################


# how to use function created below: PC_up(find_PC5(pos_controls))

# what grid.text is extpecting: e.g. : PC_up(find_PC5(df))    in this case (df = pos_controls) and derives from: 
# ```
# pos_controls <- simplified %>%
#   filter(sample_type=="pos") %>%
#   group_split(component) 
# ```
# and has as columns: 
# colnames(pos_controls[[1]])
# [1] "sample_type"   "guppied_date"  "var_explained" "component"     "PC_position"   "taxa_simple"  


######################################################################################################
######################################################################################################

# functions to find PC of interest


find_PC1 <- function(x) {
  PC <- x[[1]] %>%
    group_by(taxa_simple) %>%
    filter(n()==1) %>% # keep only taxa that appear once either up or down in PC, not both
    arrange(PC_position,taxa_simple)
  return(PC)
}
find_PC2 <- function(x) {
  PC <- x[[2]] %>%
    group_by(taxa_simple) %>%
    filter(n()==1) %>% # keep only taxa that appear once either up or down in PC, not both
    arrange(PC_position,taxa_simple)
  return(PC)
}
find_PC3 <- function(x) {
  PC <- x[[3]] %>%
    group_by(taxa_simple) %>%
    filter(n()==1) %>% # keep only taxa that appear once either up or down in PC, not both
    arrange(PC_position,taxa_simple)
  return(PC)
}
find_PC4 <- function(x) {
  PC <- x[[4]] %>%
    group_by(taxa_simple) %>%
    filter(n()==1) %>% # keep only taxa that appear once either up or down in PC, not both
    arrange(PC_position,taxa_simple)
  return(PC)
}
find_PC5 <- function(x) {
  PC <- x[[5]] %>%
    group_by(taxa_simple) %>%
    filter(n()==1) %>% # keep only taxa that appear once either up or down in PC, not both
    arrange(PC_position,taxa_simple)
  return(PC)
}

######################################################################################################
######################################################################################################

# functions to get taxa going up or down the PC

PC_down <- function(x) {
  down <- paste(as.list((x$taxa_simple[x$PC_position=="down"]),"\n"),collapse="\n")
  return(down)
}


PC_up <- function(x) {
  up <- paste(as.list((x$taxa_simple[x$PC_position=="up"]),"\n"),collapse="\n")
  return(up)
}

# interrogate the function by typing 
# PC_up(anydataframeyouwant)


######################################################################################################
######################################################################################################

# functions to get taxa going up or down the PC

get_var <- function(x) {
  var <- paste(as.list((x$var_explained)[1]))   # first item only 
  return(var)
}

# interrogate the function by typing 
# get_var(anydataframeyouwant)


######################################################################################################
######################################################################################################

######################################################################################################
######################################################################################################




# steps before this: 

##### previous steps: 

# 1 # groups for guppy are made in guppy_group.R 

# 2 # guppy is run 

# 3 # .xml conversion to .txt:        <-  MUST RUN THIS BEFORE RUNNING THIS SCRIPT 

# run forester.jar from command line to convert the .xml file to phyloXML - R readable format (.txt) : 
# this way: 
# java -cp /Users/12705859/Downloads/forester_1050.jar 
# org.forester.application.phyloxml_converter -f=dummy file.xml file.txt
# in a loop: 
# for fpath in /Users/12705859/Desktop/metapigs_base/phylosift/guppy/guppy_output/*.xml; 
# do java -cp /Users/12705859/Downloads/forester_1050.jar 
# org.forester.application.phyloxml_converter -f=dummy "$fpath" 
# "/Users/12705859/Desktop/metapigs_base/phylosift/guppy/guppy_output/$(basename "$fpath").txt"; 
# done

# 4 # .xml(txt) files are read in and parsed in guppy_XML_process.R

##### HERE : 

# 1 # .jplace files are read in and parsed

# 2 # bach effect removal 

# 3 # merges metadata

# 4 # principal component are plotted, where xml data is used to completement the taxa underlying the variation


######################################################################################################



library(vcd)
library(summarytools)
library(readr)
library(splitstackshape)
library(dplyr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(data.table)
library(sva)   # this is ComBat , careful not to install COMBAT instead which is another thing
library(openxlsx)
library(gridExtra)
#install.packages("wordspace")
library(wordspace)

library(pheatmap)

#install.packages("taxize")
library(taxize)

library(ggplot2)
library(plotrix)

source_data = "/Users/danielagaio/Gaio/github/metapigs_base/source_data/" # git 
middle_dir = "/Users/danielagaio/Gaio/github/metapigs_base/middle_dir/" # git 
guppyout_dir = "/Users/danielagaio/Desktop/metapigs_base/phylosift/guppy/guppy_output" # local 
out_dir = "/Users/danielagaio/Desktop/metapigs_base/phylosift/guppy/" # local 
out_dir_git = "/Users/danielagaio/Gaio/github/metapigs_base/out/" # git 


###########################################################################################

# manual settings 
removebatcheffect_allowed <- "yes"     # if yes, it removes the batch effect only where detected

######################################################################################################


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

####
# date formatting: 
from = c("2017-01-30",
         "t0","2017-02-01",
         "2017-02-03",
         "2017-02-06","t2","2017-02-08",
         "2017-02-10",
         "t4",
         "2017-02-16","2017-02-17",
         "t6",
         "2017-02-24", 
         "t8",
         "t9",
         "2017-03-06","2017-03-07","2017-03-08","2017-03-09","2017-03-10",
         "2017-08-14", #mock community
         "2018-01-24",  #probiotics - pos controls
         NA) # neg controls

to = c("tM",
       "t0","t0",
       "t1",
       "t2","t2","t2",
       "t3",
       "t4",
       "t5","t5",
       "t6",
       "t7", 
       "t8",
       "t9",
       "t10","t10","t10","t10","t10",
       "tNONE", #mock community
       "tNONE", #probiotics - pos controls
       "tNONE") #neg controls

# replace collection dates (date format) with groups of collection dates (character format)
mdat$collection_date <- plyr::mapvalues(as.character(mdat$collection_date), from, to)
unique(mdat$collection_date)
####

# load breed and bday data 
details <- read_excel(paste0(source_data, "pigTrial_GrowthWtsGE.hlsx.xlsx"),
                      "Piglet details")


# format details
colnames(details)[colnames(details) == 'STIG'] <- 'isolation_source'
colnames(details)[colnames(details) == 'Nursing Dam'] <- 'nurse'
colnames(details)[colnames(details) == 'STIGDAM'] <- 'stig'
colnames(details)[colnames(details) == '...8'] <- 'breed'
details$isolation_source <- gsub("G","",details$isolation_source)
details$isolation_source <- gsub("T","",details$isolation_source)

details <- details %>%
  dplyr::select(isolation_source,BIRTH_DAY,breed,stig,nurse)

###########################################################################################

# load XML data
simplified <- read_csv(paste0(middle_dir,"guppy_xml_simplified.df"))
simplified <- as.data.frame(simplified)

###########################################################################################


jplace_files = list.files(guppyout_dir,pattern=".proj")

# construct an empty dataframe to build on 
jplace_df <- data.frame(
  DNA_plate = character(),
  DNA_well = character(),
  file = character(),
  PC1 = character(),
  PC2 = character(),
  PC3 = character(),
  PC4 = character(),
  PC5 = character(),
  stringsAsFactors = FALSE
)

for (jplace_file in jplace_files) {
  
  # read in file 
  pcadat <- read_csv(file.path(guppyout_dir,jplace_file), col_names = FALSE)
  
  pcadat <- cSplit(pcadat, "X1","_")
  
  pcadat$DNA_plate <- paste0(pcadat$X1_1,"_",pcadat$X1_2)
  pcadat$DNA_well <- pcadat$X1_3
  colnames(pcadat)[1:5] <- c("PC1","PC2","PC3","PC4","PC5")
  pcadat <- pcadat %>%
    dplyr::select(DNA_plate,DNA_well,PC1,PC2,PC3,PC4,PC5)
  
  pcadat$file <- basename(jplace_file)
  
  jplace_df <- rbind(
    jplace_df, 
    pcadat
  )
  
}

# convert PC columns to numeric class 
jplace_df <- jplace_df %>%
  mutate_at('PC1',as.numeric) %>% 
  mutate_at('PC2',as.numeric) %>% 
  mutate_at('PC3',as.numeric) %>% 
  mutate_at('PC4',as.numeric) %>% 
  mutate_at('PC5',as.numeric) 

# clean the file names
jplace_df$file <- gsub('.proj', '', jplace_df$file)
jplace_df$file <- gsub('_sel.txt', '', jplace_df$file)


##############################
##############################

# run guppy_XML_process.R to get simplified df
# (I tried to load it with read.csv, read_csv, read.csv2, would not keep the same format!!!!! grrrrrrr)

##############################
##############################


# function to adjust p-value
padj_function <- function(x, na.rm = FALSE) (p.adjust(x,method="hommel"))

# determine if and where there is a batch effect
checkbatch_before <- jplace_df %>%
  group_by(file) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      PC1=kruskal.test(.$PC1, .$DNA_plate)$p.value,
      PC2=kruskal.test(.$PC2, .$DNA_plate)$p.value,
      PC3=kruskal.test(.$PC3, .$DNA_plate)$p.value,
      PC4=kruskal.test(.$PC4, .$DNA_plate)$p.value,
      PC5=kruskal.test(.$PC5, .$DNA_plate)$p.value,
      batch_removal=paste0("before"),
      stringsAsFactors=FALSE)
  }) %>%
  mutate_at(c("PC1","PC2","PC3","PC4","PC5"),padj_function) 

# picking of groups is based on (adj) p-value of batch effect being <0.05
batch_affected <- checkbatch_before %>%
  dplyr::filter(PC1<0.05|PC2<0.05|PC3<0.05|PC4<0.05|PC5<0.05) %>%
  dplyr::select(file)

to_remove_batch <- jplace_df[jplace_df$file %in% batch_affected$file,]

# splitting into multiple dataframes (by file name)
multiple_DFs <- split( to_remove_batch , f = to_remove_batch$file )

# construct an empty dataframe to build on 
unbatched <- data.frame(
  DNA_plate = character(),
  DNA_well = character(),
  file = character(),
  PC1 = character(),
  PC2 = character(),
  PC3 = character(),
  PC4 = character(),
  PC5 = character(),
  stringsAsFactors = FALSE
)


##############################
# this loop is entered if manually allowed (top of script)


if (removebatcheffect_allowed=="yes") {
  for (single_DF in multiple_DFs) {
    
    DNA_plate <- single_DF$DNA_plate
    DNA_well <- single_DF$DNA_well
    file <- single_DF$file
    
    single_DF<- data.matrix(single_DF[,3:7], rownames.force = NA)
    
    single_DF<-ComBat(dat=t(as.matrix(single_DF)),DNA_plate,mod=NULL)
    
    single_DF <- t(single_DF)
    single_DF <- as.data.frame(single_DF)
    #single_DF <- unfactor(single_DF[])
    single_DF <- cbind(DNA_plate,DNA_well,file,single_DF)
    
    unbatched <- rbind(
      unbatched,
      single_DF
    )
  }
} else {
  print("No batch effect removal allowed")
}

##############################

# check batch effect AFTER batch effect removal 
checkbatch_unbatched <- unbatched %>%
  group_by(file) %>%
  do({
    data.frame(
      sample_size=NROW(.),
      PC1=kruskal.test(.$PC1, .$DNA_plate)$p.value,
      PC2=kruskal.test(.$PC2, .$DNA_plate)$p.value,
      PC3=kruskal.test(.$PC3, .$DNA_plate)$p.value,
      PC4=kruskal.test(.$PC4, .$DNA_plate)$p.value,
      PC5=kruskal.test(.$PC5, .$DNA_plate)$p.value,
      batch_removal=paste0("after"),
      stringsAsFactors=FALSE)
  }) %>%
  dplyr::mutate_at(c("PC1","PC2","PC3","PC4","PC5"),padj_function) 


##############################


# re-join the dataframes (original with unbatched one)

rest <- anti_join(jplace_df,unbatched,by=c("DNA_plate","DNA_well","file"))
jplace_df_final <- rbind(rest,unbatched)

NROW(jplace_df_final)


##############################


# Time to plot! 


###########################################################################################


#settings for plots
theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             plot.title = element_text(),
             axis.title.x=element_text(colour="black",size=8),
             axis.title.y=element_text(colour="black",size=8),
             axis.text.x=element_text(colour="black",size=8),
             axis.text.y=element_text(colour="black",size=8),
             axis.ticks=element_line(colour="black"),
             legend.position="top",
             plot.margin=unit(c(0.3,0.3,0.3,0.3),"line"))

color_legend <- function(x, y, xlen, ylen, main, tiks, colors){
  text(x, y+.6, main, adj=c(0,0), cex=1.3)
  color.legend(x, y, x+xlen, y+ylen/4, legend=tiks, rect.col=colors, cex=0.8)
}
rbow <- rainbow(40, end=0.7, alpha=0.7)

##############################


# add one level of grouping (e.g.: all group_A* files, all timepoints, belong together)

jplace_df_final$group <- jplace_df_final$file


unique(jplace_df_final$group)

jplace_df_final$group <- gsub('pca_piggies_group_A', 'groupA', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_group_B', 'groupB', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_CTRLNEO', 'groupC', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_NEONEOD', 'groupD', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_NEONEOC', 'groupE', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_CTRLDs', 'groupF', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies_CTRLC', 'groupG', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_piggies', 'all', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_all$', 'all_tALL', jplace_df_final$group)
jplace_df_final$group <- gsub('pca_pos_controls', 'pos_tNONE', jplace_df_final$group)


unique(jplace_df_final$group)

jplace_df_final <- cSplit(jplace_df_final, "group","_")
unique(jplace_df_final$group_2)
jplace_df_final <- setnames(jplace_df_final, old = c('group_1','group_2'), new = c('sample_type','guppied_date'))

head(jplace_df_final)
unique(jplace_df_final$sample_type)
unique(jplace_df_final$guppied_date)

##############################

# merge metadata with details (breed,bday,nurse,...)
mdat_deets <- left_join(mdat,details)

# merge metadata with beta diversity data
multi_coggo <- inner_join(jplace_df_final,mdat_deets, by=c("DNA_plate","DNA_well"))
multi_coggo <- multi_coggo %>%
  dplyr::select(DNA_plate,DNA_well,sample_name,isolation_source,collection_date,Cohort,breed,BIRTH_DAY,nurse,stig,sample_type,guppied_date,PC1,PC2,PC3,PC4,PC5)


##############################

# give some order to the variables 
multi_coggo$sample_type <- factor(multi_coggo$sample_type, 
                                  levels=c("pos","all","groupA","groupB",
                                           "groupC","groupD","groupE","groupF","groupG"))

unique(multi_coggo$guppied_date)
multi_coggo$guppied_date <- factor(multi_coggo$guppied_date, 
                                   levels=c("tALL",
                                            "t0",
                                            "t2",
                                            "t4",
                                            "t6",
                                            "t8",
                                            "t9",
                                            "tNONE"))

unique(multi_coggo$guppied_date)
multi_coggo$BIRTH_DAY <- factor(multi_coggo$BIRTH_DAY, 
                                levels=c("2017-01-06", 
                                         "2017-01-07", 
                                         "2017-01-08",
                                         "2017-01-09",
                                         "2017-01-10",
                                         "2017-01-11"))

##############################

# splitting into multiple dataframes (by sample_type name)
unique(multi_coggo$sample_type)

multi_coggo <- split( multi_coggo , f = multi_coggo$sample_type )


##############################


# get dataframes 


##############################
##############################

DF_positive_controls <- as.data.frame(multi_coggo$pos)

##############################
##############################

DF_piggies <- as.data.frame(multi_coggo$all)

unique(DF_piggies$sample_type)
unique(DF_piggies$guppied_date)

DF_piggies <- DF_piggies %>%
  filter(sample_type=="all") 

##############################
##############################

DF_piggies_time <- as.data.frame(multi_coggo$all)
DF_piggies_time_t0 <- DF_piggies_time %>%
  filter(sample_type=="all") %>%
  filter(guppied_date=="t0")
DF_piggies_time_t2 <- DF_piggies_time %>%
  filter(sample_type=="all") %>%
  filter(guppied_date=="t2")
DF_piggies_time_t4 <- DF_piggies_time %>%
  filter(sample_type=="all") %>%
  filter(guppied_date=="t4")
DF_piggies_time_t6 <- DF_piggies_time %>%
  filter(sample_type=="all") %>%
  filter(guppied_date=="t6")
DF_piggies_time_t8 <- DF_piggies_time %>%
  filter(sample_type=="all") %>%
  filter(guppied_date=="t8")
DF_piggies_time_t9 <- DF_piggies_time %>%
  filter(sample_type=="all") %>%
  filter(guppied_date=="t9")
DF_piggies_time <- rbind(
  DF_piggies_time_t0,
  DF_piggies_time_t2,
  DF_piggies_time_t4,
  DF_piggies_time_t6,
  DF_piggies_time_t8,
  DF_piggies_time_t9)

##############################
##############################

groupA <- as.data.frame(multi_coggo$groupA)
groupA_t0 <- groupA %>%
  filter(guppied_date=="t0")
groupA_t2 <- groupA %>%
  filter(guppied_date=="t2")
groupA_t4 <- groupA %>%
  filter(guppied_date=="t4")
groupA_t6 <- groupA %>%
  filter(guppied_date=="t6")
groupA_t8 <- groupA %>%
  filter(guppied_date=="t8")
groupA_t9 <- groupA %>%
  filter(guppied_date=="t9")
groupA <- rbind(
  groupA_t0,
  groupA_t2,
  groupA_t4,
  groupA_t6,
  groupA_t8,
  groupA_t9)

##############################
##############################

groupB <- as.data.frame(multi_coggo$groupB)
groupB_t0 <- groupB %>%
  filter(guppied_date=="t0")
groupB_t2 <- groupB %>%
  filter(guppied_date=="t2")
groupB_t4 <- groupB %>%
  filter(guppied_date=="t4")
groupB_t6 <- groupB %>%
  filter(guppied_date=="t6")
groupB_t8 <- groupB %>%
  filter(guppied_date=="t8")
groupB_t9 <- groupB %>%
  filter(guppied_date=="t9")
groupB <- rbind(
  groupB_t0,
  groupB_t2,
  groupB_t4,
  groupB_t6,
  groupB_t8,
  groupB_t9)

##############################
##############################

groupC <- as.data.frame(multi_coggo$groupC)
groupC_t0 <- groupC %>%
  filter(guppied_date=="t0")
groupC_t2 <- groupC %>%
  filter(guppied_date=="t2")
groupC_t4 <- groupC %>%
  filter(guppied_date=="t4")
groupC_t6 <- groupC %>%
  filter(guppied_date=="t6")
groupC_t8 <- groupC %>%
  filter(guppied_date=="t8")
groupC_t9 <- groupC %>%
  filter(guppied_date=="t9")
groupC <- rbind(
  groupC_t0,
  groupC_t2,
  groupC_t4,
  groupC_t6,
  groupC_t8,
  groupC_t9)

##############################
##############################

groupD <- as.data.frame(multi_coggo$groupD)
groupD_t0 <- groupD %>%
  filter(guppied_date=="t0")
groupD_t2 <- groupD %>%
  filter(guppied_date=="t2")
groupD_t4 <- groupD %>%
  filter(guppied_date=="t4")
groupD_t6 <- groupD %>%
  filter(guppied_date=="t6")
groupD_t8 <- groupD %>%
  filter(guppied_date=="t8")
groupD_t9 <- groupD %>%
  filter(guppied_date=="t9")
groupD <- rbind(
  groupD_t0,
  groupD_t2,
  groupD_t4,
  groupD_t6,
  groupD_t8,
  groupD_t9)

##############################
##############################

groupE <- as.data.frame(multi_coggo$groupE)
groupE_t0 <- groupE %>%
  filter(guppied_date=="t0")
groupE_t2 <- groupE %>%
  filter(guppied_date=="t2")
groupE_t4 <- groupE %>%
  filter(guppied_date=="t4")
groupE_t6 <- groupE %>%
  filter(guppied_date=="t6")
groupE_t8 <- groupE %>%
  filter(guppied_date=="t8")
groupE_t9 <- groupE %>%
  filter(guppied_date=="t9")
groupE <- rbind(
  groupE_t0,
  groupE_t2,
  groupE_t4,
  groupE_t6,
  groupE_t8,
  groupE_t9)

##############################
##############################


groupF <- as.data.frame(multi_coggo$groupF)
groupF_t0 <- groupF %>%
  filter(guppied_date=="t0")
groupF_t2 <- groupF %>%
  filter(guppied_date=="t2")
groupF_t4 <- groupF %>%
  filter(guppied_date=="t4")
groupF_t6 <- groupF %>%
  filter(guppied_date=="t6")
groupF_t8 <- groupF %>%
  filter(guppied_date=="t8")
groupF_t9 <- groupF %>%
  filter(guppied_date=="t9")
groupF <- rbind(
  groupF_t0,
  groupF_t2,
  groupF_t4,
  groupF_t6,
  groupF_t8,
  groupF_t9)

##############################
##############################


groupG <- as.data.frame(multi_coggo$groupG)
groupG_t0 <- groupG %>%
  filter(guppied_date=="t0")
groupG_t2 <- groupG %>%
  filter(guppied_date=="t2")
groupG_t4 <- groupG %>%
  filter(guppied_date=="t4")
groupG_t6 <- groupG %>%
  filter(guppied_date=="t6")
groupG_t8 <- groupG %>%
  filter(guppied_date=="t8")
groupG_t9 <- groupG %>%
  filter(guppied_date=="t9")
groupG <- rbind(
  groupG_t0,
  groupG_t2,
  groupG_t4,
  groupG_t6,
  groupG_t8,
  groupG_t9)

##############################
##############################


# PLOT! 


##############################
##############################

# positive controls

DF_positive_controls$Cohort <- factor(DF_positive_controls$Cohort, 
                                      levels=c("MockCommunity",
                                               "PosControl_D-Scour",
                                               "PosControl_ColiGuard"))

unique(simplified$sample_type)
xmldata <- simplified %>%
  filter(sample_type=="pos") %>%
  group_split(component) 

PC1PC2_pos_controls <- DF_positive_controls %>%
  ggplot(., aes(x=PC1,y=PC2,color=Cohort))+
  geom_point(size=0.5)+
  theme+
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC2 (",get_var(find_PC2(xmldata)),"%)"))+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)+
  scale_color_discrete(drop=FALSE) +
  theme(plot.margin=unit(c(0.2,0.2,2.9,2.9) ,"cm"))
PC3PC4_pos_controls <- DF_positive_controls %>%
  ggplot(., aes(x=PC3,y=PC4,color=Cohort))+
  geom_point(size=0.5)+
  theme+
  xlab(paste0("PC3 (",get_var(find_PC3(xmldata)),"%)"))+
  ylab(paste0("PC4 (",get_var(find_PC4(xmldata)),"%)"))+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)+
  scale_color_discrete(drop=FALSE) +
  theme(plot.margin=unit(c(0.2,0.2,2.9,2.9) ,"cm"))
PC1PC5_pos_controls <- DF_positive_controls %>%
  ggplot(., aes(x=PC1,y=PC5,color=Cohort))+
  geom_point(size=0.5)+
  theme+
  xlab(paste0("PC1 (",get_var(find_PC1(xmldata)),"%)"))+
  ylab(paste0("PC5 (",get_var(find_PC5(xmldata)),"%)"))+
  stat_ellipse(inherit.aes = TRUE, level = 0.80)+
  scale_color_discrete(drop=FALSE) +
  theme(plot.margin=unit(c(0.2,0.2,2.9,2.9) ,"cm"))


pdf(paste0(out_dir,"guppy_pos_controls.pdf"))
### plot PC3PC4 
PC1PC2_pos_controls
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.4, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.9, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
# PC2
grid.text(PC_down(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.3, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC2(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
### plot PC3PC4 
PC3PC4_pos_controls
# PC3
grid.text(PC_down(find_PC3(xmldata)), x = unit(0.4, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC3(xmldata)), x = unit(0.9, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
# PC4
grid.text(PC_down(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.3, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(paste0(get_var(find_PC4(xmldata)),"%"), x = unit(0.1, "npc"), 
          y = unit(0.55, "npc"),
          gp = gpar(fontsize = 8, fontface = "bold"))
grid.text(PC_up(find_PC4(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
### plot PC1PC5 
PC1PC5_pos_controls
# PC1
grid.text(PC_down(find_PC1(xmldata)), x = unit(0.4, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC1(xmldata)), x = unit(0.9, "npc"), 
          y = unit(0.1, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
# PC5
grid.text(PC_down(find_PC5(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.3, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
grid.text(PC_up(find_PC5(xmldata)), x = unit(0.1, "npc"), 
          y = unit(0.8, "npc"),
          gp = gpar(fontsize = 6, fontface = "bold"))
dev.off()
