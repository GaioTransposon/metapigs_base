

library(readr)
library(dplyr)
library(ggplot2)
library(splitstackshape)
library(cowplot)


runinfo <- read_csv("~/Desktop/metapigs_base/runinfo", 
                    col_types = cols(LoadDate = col_character(), 
                                     size_MB = col_number(),
                                     avgLength = col_number(),
                                     ReleaseDate = col_character(), 
                                     bases = col_number(), 
                                     spots = col_number(), 
                                     spots_with_mates = col_number()))




runinfo <- runinfo[!grepl("Run", runinfo$Run),]

unique(runinfo$Submission) 
# "SRA1045474" "SRA1048237" "SRA1045471" "SRA1045481" "SRA901504"  "SRA879779"  "SRA900821" 
# (SRA901504 is 16S)

runinfo <- runinfo %>% 
  filter(!Submission =="SRA901504")



runinfo <- cSplit(runinfo,"ReleaseDate"," ")
unique(runinfo$ReleaseDate_1)

NROW(runinfo)
head(runinfo)

unique(runinfo$TaxID) # "1510822" 
NROW(unique(runinfo$Run)) # 3592
unique(runinfo$Sample) # 911
unique(runinfo$BioSample) # 911
hist(runinfo$spots_with_mates, breaks=100)
hist(runinfo$bases, breaks=100)
hist(runinfo$size_MB, breaks=100)


sum(runinfo$size_MB) # 2556811 Mb == 2,556811 Tb
summary(runinfo$spots_with_mates)
sum(runinfo$spots_with_mates)    # 27233308608 == 27.2 Billion reads 

summary(runinfo$spots_with_mate)

sum(runinfo$size_MB)



Dep201904_06 <- runinfo %>% filter(ReleaseDate_1=="2019-06-18"|ReleaseDate_1=="2019-04-26"|ReleaseDate_1=="2019-04-29")
NROW(unique(Dep201904_06$SampleName))


Dep202002 <- runinfo %>% filter(ReleaseDate_1=="2020-02-20")
NROW(unique(Dep202002$SampleName))


Dep20200226 <- runinfo %>% filter(ReleaseDate_1=="2020-02-26"|ReleaseDate_1=="2020-02-27")
NROW(unique(Dep20200226$SampleName))



plot_Dep201904_06 <- ggplot(Dep201904_06,aes(x=SampleName,y=spots_with_mates, color=ReleaseDate_1))+
  geom_point()+
  theme(axis.text.x=element_blank())


plot_Dep202002 <- ggplot(Dep202002,aes(x=SampleName,y=spots_with_mates, color=ReleaseDate_1))+
  geom_point()+
  theme(axis.text.x=element_blank())


plot_Dep20200226 <- ggplot(Dep20200226,aes(x=SampleName,y=spots_with_mates, color=ReleaseDate_1))+
  geom_point()+
  theme(axis.text.x=element_blank())


plot_grid(
  plot_Dep201904_06,
  plot_Dep202002,
  plot_Dep20200226,
  nrow=3
)

ggplot(runinfo,aes(x=ReleaseDate_1,y=spots_with_mates, color=ReleaseDate_1))+
  geom_boxplot()+
  geom_jitter()+
  theme(axis.text.x=element_blank())

ggplot(runinfo,aes(x=ReleaseDate_1,y=spots_with_mates, color=ReleaseDate_1))+
  geom_boxplot()+
  geom_jitter()+
  theme(axis.text.x=element_blank())+
  ylim(0,1.0e+08)


##########################
