library(gmodels)

# 
df1 <- df %>%
  dplyr::select(unrooted_pd,bwpd,pc1,pc2,pc3,pc4,pc5,Cohort,collection_date,isolation_source,
                BIRTH_DAY,cross_breed,LINE,maternal_sow,nurse_sow)
NROW(df1)

# for some reasons df1$pc2 is character and not numeric. convert: 
df1$pc2 <- as.numeric(df1$pc2)

# aggregating by avg (unique samples kept)
cols <- 1:7
df1 <- setDT(df1)[, lapply(.SD, mean), by=c(names(df1)[8:15]), .SDcols=cols]
NROW(df1)



df1

v <- piglets_factors2 %>%
  pivot_longer(cols=3:9) %>%
  dplyr::filter(value <= 0.05)
head(v)


ori <- df1 %>%
  pivot_longer(cols=9:15, names_to="variable", values_to="variable_values") 

ori
unique(v$grouping)

# while reading v 
v <- v[1:4,]
head(v)
head(ori)
datalist = list()
for (row in 1:nrow(v)) {
  
  var_to_sel <- v[row,6]
  coll_date <- v[row,1]
  group_sel <- as.data.frame(v[row,3])
  
  test <- ori[ori$collection_date %in% coll_date,]
  test <- test[test$variable %in% var_to_sel,]
  
  test <- test %>%
    dplyr::select(group_sel$grouping,isolation_source,collection_date,variable, variable_values)
  test <- as.data.frame(test)
  
  test <- test %>% 
    group_by(test[1],collection_date,variable) %>% #LINE,variable,maternal_sow,nurse_sow
    summarise(lowCI = ci(variable_values)[2],
              hiCI = ci(variable_values)[3]) 
  
  datalist[[row]] <- test # add it to your list
  

}

big_data = do.call(rbind, datalist)


# must make groupings naming the same between the two dfs: line --> LINE etc
# find way to either disply or show as table 


