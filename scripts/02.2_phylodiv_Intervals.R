library(EnvStats)


head(all_padj_Hommel)
out <- all_padj_Hommel %>%
  pivot_longer(cols=3:9) %>%
  dplyr::filter(value <= 0.05) %>%
  dplyr::filter(!grouping=="birth day - Duroc x Landrace")
head(out)

out$grouping <- gsub(pattern = "line","LINE", out$grouping)
out$grouping <- gsub(pattern = "birth day","BIRTH_DAY", out$grouping)


head(df1)
df1$BIRTH_DAY <- as.character(df1$BIRTH_DAY)
df_final <- data.frame()
for (row in 1:nrow(out)) {
  
  var_to_sel <- out[row,6]
  coll_date <- out[row,1]
  group_sel <- as.data.frame(out[row,3])
  
  test <- df1[df1$collection_date %in% coll_date,]

  test <- test %>%
    dplyr::select(group_sel$grouping,isolation_source,collection_date,var_to_sel$name)
  test <- as.data.frame(test)
  
  test$id <- paste0("set_",row)
  test$grouping <- paste0(as.character(group_sel$grouping))
  test$variable=paste0(var_to_sel)
  
  colnames(test) <- c("spec","isolation_source","collection_date","variable_values","id","grouping","variable")
  
  df_final <- rbind(df_final,test)

}


# splitting into multiple dataframes (by set)
multi_df <- split( df_final , f = df_final$id )



plot_me <- function(split_df) {
  
  # df as dataframe
  df0 <- as.data.frame(split_df)
  
  # plot
  return(print(df0 %>% 
    ggplot(., aes(x=spec,y=variable_values))+
      geom_boxplot(lwd=0.2, outlier.size = 0.5)+
      stat_n_text(size = 1.5)+
    coord_flip()+
    ylab(paste0(df0$variable))+
      xlab(paste0(df0$grouping))+
      theme(axis.text.y=element_text(size=4, vjust = 0, angle=68),
            axis.text.x=element_text(size=4),
            axis.title.x=element_text(size=8),
            axis.title.y=element_text(size=8))+
    facet_grid(~collection_date, scale="free")))
}



p1<-plot_me(df_final %>% dplyr::filter(id=="set_1"))
p2<-plot_me(df_final %>% dplyr::filter(id=="set_2"))
p3<-plot_me(df_final %>% dplyr::filter(id=="set_3"))
p4<-plot_me(df_final %>% dplyr::filter(id=="set_4"))
p5<-plot_me(df_final %>% dplyr::filter(id=="set_5"))
p6<-plot_me(df_final %>% dplyr::filter(id=="set_6"))
p7<-plot_me(df_final %>% dplyr::filter(id=="set_7"))
#plot_me(df_final %>% dplyr::filter(id=="set_8"))  # this one is already plotted in a suppl. figure
#plot_me(df_final %>% dplyr::filter(id=="set_9"))  # this one is already plotted in a suppl. figure
#plot_me(df_final %>% dplyr::filter(id=="set_10"))  # this one is already plotted in a suppl. figure


# extra plot: 

out2 <- all_padj_Hommel %>%
  pivot_longer(cols=3:9) %>%
  dplyr::filter(value <= 0.05) %>%
  dplyr::filter(grouping=="birth day - Duroc x Landrace")
head(out2)


p8<-df1 %>%
  dplyr::filter(cross_breed=="Duroc x Landrace") %>%
  dplyr::filter(collection_date=="i2.2") %>%
  dplyr::select(isolation_source,pc2,BIRTH_DAY,collection_date,cross_breed) %>%
  ggplot(., aes(x=BIRTH_DAY,y=pc2))+
  geom_boxplot(lwd=0.2, outlier.size = 0.5)+
  stat_n_text(size = 1.5)+
  coord_flip()+
  theme(axis.text.y=element_text(size=4, vjust = 0, angle=68),
        axis.text.x=element_text(size=4),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8))+
  facet_grid(~collection_date, scale="free")
  


# put them together:

sign_plots <- ggarrange(p1,p2,p3,p4,
          p5,p6,p7,p8,
          ncol = 4, nrow=2)

pdf(paste0(out_dir,"start_factors_pvalues_boxplots.pdf"))
sign_plots
dev.off()

