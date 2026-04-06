library(plotrix)
library(microbiome)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(phyloseq)
library(vegan)

# # check read depth 
options(max.print=10000)
sort(sample_sums(bacteria))

# Transpose in order to have taxa in columns and species in rows:
asv.table_Trans <- t(bacteria)
#Create rarefaction curves plot
rarecurve(asv.table_Trans, step=100, xlab = "Number of bacterial reads", ylab ="ASVs", label = FALSE, col = "red", xlim=c(0, 40000))
# rarefy samples based on above plot
bacteria_rare<- rarefy_even_depth(bacteria, sample.size = 8000, rngseed = 1)

# extract metadata from file
divMeta <- meta(bacteria_rare)

# estimate diversity measures and add diversity index data to metadata file
divMeta <- cbind(divMeta, estimate_richness(bacteria_rare, measures = c("Observed","Simpson","Shannon")))

write.csv(divMeta, "/pathtoyourworkingdirectory/adivMeta_Bacteria_rare.csv", row.names=FALSE)


# calculate "Observed" means and sems per host species per SiteType for reporting
#Uses column in metadata that has combination category for host species and site type (Restored or Natural). Repeated for fungal dataset and combined in InkScape.

divMeansHost_SiteType <- aggregate(x = divMeta["Observed"], by = divMeta["Host_SiteType"], FUN = "mean")
obsSemsHost_SiteType <- aggregate(x = divMeta["Observed"], by = divMeta["Host_SiteType"], FUN = "std.error")
# add sem values to the divMeans dataframe
divMeansHost_SiteType <- cbind(divMeansHost_SiteType, obsSemsHost_SiteType$Observed)
colnames(divMeansHost_SiteType)[3] <- "ObservedSEM"


# calculate "Simpson" means and sems per host species per SiteType for reporting
simpMeansHost_SiteType <- aggregate(x = divMeta["Simpson"], by = divMeta["Host_SiteType"], FUN = "mean")
simpSemsHost_siteType <- aggregate(x = divMeta["Simpson"], by = divMeta["Host_SiteType"], FUN = "std.error")
# add values to the divMeans dataframe
divMeansHost_SiteType <- cbind(divMeansHost_SiteType, simpMeansHost_SiteType$Simpson)
divMeansHost_SiteType <- cbind(divMeansHost_SiteType, simpSemsHost_siteType$Simpson)
colnames(divMeansHost_SiteType)[4] <- "Simpson"
colnames(divMeansHost_SiteType)[5] <- "SimpsonSEM"


# calculate "Shannon" means and sems per host species per SiteType for reporting
shanMeansHost_SiteType <- aggregate(x = divMeta["Shannon"], by = divMeta["Host_SiteType"], FUN = "mean")
shanSemsHost_SiteType <- aggregate(x = divMeta["Shannon"], by = divMeta["Host_SiteType"], FUN = "std.error")
# add values to the divMeans dataframe
divMeansHost_SiteType <- cbind(divMeansHost_SiteType, shanMeansHost_SiteType$Shannon)
divMeansHost_SiteType <- cbind(divMeansHost_SiteType, shanSemsHost_SiteType$Shannon)
colnames(divMeansHost_SiteType)[6] <- "Shannon"
colnames(divMeansHost_SiteType)[7] <- "ShannonSEM"



### Diversity Measures BoxPlot By Host species and SiteType Figure 4 ###

seed.obs <- ggplot(divMeta, aes(x=Source_SiteType, y=Observed, fill=Source_SiteType)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=1, color="black")+
  scale_fill_manual(values=c("#cccccc", "#05ec8d","#cccccc", "#05ec8d", "#cccccc", "#05ec8d", "#ec0564","#05ec8d","#cccccc", "#05ec8d"))+ lims(y=c(0,100))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("Observed ASVs")
print(seed9.obs)

seed9.sim <- ggplot(divMeta, aes(x=Source_SiteType, y=Simpson, fill=Source_SiteType)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=1, color="black")+
  scale_fill_manual(values=c("#cccccc", "#05ec8d","#cccccc", "#05ec8d", "#cccccc", "#05ec8d","#cccccc", "#05ec8d"))+ lims(y=c(0,1))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("Simpson")
print(seed9.sim)

seed.shan <- ggplot(divMeta, aes(x=Source_SiteType, y=Shannon, fill=Source_SiteType)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=1, color="black")+
  scale_fill_manual(values=c("#cccccc", "#05ec8d", "#cccccc", "#05ec8d", "#cccccc", "#05ec8d"))+ lims(y=c(0,4.5))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("Shannon")
print(seed9.shan)


grid.arrange(seed.obs,seed.sim,seed.shan, nrow=1, ncol=3)

