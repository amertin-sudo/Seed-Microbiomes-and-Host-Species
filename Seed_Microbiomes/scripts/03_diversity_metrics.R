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

# Extract ASV table
asv_table <- otu_table(bacteria)

# Convert to matrix
asv_table <- as(asv_table, "matrix")

# rarecurve expects samples as rows
if (taxa_are_rows(bacteria)) {
  asv_table <- t(asv_table)
}

# Run rarefaction curve
rarecurve(
  asv_table,
  step = 100,
  xlab = "Number of bacterial reads",
  ylab = "ASVs",
  label = FALSE,
  col = "blue",
  xlim = c(0, 15000)
)
# rarefy samples based on above plot. Bacteria = 8000. Fungi = 15000
bacteria_rare<- rarefy_even_depth(bacteria, sample.size = 8000, rngseed = 1)

# extract metadata from file
divMeta <- meta(bacteria_rare)

# estimate diversity measures and add diversity index data to metadata file
divMeta <- cbind(divMeta, estimate_richness(bacteria_rare, measures = c("Observed","Simpson","Shannon")))

write.csv(divMeta, file = "adivMeta_Bacteria_rare.csv", row.names=FALSE)

#Re-analyse using Non-rarefied dataset. This dataset is used for the bacteria due to excessive sample loss with rarefaction.

# extract metadata from file
divMeta_nonrare <- meta(bacteria)

# estimate diversity measures and add diversity index data to metadata file
divMeta_nonrare <- cbind(divMeta_nonrare, estimate_richness(bacteria, measures = c("Observed","Simpson","Shannon")))

write.csv(divMeta_nonrare, file = "adivMeta_Bacteria_nonrare.csv", row.names=FALSE)

# calculate "Observed" means and sems per host species per SiteType for reporting
#Uses column in metadata that has combination category for host species and site type (Restored or Natural). Repeated for fungal dataset and combined in InkScape.
#Need to choose which divMeta is used - divMeta = rarefied, divMeta_nonrare is non rarefied. Non-rarefied used for bacteria.

divMeansHost_SiteType <- aggregate(x = divMeta_nonrare["Observed"], by = divMeta_nonrare["Host_SiteType"], FUN = "mean")
obsSemsHost_SiteType <- aggregate(x = divMeta_nonrare["Observed"], by = divMeta_nonrare["Host_SiteType"], FUN = "std.error")
# add sem values to the divMeans dataframe
divMeansHost_SiteType <- cbind(divMeansHost_SiteType, obsSemsHost_SiteType$Observed)
colnames(divMeansHost_SiteType)[3] <- "ObservedSEM"


# calculate "Simpson" means and sems per host species per SiteType for reporting
simpMeansHost_SiteType <- aggregate(x = divMeta_nonrare["Simpson"], by = divMeta_nonrare["Host_SiteType"], FUN = "mean")
simpSemsHost_siteType <- aggregate(x = divMeta_nonrare["Simpson"], by = divMeta_nonrare["Host_SiteType"], FUN = "std.error")
# add values to the divMeans dataframe
divMeansHost_SiteType <- cbind(divMeansHost_SiteType, simpMeansHost_SiteType$Simpson)
divMeansHost_SiteType <- cbind(divMeansHost_SiteType, simpSemsHost_siteType$Simpson)
colnames(divMeansHost_SiteType)[4] <- "Simpson"
colnames(divMeansHost_SiteType)[5] <- "SimpsonSEM"


# calculate "Shannon" means and sems per host species per SiteType for reporting
shanMeansHost_SiteType <- aggregate(x = divMeta_nonrare["Shannon"], by = divMeta_nonrare["Host_SiteType"], FUN = "mean")
shanSemsHost_SiteType <- aggregate(x = divMeta_nonrare["Shannon"], by = divMeta_nonrare["Host_SiteType"], FUN = "std.error")
# add values to the divMeans dataframe
divMeansHost_SiteType <- cbind(divMeansHost_SiteType, shanMeansHost_SiteType$Shannon)
divMeansHost_SiteType <- cbind(divMeansHost_SiteType, shanSemsHost_SiteType$Shannon)
colnames(divMeansHost_SiteType)[6] <- "Shannon"
colnames(divMeansHost_SiteType)[7] <- "ShannonSEM"



### Diversity Measures BoxPlot By Host species and SiteType Figure 4 . Colours changed in InkScape to match NMDS ###

seed.obs <- ggplot(divMeta_nonrare, aes(x=Host_SiteType, y=Observed, fill=Host_SiteType)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=1, color="black")+
  scale_fill_manual(values=c("#cccccc", "#05ec8d","#cccccc", "#05ec8d", "#cccccc", "#05ec8d", "#cccccc","#05ec8d","#cccccc", "#05ec8d"))+ lims(y=c(0,100))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("Observed ASVs")
print(seed.obs)

seed.sim <- ggplot(divMeta_nonrare, aes(x=Host_SiteType, y=Simpson, fill=Host_SiteType)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=1, color="black")+
  scale_fill_manual(values=c("#cccccc", "#05ec8d","#cccccc", "#05ec8d", "#cccccc", "#05ec8d","#cccccc", "#05ec8d"))+ lims(y=c(0,1))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("Simpson")
print(seed.sim)

seed.shan <- ggplot(divMeta_nonrare, aes(x=Host_SiteType, y=Shannon, fill=Host_SiteType)) + 
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  stat_summary(fun=mean, geom="point", shape=23, size=1, color="black")+
  scale_fill_manual(values=c("#cccccc", "#05ec8d", "#cccccc", "#05ec8d", "#cccccc", "#05ec8d"))+ lims(y=c(0,4.5))+theme(axis.title.y = element_blank())+
  theme(axis.text.y = element_text(size = 14))+
  ggtitle("Shannon")
print(seed.shan)


grid.arrange(seed.obs,seed.sim,seed.shan, nrow=1, ncol=3)

