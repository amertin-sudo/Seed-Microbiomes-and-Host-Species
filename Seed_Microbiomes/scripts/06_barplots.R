library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("microbiome")
library("svglite")

##########################################
##                Bar plots             ##
##########################################


# Merge samples by HostSpecies and SiteType. There must be a column named this in the metadata

Host_SiteType<- merge_samples(bacteria, "Host_SiteType")

## Genus level relative abundance barplot - Figure 2 ##
# To represent at the Genus level 
Genus <- tax_glom(Host_SiteType,taxrank = "Genus")

# Transform counts in relative abundance and select most abundant families
Genus <- transform_sample_counts(Genus, function(x) 100 * x/sum(x))
genus <- psmelt(Genus)
genus$Genus <- as.character(genus$Genus)

#rename and pool Genera with <1% abundance
genus$Genus[genus$Abundance < 1]<- " < 1% Abundance"
write.csv(genus, "genus_barplot.csv", row.names=FALSE)

#How many levels in Genus
HowMany <- length(levels(as.factor(genus$Genus)))

#Plot
ggplot(genus, aes(x = Sample, y = Abundance, fill = Genus)) +  
  geom_bar(stat = "identity") +
  theme(legend.position="none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(vars(Genus)) +
  ylab("Relative Abundance of Bacterial Genera")
  scale_fill_manual(values = c("#D3D3D3","#1F78B4","#B2DF8A","#FFFF00","#FB9A99","#E31A1C","#808080", "#E7298A","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99", "#000000", "#E69200", "#56B4E9", "#009E73", "#FF6347", "#FFD700", "#D55F00", "#CB79A7", "#CC6583", "#FE34B4", "#EE23C1", "#FFD700", "#000566", "#FFCCCC", "#FF66FF","#FFDB6D", "#C4961A", "#00AFBB", "#009999", "#A6CEE3","#1F78B4","#B2DF8A","#FFFF00","#FB9A99","#E31A1C","#808080", "#E7298A","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99", "#000000", "#E69200", "#56B4E9", "#009E73", "#FF6347", "#FFD700", "#D55F00", "#CB79A7", "#CC6583", "#FE34B4", "#EE23C1", "#FFD700", "#000566", "#FFCCCC", "#FF66FF","#FFDB6D", "#C4961A", "#00AFBB", "#009999", "#A6CEE3","#1F78B4","#B2DF8A","#FFFF00","#FB9A99","#E31A1C","#808080", "#E7298A","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99", "#000000", "#E69200", "#56B4E9", "#009E73", "#FF6347", "#FFD700", "#D55F00", "#CB79A7", "#CC6583", "#FE34B4", "#EE23C1", "#FFD700", "#000566", "#FFCCCC", "#FF66FF","#FFDB6D", "#C4961A", "#00AFBB", "#009999"))

#-----------------------------------------------------------------------------##