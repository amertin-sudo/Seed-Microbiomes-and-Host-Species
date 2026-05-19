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
genus$Genus[genus$Abundance < 2]<- " < 2% Abundance"
write.csv(genus, "genus_barplot.csv", row.names=FALSE)

#How many levels in Genus
HowMany <- length(levels(as.factor(genus$Genus)))

#Plot - N.B. colours were changed in Inkscape. Genus IDs (when ony family was provided) were provided by sequence similarity for that ASV against UNITE v8 database.
ggplot(genus, aes(x = Sample, y = Abundance, fill = Genus)) +  
  geom_bar(stat = "identity") +
  theme(legend.position="right", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)) +  geom_col(width = 0.1) +
  ylab("Relative Abundance of Bacterial Genera \n") +  
  scale_fill_manual(values = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#1B9E77", "#E7298A","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#FFFFFF", "#6A3D9A","#FFFF99","#FFFFFF", "#6A3D9A","#FFFF99","#FFFFFF", "#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#1B9E77", "#E7298A"))


  
#-----------------------------------------------------------------------------##