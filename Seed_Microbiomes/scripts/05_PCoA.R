library(cowplot)
#Plot PCoA results using Bray Curtis dissimilarities and Hellinger transformed ASV abundance matrices from non-rarefied phyloseq object
#N.B. Data is subset for each host species##
#Load in phyloseq objects for each host#

# Load Melaleuca .rds file
bacteria_MQ <- readRDS("MQ_phyloseq_16S.rds")
composition <- transform(bacteria_MQ, transform = "log")
meta <- meta(composition)

#this script has 95% confidence ellipses around data points from each grouping 

bray <- ordinate(composition, method = "PCoA", distance = "bray")
A <- plot_ordination(physeq = composition,
                ordination = bray,
                type = "samples",
                color = "SiteType",
                shape = "SiteType") + 
  theme_bw() + theme(legend.position = "right") +
  stat_ellipse(level = 0.95) +
  geom_point(size=1.5)+ ggtitle("Melaleuca quinquenervia - Bacteria")
print(A)


# Load Microlaena .rds file
bacteria_MS <- readRDS("MS_phyloseq_16S.rds")
composition <- transform(bacteria_MS, transform = "log")
meta <- meta(composition)

#this script has 95% confidence ellipses around data points from each grouping 

bray <- ordinate(composition, method = "PCoA", distance = "bray")
B <- plot_ordination(physeq = composition,
                     ordination = bray,
                     type = "samples",
                     color = "SiteType",
                     shape = "SiteType") + 
  theme_bw() + theme(legend.position = "right") +
  stat_ellipse(level = 0.95) +
  geom_point(size=1.5)+ ggtitle("Microlaeana stipoides - Bacteria")
print(B)

# Load Themeda .rds file
bacteria_TT <- readRDS("TT_phyloseq_16S.rds")
composition <- transform(bacteria_TT, transform = "log")
meta <- meta(composition)

#this script has 95% confidence ellipses around data points from each grouping 

bray <- ordinate(composition, method = "PCoA", distance = "bray")
C <- plot_ordination(physeq = composition,
                     ordination = bray,
                     type = "samples",
                     color = "SiteType",
                     shape = "SiteType") + 
  theme_bw() + theme(legend.position = "right") +
  stat_ellipse(level = 0.95) +
  geom_point(size=1.5)+ ggtitle("Themeda triandra - Bacteria")
print(C)

#Run analysis for each host species, separately subset and each for bacteria and fungi and combine
#Construct PCoA combined Figure 3 ##
cowplot::plot_grid(A, B, C, nrow = 3)

#----------------------------------------------#