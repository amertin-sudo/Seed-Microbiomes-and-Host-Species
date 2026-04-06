library(cowplot)
#Plot PCoA results using Bray Curtis dissimilarities and Hellinger transformed ASV abundance matrices from non-rarefied phyloseq object
#N.B. Data is subset for each host species##
composition <- transform(bacteria_subset, transform = "Hellinger")
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
  geom_point(size=1.5)+ ggtitle("Bray")
print(A)

#Run analysis for each host species, separately subset and each for bacteria and fungi and combine using cowplot grid
#Construct PCoA combined Figure 3 ##
cowplot::plot_grid(A, B, C, D, E, F, nrow = 3)

#----------------------------------------------#