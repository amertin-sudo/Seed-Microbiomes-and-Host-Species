##############################
##     PERMANOVA   ##
##############################

## A PERMANOVA was run on the Hellinger transformed data set with a Bray Curtis distance matrix
library(vegan)
library(phyloseq)
library(microbiome)
library(EcolUtils)

#Compositional (Hellinger) transform abundance data
compositional <- microbiome::transform(bacteria, transform = "Hellinger")

#Generate Bray-Curtis distance matrix
bray_dist_matrix <- phyloseq::distance(compositional, method = "bray") 

#Dispersion test and plot; homogeneity of dispersion among groups is an assumption for adonis

#Bray Dispoersion tests -Betadisp
dispr.bray <- vegan::betadisper(bray_dist_matrix, phyloseq::sample_data(compositional)$Source_SiteType)
dispr.bray
plot(dispr.bray)
anova(dispr.bray) 
#reject the assumption of homogeneity of dispersion by type if p < 0.05

# Run PerMANOVA with stratified/constrained permutations to account for non-independence between site pairs.

# Extract metadata
meta <- as(sample_data(compositional), "data.frame")

# Ensure factors are set correctly
meta$Source <- as.factor(meta$Source)
meta$SiteType <- as.factor(meta$SiteType)  # Restored vs Natural
meta$SitePair <- as.factor(meta$SitePair)  # column to identify site pairs i.e. SitePair1, SitePair2, for each sample.

# PERMANOVA with constrained permutations
adonis2(
  bray_dist_matrix ~ Source * SiteType,
  data = meta,
  permutations = 999,
  strata = meta$SitePair
)

# Pairwise comparisons constrained by site pairs


# Use function by Pedro Martine Zarbizu https://github.com/pmartinezarbizu/pairwiseAdonis
pairwiseAdonis::pairwise.adonis2(bray_dist_matrix ~ Source*SiteType, data = meta, strata = meta$SitePair)

####-----------------------------------------------------------------###
