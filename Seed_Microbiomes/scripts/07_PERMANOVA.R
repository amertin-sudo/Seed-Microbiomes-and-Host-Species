##############################
##     PERMANOVA   ##
##############################

## A PERMANOVA was run on the log transformed data set with a Bray Curtis distance matrix
library(vegan)
library(phyloseq)
library(microbiome)
library(EcolUtils)

#Compositional (Hellinger) transform abundance data
compositional <- microbiome::transform(bacteria, transform = "log")

#Generate Bray-Curtis distance matrix
bray_dist_matrix <- phyloseq::distance(compositional, method = "bray") 

#Dispersion test and plot; homogeneity of dispersion among groups is an assumption for adonis

#Bray Dispersion tests -Betadisp
#reject the assumption of homogeneity of dispersion by type if p < 0.05

# Extract metadata
meta <- as(sample_data(compositional), "data.frame")

meta$Source <- as.factor(meta$Source)
meta$SiteType <- as.factor(meta$SiteType)
meta$Source_SiteType <- interaction(meta$Source, meta$SiteType)

####################################
# Dispersion tests
####################################

# Source
disp_source <- betadisper(bray_dist_matrix, meta$Source)
anova(disp_source)
permutest(disp_source, permutations = 999)

# SiteType
disp_sitetype <- betadisper(bray_dist_matrix, meta$SiteType)
anova(disp_sitetype)
permutest(disp_sitetype, permutations = 999)

# Interaction grouping
disp_interaction <- betadisper(bray_dist_matrix, meta$Source_SiteType)
anova(disp_interaction)
permutest(disp_interaction, permutations = 999)


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
pairwiseAdonis::pairwise.adonis2(bray_dist_matrix ~ Host_SiteType, data = meta)

####-----------------------------------------------------------------###
