library(nlme)
library(emmeans)
library(dplyr)
library(ggplot2)

###Diversity ANOVAs###

##Need to first run "01_import data," "02_decontam,", and "03_diversity_metrics" 
#R scripts to get divMeta file

#Observed ASVs_ANOVA
Obs <- lm(Observed ~ Source*SiteType, data = divMeta_nonrare, na.action = na.omit)
summary(Obs)
anova(Obs)
#Pairwise tests
pairs(emmeans(Obs, "SiteType", by = "Source"))
#pairs(emmeans(Obs, "Source")) # run this for when only one factor significant

#Simpsons_ANOVA
Sim <- lm(Simpson ~ Source*SiteType, data = divMeta_nonrare, na.action = na.omit)
summary(Sim)
anova(Sim)
#Pairwise tests
pairs(emmeans(Sim, "SiteType", by = "Source")) #run this for when interaction term significant
#pairs(emmeans(Sim, "Source")) # run this for when only one factor significant

#Shannon_ANOVA
Shan <- lm(Shannon ~ Source*SiteType, data = divMeta_nonrare, na.action = na.omit)
summary(Shan)
anova(Shan)

#Pairwise tests
pairs(emmeans(Shan, "SiteType", by = "Source"))
#pairs(emmeans(Shan, "Source"))

rm(Obs)
rm(Shan)
rm(Sim)

