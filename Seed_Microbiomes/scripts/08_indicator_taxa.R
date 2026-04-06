
###################### Indicator Species Analysis - Indicspecies ###############################################

library(indicspecies)
library(plyr)
library(phylosmith)
library(microbiome)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(svglite)
remotes::install_github("schuyler-smith/phylosmith")
library(tidyverse)
library(tidytext)
library(patchwork)

###For this to work correctly, the order of your samples MUST MATCH the order of the classifications you assign.
###You are working from the asv table which does not have the corresponding metadata.
###To get around this I pre-organize my data by assigning each treatment a numerical value and sorting to keep each group together.
#Analysis is run on a subset of samples that correspond to a single host species. This is then repeated separately for bacteria and fungal data.


##ASV LEVEL###

# This package does not work with data in Phyloseq format, so we will have to adjust that.
rel <- transform_sample_counts(bacteria_subset, function(x) 100 * x/sum(x))

#order metadata by variable of interest and check metadata table to confirm
rel <- set_sample_order(rel, 'SiteType')
meta <- microbiome::meta(rel)

# Community data matrix: convert from phyloseq object. This asv.table must have no asvs with zero sums.
asv.table <- as(otu_table(rel), "matrix")

# Convert in order to have taxa in columns and species in rows:
asv.table <- t(asv.table)

# save row names to file
sampNames <- row.names(asv.table)

# remove rownames from otu table
row.names(asv.table) <- NULL

# convert otu table to data frame
asv.table <- as.data.frame(asv.table)

# add row names to data frame as column of data
asv.table <- cbind(sampNames, asv.table)


# Defining classification of samples (groups have to be in discrete categories)
count(meta$SiteType)
group <- c(rep("Natural", 56), rep("Restored", 63))

# Performing the indicator value analysis (omit sampNames column in asv table, otherwise indval fails)
indval = multipatt(asv.table[,-1], group, control = how(nperm=999))

# Displaying the data: we can set the thresholds of specificity, sensitivity and significance
summary(indval, At = 0.2, Bt = 0.2, alpha = 0.05, func = "IndVal")

# The indicator value index is the product of two components, called "At" and "Bt"
# Component "A" is the probability that a sample belongs to its target group
# given the fact that the ASV has been found. This conditional probability is called the specificity
# of the ASV as indicator of the sample group. Component "B" is the probability of finding the ASV
# in samples belonging to the sample group. This second conditional probability is called the fidelity
# of the ASV as indicator of the target sample group.

# For easier reporting, transfer the output in a file:
capture.output(summary(indval, invdalcomp = TRUE, At = 0.2, Bt = 0.2, alpha = 0.05), file = "IndvalOutput_hostSpecies.txt")

## Horizontal Bar Chart for Indicator species ##########

#run once for bacteria, again for fungi, then combine together at end.


df_bacteria <- read_csv("indicator_species_Bacterial.csv") #needs to have guild data

# Define consistent guild levels & colours
guild_levels <- c("LichenSymbiont", "PlantSymbiont", "PlantPathogen",
                  "Rhizosphere", "Soil", "Unassigned")

guild_colors <- c(
  "LichenSymbiont" = "#8ff2d1",      # red
  "PlantSymbiont" = "#f59784",      # blue
  "PlantPathogen" = "#8fb0f2", # green
  "Rhizosphere" = "#f2d18f",    # purple
  "Soil" = "#F28fb0",
  "Unassigned" = "darkgrey"
)

df_bacteria$Guild <- factor(df_bacteria$Guild, levels = guild_levels)

plot_list <- list()

for(host in unique(df_bacteria$HostSpecies)) {
  
  df_sub <- df_bacteria %>% filter(HostSpecies == host)
  
  # Order species: SiteType first, then IndVal
  df_sub <- df_sub %>%
    mutate(SiteType = factor(SiteType, levels = c("Restored", "Natural"))) %>%
    arrange(SiteType, IndVal) %>%
    mutate(Species_order = row_number())
  
  # Dotted line between Restored and Natural
  line_pos <- df_sub %>%
    filter(SiteType == "Restored") %>%
    summarise(yint = max(Species_order) + 0.5)
  
  # Midpoint label positions
  label_pos <- df_sub %>%
    group_by(SiteType) %>%
    summarise(ypos = mean(Species_order))
  
  p_bacteria <- ggplot(df_sub, aes(
    x = IndVal,
    y = reorder(Species, Species_order),
    fill = Guild
  )) +
    geom_col() +
    scale_fill_manual(values = guild_colors, drop = FALSE) +  # keep all levels
    scale_x_continuous(limits = c(0, 1.1),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = "IndVal test statistic", y = "Species", fill = "Guild",
         title = host) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 7),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "bottom"
    ) +
    geom_hline(
      data = line_pos,
      aes(yintercept = yint),
      linetype = "dotted",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = label_pos,
      aes(
        x = 1.05,
        y = ypos,
        label = SiteType
      ),
      inherit.aes = FALSE,
      hjust = 0,
      fontface = "bold",
      size = 3.5
    )
  
  plot_list[[host]] <- p_bacteria
}

# Combine plots with a single shared legend
combined_plot_bacteria <- wrap_plots(plot_list, ncol = 1, guides = "collect") &
  theme(legend.position = "right")

print(combined_plot_bacteria)

# --- Combine fungi and bacteria outputs ---
final_plot <- (combined_plot_fungi | combined_plot_bacteria)

print(final_plot)