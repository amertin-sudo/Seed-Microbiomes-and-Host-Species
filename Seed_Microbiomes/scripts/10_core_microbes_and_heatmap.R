
# Identify core microbial taxa using a prevalence-abundance framework
# Requires: phyloseq (optional), tidyverse, microbiome (optional), pheatmap
#
# Core definition used here:
# - A taxon is considered "core" if:
#     prevalence >= prevalence_threshold    (prevalence = percent of samples where taxon >= detection_threshold)
# AND mean_relative_abundance >= abundance_threshold
# This dual threshold reduces calling very-rare but widespread taxa or very-abundant but rare taxa "core".

# Input: either provide a phyloseq object file (RDS), or OTU table + metadata CSVs.Need to subset phyloseq per Host species.

phyloseq_rds <- "/pathtords/sample.rds"  # e.g. "data/my_phyloseq.rds" OR set to NULL to use CSV inputs
otu_csv <- "data/otu_table.csv"    # rows = taxa, columns = samples OR opposite (script will detect)
meta_csv <- "data/sample_metadata.csv" # sample metadata, with column "SampleID" matching sample names

# Thresholds (customise)
prevalence_threshold <- 0.5   # fraction of samples (0-1) : present in >=50% of samples
detection_threshold  <- 0.0001 # relative abundance detection threshold (e.g. 0.0001 = 0.01%)
abundance_threshold  <- 0.001  # mean relative abundance threshold (e.g. 0.001 = 0.1%)

# Output paths
out_dir <- "output_dir"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Load libraries
# ----------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)    
  library(microbiome)  
  library(pheatmap)
})


# Function to compute prevalence and mean relative abundance
compute_prevalence_abundance <- function(relab_mat, detection_thresh = detection_threshold) {
  # relab_mat: taxa x samples (rows taxa, columns samples) numeric matrix of relative abundances
  # returns a tibble: Taxon, prevalence (0-1), mean_rel_abundance
  taxa <- rownames(relab_mat)
  present_mat <- relab_mat >= detection_thresh
  prevalence <- rowSums(present_mat, na.rm = TRUE) / ncol(relab_mat)
  mean_abund <- rowMeans(relab_mat, na.rm = TRUE)
  tibble(
    Taxon = taxa,
    Prevalence = prevalence,
    MeanRelAbundance = mean_abund
  )
}

# ----------------------------
# Load data (phyloseq or CSV)
# ----------------------------
if (!is.null(phyloseq_rds) && file.exists(phyloseq_rds)) {
  message("Loading phyloseq object from RDS...")
  ps <- readRDS(phyloseq_rds)
  if (!inherits(ps, "phyloseq")) stop("RDS does not contain a phyloseq object.")
  # Extract OTU table
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps) == FALSE) {
    # phyloseq stores taxa as columns sometimes; ensure taxa are rows
    otu <- t(otu)
  }
  sample_meta <- as.data.frame(sample_data(ps))
  sample_meta$SampleID <- rownames(sample_meta)
  taxa_meta <- if (!is.null(tax_table(ps, errorIfNULL = FALSE))) as.data.frame(tax_table(ps)) else NULL
} else {
  message("Loading OTU table and metadata from CSVs...")
  if (!file.exists(otu_csv)) stop("OTU CSV not found: ", otu_csv)
  if (!file.exists(meta_csv)) stop("Metadata CSV not found: ", meta_csv)
  otu_df <- read_csv(otu_csv, col_types = cols())
  meta_df <- read_csv(meta_csv, col_types = cols())
  # Expect first column is Taxon or feature ID OR first row are taxa names.
  # Detect orientation: if first column is "Taxon" or "FeatureID", assume rows = taxa
  if ("Taxon" %in% colnames(otu_df) | "FeatureID" %in% colnames(otu_df) | "ASV" %in% colnames(otu_df)) {
    # taxa in first column
    first_col_name <- colnames(otu_df)[1]
    otu_mat <- otu_df %>%
      rename(Taxon = !!sym(first_col_name)) %>%
      column_to_rownames("Taxon") %>%
      as.matrix()
  } else {
    # Maybe rows are samples and columns taxa: detect by matching sample IDs with metadata
    # If a column in otu_df matches SampleID in metadata, assume samples are rows
    if ("SampleID" %in% colnames(meta_df)) {
      sample_ids <- meta_df$SampleID
      common <- intersect(sample_ids, colnames(otu_df))
      if (length(common) > 0) {
        # columns include sample IDs -> assume taxa are rows with rownames in first column
        # If there is no explicit taxa column, create sequential IDs
        if (!("Taxon" %in% colnames(otu_df))) {
          otu_df <- otu_df %>% mutate(Taxon = row_number())
          otu_df <- otu_df %>% select(Taxon, everything())
        }
        otu_mat <- otu_df %>% rename(Taxon = !!sym(colnames(otu_df)[1])) %>%
          column_to_rownames("Taxon") %>% as.matrix()
      } else {
        # fallback: treat first column as taxa
        otu_mat <- otu_df %>%
          rename(Taxon = !!sym(colnames(otu_df)[1])) %>%
          column_to_rownames("Taxon") %>%
          as.matrix()
      }
    } else {
      # fallback: treat first column as taxa
      otu_mat <- otu_df %>%
        rename(Taxon = !!sym(colnames(otu_df)[1])) %>%
        column_to_rownames("Taxon") %>%
        as.matrix()
    }
  }
  otu <- otu_mat
  sample_meta <- meta_df
  if (!("SampleID" %in% colnames(sample_meta))) {
    # try to infer sample names from rownames
    sample_meta$SampleID <- sample_meta[[1]]
  }
  taxa_meta <- NULL
}

# Ensure sample names align between otu and metadata
sample_names_in_otu <- colnames(otu)
if (is.null(sample_names_in_otu)) stop("Cannot read sample names from OTU table.")
# If metadata has SampleID column, use it to align; otherwise assume metadata rows correspond to sample order
if ("SampleID" %in% colnames(sample_meta)) {
  missing_meta <- setdiff(sample_names_in_otu, sample_meta$SampleID)
  if (length(missing_meta) > 0) {
    warning("Some samples in OTU table not found in metadata: ", paste(missing_meta[1:min(5,length(missing_meta))], collapse = ", "),
            if (length(missing_meta) > 5) paste0(" (+", length(missing_meta)-5, " more)"))
    # proceed anyway; metadata subset will be used for grouping only when available
  }
} else {
  # attempt to assign sample_meta rownames to SampleID
  if (!is.null(rownames(sample_meta))) {
    sample_meta <- sample_meta %>% rownames_to_column("SampleID")
  } else {
    # create SampleID from OTU table order (dangerous)
    sample_meta$SampleID <- sample_names_in_otu
  }
}

# Trim OTU to samples present
common_samples <- intersect(sample_names_in_otu, sample_meta$SampleID)
if (length(common_samples) == 0) stop("No overlapping sample names between OTU table and metadata.")
if (length(common_samples) < ncol(otu)) {
  message("Subsetting OTU table to samples with metadata (", length(common_samples), " samples).")
  otu <- otu[, common_samples, drop = FALSE]
}
# Ensure taxa have names
if (is.null(rownames(otu)) || any(rownames(otu) == "")) {
  rownames(otu) <- paste0("Taxon", seq_len(nrow(otu)))
}

# ----------------------------
# Convert to relative abundance
# ----------------------------
# If otu is integer counts, convert each sample column to counts / sum(counts)
otu_numeric <- apply(otu, 2, as.numeric)
rownames(otu_numeric) <- rownames(otu)
colnames(otu_numeric) <- colnames(otu)
# Replace NA with 0
otu_numeric[is.na(otu_numeric)] <- 0

# sample sums
sample_sums <- colSums(otu_numeric, na.rm = TRUE)
if (any(sample_sums == 0)) {
  warning("Some samples have zero total counts; their relative abundances will be NA. Consider removing these samples.")
}
relab <- sweep(otu_numeric, 2, sample_sums, FUN = "/")  # taxa x samples matrix of relative abundance

# ----------------------------
# Compute prevalence & mean abundance overall
# ----------------------------
prev_abund_tbl <- compute_prevalence_abundance(relab, detection_thresh = detection_threshold)

# Mark core taxa using thresholds
prev_abund_tbl <- prev_abund_tbl %>%
  mutate(
    IsCore = (Prevalence >= prevalence_threshold) & (MeanRelAbundance >= abundance_threshold)
  ) %>%
  arrange(desc(IsCore), desc(Prevalence), desc(MeanRelAbundance))

# Save results
write_csv(prev_abund_tbl, file.path(out_dir, "core_prevalence_abundance_overall.csv"))
message("Saved overall prevalence/abundance table to ", file.path(out_dir, "core_prevalence_abundance_overall.csv"))

# ----------------------------
# Identify core taxa list
# ----------------------------
core_taxa_overall <- prev_abund_tbl %>% filter(IsCore) %>% pull(Taxon)
message("Number of core taxa (overall): ", length(core_taxa_overall))

# Save core taxa list
write_lines(core_taxa_overall, file.path(out_dir, "core_taxa_overall_HostSpecies_Bacteria.txt"))


# ----------------------------
# General Plots to investigate trends
# ----------------------------
# 1) Prevalence vs mean abundance scatter, core highlighted
plot_df <- prev_abund_tbl %>% mutate(LogMean = log10(MeanRelAbundance + 1e-12))

p1 <- ggplot(plot_df, aes(x = Prevalence, y = LogMean, label = Taxon)) +
  geom_point(aes(color = IsCore), alpha = 0.7) +
  labs(x = "Prevalence (fraction of samples)",
       y = "log10(Mean relative abundance)",
       title = "Prevalence vs mean relative abundance",
       subtitle = paste0("Core defined as prevalence >= ", prevalence_threshold,
                         " and mean abundance >= ", abundance_threshold)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey60")) +
  theme_minimal()

ggsave(filename = file.path(out_dir, "prevalence_vs_abundance.png"), plot = p1, width = 7, height = 5, dpi = 300)
message("Saved plot: prevalence_vs_abundance.png")

# 2) Heatmap of core taxa across samples (only if any cores)
if (length(core_taxa_overall) > 0) {
  # subset relab to core taxa and order samples by group if available
  mat_core <- relab[core_taxa_overall, , drop = FALSE]
  # convert to percent for display
  mat_core_pct <- mat_core * 100
  # optional sample annotation
  sample_ann <- NULL
  if ("SampleID" %in% colnames(sample_meta)) {
    sample_meta_sub <- sample_meta %>% filter(SampleID %in% colnames(mat_core)) %>%
      column_to_rownames("SampleID")
    if (!is.null(grouping_variable) && grouping_variable %in% colnames(sample_meta_sub)) {
      sample_ann <- data.frame(Group = sample_meta_sub[[grouping_variable]])
      rownames(sample_ann) <- rownames(sample_meta_sub)
    }
  }
  # pheatmap
  pheatmap(mat_core_pct,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           annotation_col = sample_ann,
           main = "Heatmap of core taxa (relative abundance, %)",
           filename = file.path(out_dir, "heatmap_core_taxa.png"),
           width = 9, height = 6)
  message("Saved heatmap: heatmap_core_taxa.png")
} else {
  message("No core taxa found with the current thresholds; skipping heatmap.")
}

# ----------------------------
# Summary output
# ----------------------------
summary_out <- list(
  n_samples = ncol(relab),
  n_taxa = nrow(relab),
  prevalence_threshold = prevalence_threshold,
  detection_threshold = detection_threshold,
  abundance_threshold = abundance_threshold,
  n_core_taxa_overall = length(core_taxa_overall)
)
write_lines(capture.output(print(summary_out)), file.path(out_dir, "summary_MQ.txt"))
message("Saved summary to ", file.path(out_dir, "summary-final.txt"))

# ----------------------------
# OPTIONAL: Using microbiome::core (if installed & phyloseq provided)
# ----------------------------
if (exists("ps") && inherits(ps, "phyloseq")) {
  # microbiome::core can compute core across prevalence cutoffs easily
  if (requireNamespace("microbiome", quietly = TRUE)) {
    message("Also computing microbiome::core standard output (prevalence sweep)...")
    core_sweep <- microbiome::core(ps, detection = detection_threshold, prevalences = seq(0,1,by=0.01))
    # core_sweep is a list: you can inspect how many taxa remain core at different prevalence cutoffs
    # Save a simple table: prevalence cutoff -> n taxa core
    core_counts <- sapply(core_sweep, function(x) sum(x))
    core_counts_tbl <- tibble(Prevalence = seq(0,1,by=0.01), nCore = as.integer(core_counts))
    write_csv(core_counts_tbl, file.path(out_dir, "microbiome_core_sweep_counts.csv"))
    message("Saved microbiome::core sweep counts to microbiome_core_sweep_counts.csv")
  }
}

message("Core detection complete. Results are in folder: ", normalizePath(out_dir))

# ----------------------------
# Run Anlysis: Identify core taxa per group (e.g., per Site or SiteType)
# ----------------------------

grouping_variable <- "SiteType"  # must be column in metadata, or NULL

# ----------------------------
# Ensure sample_meta is a plain data.frame
# ----------------------------
if (inherits(sample_meta, "sample_data")) {
  sample_meta <- as(sample_meta, "data.frame")
}
sample_meta <- as.data.frame(sample_meta)  # in case it's still a tibble
if (!"SampleID" %in% colnames(sample_meta)) {
  sample_meta$SampleID <- rownames(sample_meta)
}

# ----------------------------
# Compute core taxa per group
# ----------------------------
core_by_group_tbl <- NULL

if (!is.null(grouping_variable) && grouping_variable %in% colnames(sample_meta)) {
  message("Computing core taxa per group variable: ", grouping_variable)
  
  # Ensure grouping variable is treated as factor/character
  sample_meta[[grouping_variable]] <- as.character(sample_meta[[grouping_variable]])
  groups <- unique(sample_meta[[grouping_variable]])
  
  core_by_group_list <- list()
  
  for (g in groups) {
    # Subset samples in this group
    samps_g <- sample_meta %>%
      filter(.data[[grouping_variable]] == g) %>%
      pull(SampleID)
    
    # Keep only samples present in the relative abundance matrix
    samps_g <- intersect(samps_g, colnames(relab))
    
    if (length(samps_g) < 2) {
      warning("Group ", g, " has <2 samples (", length(samps_g), "). Skipping group-level core computation.")
      next
    }
    
    # Subset relative abundance matrix for this group
    relab_g <- relab[, samps_g, drop = FALSE]
    
    # Compute prevalence & mean relative abundance
    pf <- compute_prevalence_abundance(relab_g, detection_thresh = detection_threshold) %>%
      mutate(
        Group = g,
        IsCore = (Prevalence >= prevalence_threshold) & (MeanRelAbundance >= abundance_threshold)
      )
    
    core_by_group_list[[as.character(g)]] <- pf
  }
  
  # Combine all groups into one table
  if (length(core_by_group_list) > 0) {
    core_by_group_tbl <- bind_rows(core_by_group_list) %>%
      select(Group, everything())
    
    # Export CSV
    out_file <- file.path(out_dir, paste0("core_prevalence_abundance_by_", grouping_variable, ".csv"))
    write_csv(core_by_group_tbl, out_file)
    message("Saved group-level core table to ", out_file)
  } else {
    message("No groups had enough samples to compute core taxa.")
  }
}

# ------------------------Heatmap using core_by_group file-----------------------
library(dplyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)

# ----------------------------
# Ensure sample_meta is plain data.frame without rownames
# ----------------------------
if (inherits(sample_meta, "sample_data")) {
  sample_meta <- as(sample_meta, "data.frame")
} else {
  sample_meta <- as.data.frame(sample_meta)
}
rownames(sample_meta) <- NULL  # remove any existing rownames

# Ensure SampleID exists
if (!"SampleID" %in% colnames(sample_meta)) {
  sample_meta$SampleID <- rownames(sample_meta)
}

# ----------------------------
# Subset core taxa
# ----------------------------
core_taxa_group <- core_by_group_tbl %>% 
  filter(IsCore) %>% 
  pull(Taxon) %>% unique()

core_relab <- relab[core_taxa_group, , drop = FALSE]

# ----------------------------
# Aggregate relative abundance by Site
# ----------------------------
site_map <- sample_meta %>%
  select(SampleID, Site, SiteType) %>%
  column_to_rownames("SampleID")

core_relab_site <- sapply(unique(site_map$Site), function(s) {
  cols <- rownames(site_map)[site_map$Site == s]
  if (length(cols) == 1) {
    core_relab[, cols]
  } else {
    rowMeans(core_relab[, cols, drop = FALSE], na.rm = TRUE)
  }
})
core_relab_site <- as.matrix(core_relab_site)

# Convert matrix to data frame with rownames as first column
core_relab_site_df <- as.data.frame(core_relab_site)
core_relab_site_df <- tibble::rownames_to_column(core_relab_site_df, var = "Taxon")

# Save as CSV
write.csv(core_relab_site_df,
          file = file.path(out_dir, "core_relab_site_MQ_Bacteria.csv"),
          row.names = FALSE)

# Read exported .csv and add taxon names
# Read in the csv file
core_relab <- read.csv("core_relab_site_TT_Fungi.csv", 
                       header = TRUE, 
                       row.names = 1)  # row.names=1 if first column is IDs

# Convert to matrix
core_relab_site <- as.matrix(core_relab)

# Check structure
str(core_relab_site)


# ----------------------------
# Prepare Site-level annotation
# ----------------------------
site_info <- sample_meta %>%
  select(Site, SiteType) %>%
  distinct()  # one row per Site

# Ensure SiteType is factor with correct order
site_info$SiteType <- factor(site_info$SiteType, levels = c("Restored", "Natural"))

# Order Sites by SiteType
site_info <- site_info[order(site_info$SiteType), ]
core_relab_site <- core_relab_site[, site_info$Site, drop = FALSE]

# Build annotation_col
annotation_col <- data.frame(
  Site = colnames(core_relab_site),
  SiteType = site_info$SiteType,
  row.names = colnames(core_relab_site)
)

# Compute gap column index between SiteTypes
site_type_numeric <- as.numeric(annotation_col$SiteType)
gap_col_index <- which(diff(site_type_numeric) != 0)

# Annotation colors
site_type_colors <- c(Restored = "#66c2a5", Natural = "#fc8d62")
site_levels <- annotation_col$Site
n_colors <- length(site_levels)
site_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_colors)
names(site_colors) <- site_levels

ann_colors <- list(
  SiteType = site_type_colors,
  Site = site_colors
)

# Heatmap color palette
abundance_colors <- colorRampPalette(c("white","magenta", "#F078D2", "red"))(100)

# ----------------------------
# Export heatmap as SVG
# ----------------------------

svg(file.path(out_dir, "core_taxa_heatmap_by_site.svg"), width = 12, height = 6)

TT_FUN<- pheatmap(
  core_relab_site,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cellwidth = 15,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  gaps_col = gap_col_index,
  fontsize_row = 8,
  color = abundance_colors,
  border_color = "grey60",
  main = "Core Taxa Heatmap (Site-level mean abundance, Restored to Natural)"
)

dev.off()

# Create combined plots
library(gridExtra)

# Extract the gtable objects from each existing pheatmap, where MQ, MS, TT represent each of the host species
g1 <- MQ_BAC$gtable
g2 <- MQ_FUN$gtable
g3 <- MS_BAC$gtable
g4 <- MS_FUN$gtable
g5 <- TT_BAC$gtable
g6 <- TT_FUN$gtable

svg("combined_heatmaps.svg", width = 20, height = 20)

# Arrange them in a single output (e.g., side by side)
final_plot <- grid.arrange(g1, g2, g3, g4, g5, g6, nrow = 3)

dev.off()



##------------------------------------------------------------------------###