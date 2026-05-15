
#Importing data and converting to a phyloseq object for further analyses.
#Bacterial 16S rRNA data and ITS2 data were imported and analysed separately.

#set working directory
setwd("/pathtofolder/Inputdatafolder")

# load useful function
KillZeroRCs <- function(x) {
  x[ which( rowSums(x) != 0) , ] -> x
  x[ , which( colSums(x) != 0) ] -> x
  return(x)
}

# read in OTU table
ASV <- read.table(
  "ASV_16S_All_Hosts_Subset_Final.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

# remove zero-sum row/columns
ASV <- KillZeroRCs(ASV)

# read in tax table
# fill logical indicates that rows have unequal lengths due to blank fields
TAXA <- read.table(
  "Taxonomy_Filtered_16S_Revised.txt",
  sep = "\t",
  fill = TRUE,
  row.names = 1
)

# add levels of taxonomy to tax table
colnames(TAXA) <- c("Domain","Phylum", "Class", "Order", "Family", "Genus")

# read in METAdata
META <- read.table(
  "Metadata_Run2R3R4Final_v2_16S.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1
)


# convert OTU and tax tables to matrices for phyloseq
ASV <- as.matrix(ASV)
TAXA <- as.matrix(TAXA)

# install bioconductor to install BiocManager then install the packages phyloseq and Microbiome
#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("phyloseq", force = TRUE)
# had to install 'RCurl' using the source binaries at the Github page of RCurl.

#if (!require("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("microbiome")
#BiocManager::install("Phylosmith")

# load libraries
library(ape)
library(phyloseq)
library(microbiome)
library(phylosmith)

# combine OTU, taxa and METAdata files into a phyloseq object
phy <- phyloseq(
  otu_table(ASV,taxa_are_rows = T),
  tax_table(TAXA),
  sample_data(META)
)

# remove zero sum ASVs
phy <- prune_taxa((taxa_sums(phy) > 0), phy)


# remove non-bacterial reads. Repeat the same for fungal dataset.
bacteria <- subset_taxa(phy, Domain == "Bacteria") 

# remove zero sum ASV
bacteria <- prune_taxa((taxa_sums(bacteria) > 0),bacteria) 
