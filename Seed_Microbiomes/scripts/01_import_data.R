
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
  "asv.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

# remove zero-sum row/columns
ASV <- KillZeroRCs(ASV)

# read in tax table
# fill logical indicates that rows have unequal lengths due to blank fields
TAXA <- read.table(
  "taxonomy.txt",
  sep = "\t",
  fill = TRUE,
  row.names = 1
)

# add levels of taxonomy to tax table
colnames(TAXA) <- c("Domain","Phylum", "Class", "Order", "Family", "Genus")

# read in METAdata
META <- read.table(
  "Metadata.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1
)


# convert OTU and tax tables to matrices for phyloseq
ASV <- as.matrix(ASV)
TAXA <- as.matrix(TAXA)

# install bioconductor to install BiocManager then install the packages phyloseq and Microbiome
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq", force = TRUE)
# had to install 'RCurl' using the source binaries at the Github page of RCurl.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")
BiocManager::install("Phylosmith")

# load libraries
library(ape)
library(phyloseq)
library(microbiome)

# combine OTU, taxa and METAdata files into a phyloseq object
phy <- phyloseq(
  otu_table(ASV,taxa_are_rows = T),
  tax_table(TAXA),
  sample_data(META)
)

# remove zero sum ASV
phy <- prune_taxa((taxa_sums(phy) > 0), phy) #2302 ASVs in 156 samples


# remove non-bacterial reads. Repeat the same for fungal dataset.
bacteria <- subset_taxa(phy, Domain == "Bacteria") 

# remove zero sum ASV
bacteria <- prune_taxa((taxa_sums(bacteria) > 0),bacteria) 
