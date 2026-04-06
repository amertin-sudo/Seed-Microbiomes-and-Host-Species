############################################################
# Assign fungal functional traits using microeco + FungalTraits
# Bacterial traits are assigned manually from "https://plabase.cs.uni-tuebingen.de/pb/plaba_db.php"
##

library(phyloseq)
library(microeco)

#-----------------------------
# 1. Extract data from phyloseq object
#-----------------------------

# OTU table
otu_table_fungi <- as.data.frame(otu_table(fungi))

# Ensure taxa are rows
if (taxa_are_rows(fungi) == FALSE) {
  otu_table_fungi <- t(otu_table_fungi)
}

# Taxonomy table
taxonomy_table_fungi <- as.data.frame(tax_table(fungi))

# Sample metadata
sample_info_fungi <- as.data.frame(sample_data(fungi))

#-----------------------------
# 2. Create microeco object
#-----------------------------

dataset_fungi <- microtable$new(
  sample_table = sample_info_fungi,
  otu_table = otu_table_fungi,
  tax_table = taxonomy_table_fungi
)

#-----------------------------
# 3. Assign fungal traits
#-----------------------------

t2_fungi <- trans_func$new(dataset_fungi)

t2_fungi$cal_spe_func(fungi_database = "FungalTraits")

#-----------------------------
# 4. Export results
#-----------------------------

write.csv(t2_fungi$res_spe_func,
          file = "fungi_functional_traits.csv",
          row.names = FALSE)

# Select indicator ASVs and add trait info to dataframe used to construct horizontal barchart in previous script.
