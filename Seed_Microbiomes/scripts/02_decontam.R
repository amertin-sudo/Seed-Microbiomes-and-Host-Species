#02_Remove PCR control and DNA extraction control samples and ASVs from the phyloseq object

#BiocManager::install("decontam")
library(decontam)
library(microbiome)
#identify contaminants first from PCR negatives
consList <- isContaminant(seqtab = bacteria, neg = "PCR_Control", method = "prevalence")

# pull out the names of contaminants
cons <- rownames(consList)[consList$contaminant=="TRUE"]
cons <- as.character(cons) #6 contaminants identified from PCR_controls
write.csv(x = consList, row.names = TRUE, file = "consList.csv")
# - - - - - - - - - - - - - - - - - - - - - - - - - #

# to get info on the contaminants, uncomment the following code to
# run it ON THE FILE WITH THE CONTAMINANT ASVs IN-SITU,
# then combine the consPer.csv and taxonomy.csv file data

# subset the non neg-control samples
subset_samples <- subset_samples(bacteria, PCR_Control == "FALSE")
# merge the samples
merged_samples <- merge_samples(subset_samples, "PCR_Control", fun = sum)
# transform counts to percentages
merged_samples <- transform_sample_counts(merged_samples, function(x) 100 * x/sum(x))
# extract the cons percentage data
pruned_samples <- prune_taxa(x = merged_samples, taxa = cons)
# write otu table to dataframe
pruned <- data.frame(t(pruned_samples@otu_table))
# write pruned to csv
write.csv(x = pruned, row.names = TRUE, file = "consPCRcontrol.csv")
# subset the contaminant ASVs
bacteriaCons <- prune_taxa(bacteria, taxa = cons)
# write the contaminants to a file for reference
contaminants <-bacteriaCons@tax_table
contaxa <- contaminants@.Data
write.csv(contaxa, "contaxaPCR_control.csv")


# - - - - - - - - - - - - - - - - - - - - - - - - - #

#remove the contaminants from the main bacteria phyloseq object
bacteria <- remove_taxa(bacteria, taxa = cons)


#Remove PCR_controls. must have column in metadata that specifies PCR_Control status for each sample
bacteria <- subset_samples(bacteria, PCR_Control == "FALSE") 


# identify contaminants from extraction blanks
consList2 <- isContaminant(seqtab = bacteria, neg = "Extract_Control", method = "prevalence")

# pull out the names of contaminants
cons2 <- rownames(consList2)[consList2$contaminant=="TRUE"]
cons2 <- as.character(cons2) #gives you count of contaminants identified from extraction blanks

# - - - - - - - - - - - - - - - - - - - - - - - - - #

# following code will provide taxonomic info for the identified contaminants
# run it ON THE FILE WITH THE CONTAMINANT ASVs IN-SITU,
# then combine the consPer.csv and taxonomy.csv file data

#subset the non neg-control samples
subset_samples2 <- subset_samples(bacteria, Extract_Control == "FALSE")
#merge the samples
merged_samples2 <- merge_samples(subset_samples2, group = "Extract_Control", fun = sum)
#transform counts to percentages
merged_samples2 <- transform_sample_counts(merged_samples2, function(x) 100 * x/sum(x))
#extract the cons percentage data
pruned_samples2 <- prune_taxa(x = merged_samples2, taxa = cons2)
#write otu table to dataframe
final2 <- data.frame(t(pruned_samples2@otu_table))
#write final2 to csv
write.csv(x = final2, row.names = TRUE, file = "consExtraction.csv")
#subset the contaminant ASVs
bacteriaCons2 <- prune_taxa(bacteria, taxa = cons2)
#write the contaminants to a file for reference
contaminants2 <-bacteriaCons2@tax_table
contaxa2 <- contaminants2@.Data
write.csv(contaxa2, "contaxaExtraction.csv")

# remove the extraction contaminants from the main bacteria phyloseq object
bacteria <- remove_taxa(bacteria, taxa = cons2)

# Remove extraction controls from bacteria 
bacteria <- subset_samples(bacteria, Extract_Control=="FALSE")

#subset mock community sample and remove them from bacteria phyloseq object
MC <- subset_samples(bacteria, Site == "Mock_Community")
MC <- prune_taxa((taxa_sums(MC)>0),MC) #9 ASVs in 1 sample
bacteria <- remove_samples(c("Mock_4_FWD24", "Mock_5_FWD24", "Mock_6_FWD24", "Mock_7_FWD16", "Mock_8_FWD16", "Mock_Fungi_1_Control_FWD24", "Mock_Fungi_2_Control_FWD24", "Mock_Fungi_3_Control_FWD24"), bacteria)
bacteria <- prune_taxa((taxa_sums(bacteria) > 0), bacteria)

# Inspect values first
unique(sample_data(bacteria)$Site)

# Convert to character in case it is stored as a factor
sample_data(bacteria)$site <- as.character(sample_data(bacteria)$Site)

# Remove Mock Community samples
bacteria <- subset_samples(
  bacteria,
  Site != "Mock_Community"
)

# Remove taxa that now have zero counts
bacteria <- prune_taxa(taxa_sums(bacteria) > 0, bacteria)

#view reads by sample
sort(sample_sums(bacteria)) 

#remove samples with less than n reads
bacteria <- prune_samples((sample_sums(bacteria) > 0), bacteria) 
sort(sample_sums(bacteria)) 
#bacteria <- prune_samples((sample_sums(bacteria) > 400), bacteria) # for this step must not have samples with 0 reads.


# - - - - - - - - - - - - - - - - - - - - - - - - - #

# Extract metadata from phyloseq object for correct sample and ASV sums per host etc. after QC steps
metadata <- data.frame(sample_data(bacteria))

# View metadata
head(metadata)

# Export as CSV
write.csv(
  metadata,
  file = "bacteria_sample_metadata_postQC.csv",
  row.names = TRUE
)

# Extract ASV table for ASV counts
asv_mat <- otu_table(bacteria)

if (!taxa_are_rows(bacteria)) {
  asv_mat <- t(asv_mat)
}

asv_table <- as.data.frame(asv_mat)

# Preview
head(asv_table)

write.csv(
  asv_table,
  file = "bacteria_ASV_table_postQC.csv",
  row.names = TRUE
)


