library(dplyr)



input_file = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/integrity_analysis/nc_coverage.tsv"
out_folder = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/"
outgroup_coverage_threshold <- 80


# Read the data
data <- read.table(input_file, header = TRUE, sep = "\t")


# Start with last outgroup
data_last <- data[data$outgroup == "last", ]
# Get the % of genes for which either qcov or qcov_orf <= 70
data_last <- data_last %>%
  group_by(denovo) %>%
  summarise(
    total_genes = n(),
    low_qcov_genes = sum(qcov <= 70 | qcov_orf <= 70),
    low_qcov_percentage = (low_qcov_genes / total_genes) * 100
  )
# Write the results to a file
write.table(data_last, file = paste0(out_folder, "integrity_analysis/low_qcov_last_outgroup.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


# Same for the before last outgroup
data_before_last <- data[data$outgroup == "before_last", ]
data_before_last <- data_before_last %>%
  group_by(denovo) %>%
  summarise(
    total_genes = n(),
    low_qcov_genes = sum(qcov <= 70 | qcov_orf <= 70),
    low_qcov_percentage = (low_qcov_genes / total_genes) * 100
  )
write.table(data_before_last, file = paste0(out_folder, "integrity_analysis/low_qcov_before_last_outgroup.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)



# Get the genes with percentage of no integrity >= threshold
genes_last_below_threshold <- data_last[data_last$low_qcov_percentage > outgroup_coverage_threshold, "denovo"]
genes_before_last_below_threshold <- data_before_last[data_before_last$low_qcov_percentage > outgroup_coverage_threshold, "denovo"]
# Get the intersection of the two sets of genes
good_candidates <- intersect(genes_last_below_threshold, genes_before_last_below_threshold)
# Write the good candidates to a file
write.table(good_candidates, file = paste0(out_folder, "good_candidates.txt"), sep = "", row.names = FALSE, col.names = FALSE, quote = FALSE)
