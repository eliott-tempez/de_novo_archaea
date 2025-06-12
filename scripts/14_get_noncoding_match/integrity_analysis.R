library(dplyr)



input_file = "/home/eliott/Documents/UNI/M2/Stage/M2_stage_I2BC/results/14_get_noncoding_match/nc_coverage.tsv"


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


# Before last outgroup
data_before_last <- data[data$outgroup == "before_last", ]
data_before_last <- data_before_last %>%
  group_by(denovo) %>%
  summarise(
    total_genes = n(),
    low_qcov_genes = sum(qcov <= 70 | qcov_orf <= 70),
    low_qcov_percentage = (low_qcov_genes / total_genes) * 100
  )
