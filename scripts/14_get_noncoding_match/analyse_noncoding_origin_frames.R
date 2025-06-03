library(dplyr)

ggsave <- function(..., bg = "white",
                   width = 1000, height = 1000,
                   units = "px", dpi = 100) {
  ggplot2::ggsave(..., bg = bg,
                  width = width,
                  height = height,
                  units = units,
                  dpi = dpi)
}


input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/denovo_noncoding_status.tsv"


data <- read.table(input_file, header = TRUE, sep = "\t")
# Remove columns that have 3 codons or less
#data[, c("intergenic", "f.0", "f.1", "f.2")] <- lapply(data[, c("intergenic", "f.0", "f.1", "f.2")], function(x) ifelse(x < 10, 0, x))




# Group by denovo and print the origins for each nc match
data %>%
  group_by(gene_id)