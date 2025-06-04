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


#input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/denovo_noncoding_status.tsv"
input_file <- "/home/eliott/Documents/UNI/M2/Stage/M2_stage_I2BC/results/14_get_noncoding_match/denovo_noncoding_status.tsv"


# Read the data
data <- read.table(input_file, header = TRUE, sep = "\t")
# Remove columns that have 3 codons or less
data[, c("intergenic", "f0", "f1", "f2")] <- lapply(data[, c("intergenic", "f0", "f1", "f2")], function(x) ifelse(x < 10, 0, x))

# Extract the genes which don't have the same origin across all neighbours
data_origin <- data %>%
  group_by(gene_id) %>%
  summarise(
    intergenic = sum(intergenic > 0),
    f0 = sum(f0 > 0),
    f1 = sum(f1 > 0),
    f2 = sum(f2 > 0),
    n = n()
  )
print("Genes with different origins across neighbours:")
different_origin <- data_origin %>%
  filter(!(intergenic %in% c(0, n)) |
           !(f0 %in% c(0, n)) |
           !(f1 %in% c(0, n)) |
           !(f2 %in% c(0, n)))
same_origin <- data_origin %>%
  filter(intergenic %in% c(0, n) &
           f0 %in% c(0, n) &
           f1 %in% c(0, n) &
           f2 %in% c(0, n))
