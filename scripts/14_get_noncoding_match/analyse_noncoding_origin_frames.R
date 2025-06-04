library(dplyr)
library(ggVennDiagram)
library(ggplot2)

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



##### Origin of the furthest outgroup #####
# Compute furthest_outgroup for each gene_id
furthest_outgroup_df <- data %>%
  group_by(gene_id) %>%
  summarise(furthest_outgroup = max(outgroup_nb, na.rm = TRUE)) %>%
  ungroup()
# Add furthest_outgroup column to data
data <- data %>%
  left_join(furthest_outgroup_df, by = "gene_id")
# Keep only the furthest outgroup for each gene_id
furthest_outgroup_data <- data %>%
  rowwise() %>%
  filter(outgroup_nb == furthest_outgroup) %>%
  ungroup() %>%
  # Keep only 1 row per gene_id
  distinct(gene_id, .keep_all = TRUE)


## Venn diagram
intergenic_genes <- furthest_outgroup_data %>%
  filter(intergenic > 0) %>%
  pull(gene_id)
f0_genes <- furthest_outgroup_data %>%
  filter(f0 > 0) %>%
  pull(gene_id)
f1_genes <- furthest_outgroup_data %>%
  filter(f1 > 0) %>%
  pull(gene_id)
f2_genes <- furthest_outgroup_data %>%
  filter(f2 > 0) %>%
  pull(gene_id)
venn_list = list(
  intergenic = intergenic_genes,
  f0 = f0_genes,
  f1 = f1_genes,
  f2 = f2_genes
)
# Plot
ggVennDiagram(venn_list, label_alpha = 0) +
  labs(title = "Origin of the furthest outgroup for de novo genes in NC match") +
  theme(legend.position = "none") +
    scale_fill_gradientn(
    colors = c("white", "#ffbec8", "#ff7070"),
    values = scales::rescale(c(0, 7, 54)))










##### Analyse all NC matches #####
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
  rowwise() %>%
  filter(
    !(intergenic %in% c(0, n)) |
      !(f0 %in% c(0, n)) |
      !(f1 %in% c(0, n)) |
      !(f2 %in% c(0, n))
  ) %>%
  ungroup()
print(different_origin)

same_origin <- data_origin %>%
  rowwise() %>%
  filter(
    intergenic %in% c(0, n) &
      f0 %in% c(0, n) &
      f1 %in% c(0, n) &
      f2 %in% c(0, n)
  ) %>%
  ungroup()

