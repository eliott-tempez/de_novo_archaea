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


input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/denovo_noncoding_status.tsv"
out_folder <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match"
#input_file <- "/home/eliott/Documents/UNI/M2/Stage/M2_stage_I2BC/results/14_get_noncoding_match/denovo_noncoding_status.tsv"



# Read the data
data <- read.table(input_file, header = TRUE, sep = "\t")
# Remove columns that have 3 codons or less
data[, c("intergenic", "f0", "f1", "f2")] <- lapply(data[, c("intergenic", "f0", "f1", "f2")], function(x) ifelse(x < 10, 0, x))




##################################
##### Analyse all NC matches #####
##################################
# Number of different outgroups for each gene
outgroups <- data %>%
  group_by(gene_id) %>%
  summarise(outgroup_nb = n_distinct(outgroup_nb)) %>%
  ungroup()
# Violin plot
ggplot(outgroups, aes(x = "", y = outgroup_nb)) +
  geom_violin(fill = "#ff7070", alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "#a7a7a7", alpha = 0.8) +
  labs(title = "Number of different outgroups with NC match\nfor each de novo gene") +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(breaks = seq(0, max(outgroups$outgroup_nb), 2),
                     limits = c(0, max(outgroups$outgroup_nb)))
#ggsave(file.path(out_folder, "outgroup_nb_violin_plot.png"))


# Print outgroup summary to file
outgroups_count <- data %>%
  group_by(gene_id, outgroup_nb) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = outgroup_nb, values_from = count, values_fill = 0)
outgroup_count <- outgroups_count[c("gene_id", as.character(seq(1, 15)))]
#write.table(outgroup_count, file = file.path(out_folder, "outgroup_count.tsv"),
#            sep = "\t", row.names = FALSE, quote = FALSE)

















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
different_origin <- data_origin %>%
  rowwise() %>%
  filter(
    !(intergenic %in% c(0, n)) |
      !(f0 %in% c(0, n)) |
      !(f1 %in% c(0, n)) |
      !(f2 %in% c(0, n))
  ) %>%
  ungroup()
print(paste("there are", nrow(different_origin), "genes with different origins across neighbours."))

# Genes with the same origin accross neighbours
same_origin <- data_origin %>%
  rowwise() %>%
  filter(
    intergenic %in% c(0, n) &
      f0 %in% c(0, n) &
      f1 %in% c(0, n) &
      f2 %in% c(0, n)
  ) %>%
  ungroup()
print(paste("there are", nrow(same_origin), "genes with the same origin across all neighbours."))
















##################################
####### Furthest outgroup ########
##################################
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
  ungroup()

# Origin of the furthest outgroup
furthest_data_origin <- furthest_outgroup_data %>%
  group_by(gene_id) %>%
  summarise(
    intergenic = sum(intergenic > 0),
    f0 = sum(f0 > 0),
    f1 = sum(f1 > 0),
    f2 = sum(f2 > 0),
    n = n()
  )
furthest_different_origin <- furthest_data_origin %>%
  rowwise() %>%
  filter(
    !(intergenic %in% c(0, n)) |
      !(f0 %in% c(0, n)) |
      !(f1 %in% c(0, n)) |
      !(f2 %in% c(0, n))
  ) %>%
  ungroup()







# For the rest of the analysis, keep only the first row for each unique gene
furthest_outgroup_data <- furthest_outgroup_data %>%
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
#ggsave(file.path(out_folder, "furthest_outgroup_venn_diagram.png"))



## Get the genes with the NC match straddling a gene
straddling_genes <- furthest_outgroup_data %>%
  filter(f0 > 0 | f1 > 0 | f2 > 0)
write.table(straddling_genes, file = file.path(out_folder, "straddling_genes.csv"),
            row.names = FALSE, quote = FALSE, sep = "\t")
