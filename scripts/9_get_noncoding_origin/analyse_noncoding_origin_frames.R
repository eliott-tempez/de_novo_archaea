# This script analyses the denovo_noncoding_status.tsv file
# created by get_all_origin_frames_denovo.py.


library(dplyr)
library(tidyr)
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


input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/9_get_noncoding_origin/denovo_noncoding_status.tsv"
out_folder <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/9_get_noncoding_origin"
cluster_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/5_de_novo_clustering/good_candidates_clustering.tsv"



# Read the data
data <- read.table(input_file, header = TRUE, sep = "\t")
# Remove columns that have 3 codons or less
data[, c("intergenic", "f0", "f1", "f2")] <- lapply(data[, c("intergenic", "f0", "f1", "f2")], function(x) ifelse(x < 10, 0, x))
# Add the clusters
cluster_df <- read.table(cluster_file, header = FALSE, sep = "\t")
colnames(cluster_df) <- c("cluster", "gene_id")
data <- data %>%
  left_join(cluster_df, by = "gene_id") %>%
  mutate(cluster = ifelse(is.na(cluster), "no_cluster", cluster))



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
  group_by(gene_id, cluster) %>%
  summarise(
    intergenic = sum(intergenic > 0),
    f0 = sum(f0 > 0),
    f1 = sum(f1 > 0),
    f2 = sum(f2 > 0),
    n = n(),
    .groups = "drop"
  )

# Print the ones that are 100% intergenic
furthest_origin_intergenic <- furthest_data_origin %>%
  filter(f0 == 0 & f1 == 0 & f2 == 0) %>%
  arrange(cluster)
print(paste("there are", nrow(furthest_origin_intergenic), "genes that are 100% intergenic for the furthest outgroup."))

# Print the ones that are not 100% intergenic
furthest_origin_not_intergenic <- furthest_data_origin %>%
  filter(f0 > 0 | f1 > 0 | f2 > 0) %>%
  arrange(cluster)
print(furthest_origin_not_intergenic)
print(furthest_outgroup_data[furthest_outgroup_data$gene_id %in% furthest_origin_not_intergenic$gene_id,])
furthest_outgroup_data %>%
  filter(gene_id %in% furthest_origin_not_intergenic$gene_id) %>%
  select(gene_id, cluster, intergenic, f0, f1, f2) %>%
  arrange(cluster) %>%
  print(n = Inf)


# Get the different/same origins
furthest_different_origin <- furthest_data_origin %>%
  rowwise() %>%
  filter(
    !(intergenic %in% c(0, n)) |
      !(f0 %in% c(0, n)) |
      !(f1 %in% c(0, n)) |
      !(f2 %in% c(0, n))
  ) %>%
  ungroup()
print(paste("there are", nrow(furthest_different_origin), "genes with different origins for the furthest outgroup."))
furthest_same_origin <- furthest_data_origin %>%
  rowwise() %>%
  filter(
    intergenic %in% c(0, n) &
      f0 %in% c(0, n) &
      f1 %in% c(0, n) &
      f2 %in% c(0, n)
  ) %>%
  ungroup()
print(paste("there are", nrow(furthest_same_origin), "genes with the same origin for the furthest outgroup."))


# Genes that both have the same origin, and are not 100% intergenic
same_origin_not_intergenic <- inner_join(furthest_origin_not_intergenic, furthest_same_origin, by = "gene_id") %>%
  select(gene_id) %>%
  left_join(furthest_outgroup_data, by = "gene_id") %>%
  arrange(cluster)
# Genes that both have different origins, and are not 100% intergenic
different_origin_not_intergenic <- inner_join(furthest_origin_not_intergenic, furthest_different_origin, by = "gene_id") %>%
  select(gene_id) %>%
  left_join(furthest_outgroup_data, by = "gene_id") %>%
  arrange(cluster)
# Write to file
write.table(different_origin_not_intergenic,
            file = file.path(out_folder, "straddling_genes_discordant.csv"),
            row.names = FALSE, quote = FALSE, sep = "\t")






# For the rest of the analysis, keep only the first row for each unique gene
# that both have the same origin, and are not 100% intergenic
# and keep only first row for each unique cluster
straddling_data <- same_origin_not_intergenic %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  distinct(cluster, .keep_all = TRUE)
write.table(straddling_data,
            file = file.path(out_folder, "straddling_genes_concordant.csv"),
            row.names = FALSE, quote = FALSE, sep = "\t")


