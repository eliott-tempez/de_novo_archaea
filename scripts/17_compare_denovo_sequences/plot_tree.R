library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
library(yaml)
ggsave <- function(..., bg = "white",
                   width = 1000, height = 1000,
                   units = "px", dpi = 100) {
  ggplot2::ggsave(..., bg = bg,
                  width = width,
                  height = height,
                  units = units,
                  dpi = dpi)
}



# Files
config <- yaml.load_file("/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/my_functions/filepaths.yaml")
paths <- config$local_paths

tree_file <- paths$tree_file
input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/17_compare_denovo_sequences/sequence_features_good_candidates_all.csv"
good_candidates_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/good_candidates.txt"
out_folder <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/17_compare_denovo_sequences/"


# Read data
data <- read.table(input_file, header = TRUE, sep = "\t")
denovo <- scan(good_candidates_file, what = "", sep = "\n")
# Remove single aa use
data <- data[, !grepl("^[A-Z]_use$", colnames(data))]
# Remove iorfs
data <- data[!grepl("iorf", data$cds), ]
# Calculate number of denovo per genome
data_ndenovo <- data[data$cds %in% denovo, ] %>%
  group_by(genome) %>%
  summarise(n_denovo = n(), .groups = "drop")
# Ensure all genomes are represented, even those with no denovo
data_ndenovo <- data.frame(genome = unique(data$genome)) %>%
  left_join(data_ndenovo, by = "genome") %>%
  mutate(n_denovo = replace_na(n_denovo, 0))
# Pivot to longer format
descriptors <- setdiff(colnames(data), c("genome", "cds", "type"))
data <- pivot_longer(data, cols = all_of(descriptors), names_to = "feature", values_to = "value")
# Calculate median per genome for each descriptor
data_median <- data %>%
  group_by(genome, feature) %>%
  summarise(median = median(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = feature, values_from = median) %>%
  as.data.frame()
# Add number of denovo to data_median
data_median <- merge(data_median, data_ndenovo, by = "genome", all.x = TRUE)
# Add genomes as rownames
rownames(data_median) <- data_median$genome
data_median <- data_median[, -1]
# Apply signif to all columns except n_denovo
data_median[, -which(colnames(data_median) %in% c("n_denovo"))] <- 
  apply(data_median[, -which(colnames(data_median) %in% c("n_denovo"))], 2, signif, digits = 3)




# Plot tree
tree <- read.tree(tree_file)
p <- ggtree(tree, layout = "circular", branch.length = "none")

#### Basic tree
# Number of denovo
p <- p + new_scale_fill()
p <- gheatmap(p, data_median[, "n_denovo", drop = FALSE],
              width = 0.05, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "black",
                      name = "# denovo", guide = guide_colorbar(order = 1),
                      breaks = range(data_median[, "n_denovo", drop = TRUE], na.rm = TRUE))
# GC %
p <- p + new_scale_fill()
p <- gheatmap(p, data_median[, "gc_species", drop = FALSE],
              width = 0.05, colnames = FALSE, offset = 1) +
  scale_fill_gradient(low = "#009E73", high = "#6b00b2",
                      name = "GC %", guide = guide_colorbar(order = 2),
                      breaks = range(data_median[, "gc_species", drop = TRUE], na.rm = TRUE))


color_gradients <- c("#0072B2", "#CC79A7", "#D55E00", "#56B4E9", "#F0E442", "#E69F00", "#009E73")


#### Basic descriptors
basic_descriptors <- c("length", "gc_rate", "inter_gc_rate")
key_heigth <- 1 / (length(basic_descriptors))
i <- 0
q <- p
for (descriptor in basic_descriptors) {
  i <- i + 1
  q <- q + new_scale_fill()
  q <- gheatmap(q, data_median[, descriptor, drop = FALSE],
        width = 0.05, colnames = FALSE, offset = i + 3) +
  scale_fill_gradient(low = "white", high = color_gradients[i],
            name = descriptor, guide = guide_colorbar(order = i + 2),
            breaks = range(data_median[, descriptor, drop = TRUE], na.rm = TRUE))
}

# Reduce legend size
q <- q + theme(
  legend.box = "vertical",
  legend.position = "right",
  legend.key.height = unit(key_heigth, "cm"),
)
q
ggsave(file.path(out_folder, "tree_basic_descriptors.png"))


### HCA, disorder, aggregation
orfmine_descriptors <- c("hca", "disord", "aggreg")
key_heigth <- 1 / (length(orfmine_descriptors))
i <- 0
q <- p
for (descriptor in orfmine_descriptors) {
  i <- i + 1
  q <- q + new_scale_fill()
  q <- gheatmap(q, data_median[, descriptor, drop = FALSE],
                width = 0.05, colnames = FALSE, offset = i + 3) +
    scale_fill_gradient(low = "white", high = color_gradients[i],
                        name = descriptor, guide = guide_colorbar(order = i + 2),
                        breaks = range(data_median[, descriptor, drop = TRUE], na.rm = TRUE))
}
# Reduce legend size
q <- q + theme(
  legend.box = "vertical",
  legend.position = "right",
  legend.key.height = unit(key_heigth, "cm")
)
q
ggsave(file.path(out_folder, "tree_orfmine_descriptors.png"))


### biopython descriptors
biopython_descriptors <- c("aromaticity", "instability", "mean_flexibility", "hydropathy")
key_heigth <- 1 / (length(biopython_descriptors))
i <- 0
q <- p
for (descriptor in biopython_descriptors) {
  i <- i + 1
  q <- q + new_scale_fill()
  q <- gheatmap(q, data_median[, descriptor, drop = FALSE],
        width = 0.05, colnames = FALSE, offset = i + 3) +
  scale_fill_gradient(low = "white", high = color_gradients[i],
            name = descriptor, guide = guide_colorbar(order = i + 2),
            breaks = range(data_median[, descriptor, drop = TRUE], na.rm = TRUE))
}
# Reduce legend size
q <- q + theme(
  legend.box = "vertical",
  legend.position = "right",
  legend.key.height = unit(key_heigth, "cm")
)
q
ggsave(file.path(out_folder, "tree_biopython_descriptors.png"))


### AA use
aa_use_descriptors <- c("polar_use", "hydrophobic_use", "positive_use", "negative_use", "proline.glycine_use", "alanine_use", "cysteine_use")
key_heigth <- 1 / (length(aa_use_descriptors))
i <- 0
q <- p
for (descriptor in aa_use_descriptors) {
  i <- i + 1
  q <- q + new_scale_fill()
  q <- gheatmap(q, data_median[, descriptor, drop = FALSE],
        width = 0.05, colnames = FALSE, offset = i + 3) +
  scale_fill_gradient(low = "white", high = color_gradients[i],
            name = descriptor, guide = guide_colorbar(order = i + 2),
            breaks = range(data_median[, descriptor, drop = TRUE], na.rm = TRUE))
}
# Reduce legend size
q <- q + theme(
  legend.box = "vertical",
  legend.position = "right",
  legend.key.height = unit(key_heigth, "cm")
)
q
ggsave(file.path(out_folder, "tree_aa_descriptors.png"))
