library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)

# Filenames
tree_file <- "/home/eliott.tempez/Documents/archaea_data/whole_tree.nwk"
conservation_db_dir <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/calculate_conservation/"
focal_species <- "GCA_001433455@Thermococcus_barophilus_CH5"

# Get the conservation db
conservation_db_file <- paste0(conservation_db_dir, "conservation_df_", focal_species, ".tsv")
conservation_db <- read.table(conservation_db_file, header = TRUE, sep = "\t", row.names = 1)
conservation_db[is.na(conservation_db)] <- 0
# Change to percents
default_val <- conservation_db[focal_species, ]
conservation_db_p <- conservation_db / default_val * 100



# Read the tree
tree <- read.tree(tree_file)
# Plot the tree
p <- ggtree(tree, layout = "circular", branch.length = "none") +
  geom_tippoint(aes(subset = (label == focal_species)), shape = 23,
                color = "black", fill = "yellow",
                size = 3, stroke = 1, alpha = 1)
# Add the heatmap
p <- gheatmap(p, conservation_db_p, offset = .8,
              width = .1, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "red", name = "TRGs")
p <- p + ggtitle(paste0("Conservation rate for\n", focal_species)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -10))
p
