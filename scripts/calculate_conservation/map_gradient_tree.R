library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
ggsave <- function(..., bg = "white", width = 1000, height = 1000, units = "px", dpi = 100) {
  ggplot2::ggsave(..., bg = bg, width = width, height = height, units = units, dpi = dpi)
}

# Filenames
tree_file <- "/home/eliott.tempez/Documents/archaea_data/whole_tree.nwk"
conservation_db_dir <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/calculate_conservation/"
focal_species <- "GCA_001433455@Thermococcus_barophilus_CH5"
output_dir <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/calculate_conservation/"

# Get the conservation db
conservation_db_file <- paste0(conservation_db_dir, "conservation_df_", focal_species, ".tsv")
conservation_db <- read.table(conservation_db_file, header = TRUE, sep = "\t", row.names = 1)
conservation_db[is.na(conservation_db)] <- 0
# Change to percents
default_val <- as.numeric(conservation_db[focal_species, ])
# To remove : avoiding division by 0
default_val[default_val == 0] <- 1
conservation_db_p <- sweep(conservation_db, 2, default_val, FUN = "/") * 100



# Read the tree
tree <- read.tree(tree_file)
# Plot the tree
p <- ggtree(tree, layout = "circular", branch.length = "none") +
  geom_tippoint(aes(subset = (label == focal_species)), shape = 23,
                color = "black", fill = "yellow",
                size = 3, stroke = 1, alpha = 1)

# Add the heatmaps
## TRGs ##
p <- p + new_scale_fill()
p <- gheatmap(p, conservation_db_p[, "n_trg", drop = FALSE],
        width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "red", name = "TRGs",
            guide = guide_colorbar(order = 1), limits = c(0, 100))
## Noncoding ##
p <- p + new_scale_fill()
p <- gheatmap(p, conservation_db_p[, "n_noncoding", drop = FALSE], offset = 1,
              width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "#6b00b2", name = "Noncoding",
                      guide = guide_colorbar(order = 2), limits = c(0, 100))
## Frame 1 ##
p <- p + new_scale_fill()
p <- gheatmap(p, conservation_db_p[, "n_f1", drop = FALSE], offset = 2,
              width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "#009E73", name = "+1",
                      guide = guide_colorbar(order = 3), limits = c(0, 100))
## Frame 2 ##
p <- p + new_scale_fill()
p <- gheatmap(p, conservation_db_p[, "n_f2", drop = FALSE], offset = 3,
              width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "#0072B2", name = "+2",
                      guide = guide_colorbar(order = 4), limits = c(0, 100))
## CDSs ##
p <- p + new_scale_fill()
p <- gheatmap(p, conservation_db_p[, "n_cds", drop = FALSE], offset = 4,
              width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "#e67b00", name = "CDS",
                      guide = guide_colorbar(order = 5), limits = c(0, 100))


# Add the title
p <- p + ggtitle(paste0("Conservation rate for\n", focal_species)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -10))
p

# Save the image
ggsave(paste0(output_dir, "conservation_map_", focal_species, ".png"))