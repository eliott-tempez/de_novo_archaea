library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
library(argparse)
ggsave <- function(..., bg = "white",
                   width = 1000, height = 1000,
                   units = "px", dpi = 100) {
  ggplot2::ggsave(..., bg = bg,
                  width = width,
                  height = height,
                  units = units,
                  dpi = dpi)
}

# Filenames
tree_file <- "/home/eliott.tempez/Documents/archaea_data/whole_tree.nwk"
conservation_db_dir <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/calculate_conservation/"
output_dir <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/calculate_conservation/"

# Get the focal species
parser <- ArgumentParser(description = "Map gradient tree")
parser$add_argument('--focal_species', type = 'character', required = TRUE, help = 'Focal species identifier')
args <- parser$parse_args()
focal_species <- args$focal_species


# Get the conservation db
conservation_db_file <- paste0(conservation_db_dir,
                               "conservation_df_",
                               focal_species, ".tsv")
conservation_db <- read.table(conservation_db_file,
                              header = TRUE, sep = "\t",
                              row.names = 1)
# Get the ssearch result
ssearch_db <- conservation_db[, c("ssearch_f0", "ssearch_f1", "ssearch_f2", "ssearch_f0_comp", "ssearch_f1_comp", "ssearch_f2_comp")]
colnames(ssearch_db) <- c("+0 ", "+1 ", "+2 ", "-0 ", "-1 ", "-2 ")
conservation_db <- conservation_db[, c("n_denovo", "n_trg", "n_cds", "n_intergenic")]
# Avoid division by 0 if no de novo
if (conservation_db[focal_species, "n_denovo"] == 0) {
  conservation_db[focal_species, "n_denovo"] <- 1
}
# Change to percents
default_val <- as.numeric(conservation_db[focal_species, ])
conservation_db_p <- sweep(conservation_db, 2, default_val, FUN = "/") * 100
default_val <- as.numeric(ssearch_db[focal_species, ])
ssearch_db_p <- sweep(ssearch_db, 2, default_val, FUN = "/") * 100



####### Conservation ########
# Read the tree
tree <- read.tree(tree_file)
# Plot the tree
p <- ggtree(tree, layout = "circular", branch.length = "none") +
  geom_tippoint(aes(subset = (label == focal_species)), shape = 23,
                color = "black", fill = "yellow",
                size = 3, stroke = 1, alpha = 1)

# Add the heatmaps
## De novo ##
p <- p + new_scale_fill()
p <- gheatmap(p, conservation_db_p[, "n_denovo", drop = FALSE],
              width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "black", name = "De novo",
                      guide = guide_colorbar(order = 1), limits = c(0, 100))
## TRGs ##
p <- p + new_scale_fill()
p <- gheatmap(p, conservation_db_p[, "n_trg", drop = FALSE], offset = 1,
              width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "red", name = "TRGs",
                      guide = guide_colorbar(order = 2), limits = c(0, 100))
## Noncoding ##
p <- p + new_scale_fill()
p <- gheatmap(p, conservation_db_p[, "n_intergenic", drop = FALSE], offset = 2,
              width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "#6b00b2", name = "Intergenic",
                      guide = guide_colorbar(order = 3), limits = c(0, 100))
## CDSs ##
p <- p + new_scale_fill()
p <- gheatmap(p, conservation_db_p[, "n_cds", drop = FALSE], offset = 3,
              width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "#e67b00", name = "CDS",
                      guide = guide_colorbar(order = 4), limits = c(0, 100))


# Add the title
p <- p + ggtitle(paste0("Conservation rate for\n", focal_species)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -10))

# Save the image
ggsave(paste0(output_dir, "conservation_map_", focal_species, ".png"))







######## SSEARCH ########
# Plot the tree
p <- ggtree(tree, layout = "circular", branch.length = "none") +
  geom_tippoint(aes(subset = (label == focal_species)), shape = 23,
                color = "black", fill = "yellow",
                size = 3, stroke = 1, alpha = 1)
# Add each reading frame in the same heatmap color
## Frame 0 ##
p <- p + new_scale_fill()
p <- gheatmap(p, ssearch_db_p[, "+0 ", drop = FALSE],
              width = .05, colnames = TRUE, colnames_angle = 90, font.size = 3) +
  scale_fill_gradient(low = "white", high = "#009E73", name = "conservation %",
                      guide = guide_colorbar(), limits = c(0, 100))
## Frame 1 ##
p <- p + new_scale_fill()
p <- gheatmap(p, ssearch_db_p[, "+1 ", drop = FALSE], offset = 1,
              width = .05, colnames = TRUE, colnames_angle = 90, font.size = 3) +
  scale_fill_gradient(low = "white", high = "#009E73",
                      guide = "none", limits = c(0, 100))
## Frame 2 ##
p <- p + new_scale_fill()
p <- gheatmap(p, ssearch_db_p[, "+2 ", drop = FALSE], offset = 2,
              width = .05, colnames = TRUE, colnames_angle = 90, font.size = 3) +
  scale_fill_gradient(low = "white", high = "#009E73",
                      guide = "none", limits = c(0, 100))
## Frame -0 ##
p <- p + new_scale_fill()
p <- gheatmap(p, ssearch_db_p[, "-0 ", drop = FALSE], offset = 3,
              width = .05, colnames = TRUE, colnames_angle = 90, font.size = 3) +
  scale_fill_gradient(low = "white", high = "#009E73",
                      guide = "none", limits = c(0, 100))
## Frame -1 ##
p <- p + new_scale_fill()
p <- gheatmap(p, ssearch_db_p[, "-1 ", drop = FALSE], offset = 4,
              width = .05, colnames = TRUE, colnames_angle = 90, font.size = 3) +
  scale_fill_gradient(low = "white", high = "#009E73",
                      guide = "none", limits = c(0, 100))
## Frame -2 ##
p <- p + new_scale_fill()
p <- gheatmap(p, ssearch_db_p[, "-2 ", drop = FALSE], offset = 5,
              width = .05, colnames = TRUE, colnames_angle = 90, font.size = 3) +
  scale_fill_gradient(low = "white", high = "#009E73",
                      guide = "none", limits = c(0, 100))

# Add the title
p <- p + ggtitle(paste0("Conservation rate for the 6 reading frames\n", focal_species)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -10))

# Save the image
ggsave(paste0(output_dir, "conservation_map_frames_", focal_species, ".png"))