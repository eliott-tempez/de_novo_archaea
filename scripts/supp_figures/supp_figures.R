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


config <- yaml.load_file("/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/my_functions/filepaths.yaml")
paths <- config$local_paths


# Filenames
tree_file <- paths$tree_file
output_dir <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/supp_figures/"


#### Tree with branch lengths ####
# Outlier
outlier <- c("GCA_001484685@Thermococcus_sp_2319x1")

tree <- read.tree(tree_file)
ggtree(tree, layout = "circular") +
  geom_tippoint(aes(subset = (label %in% outlier)), shape = 21,
                color = "black", fill = "#1273e2",
                size = 2, stroke = 1, alpha = 1,
                position = position_nudge(x = 0.015, y = 0.015))
ggsave(paste0(output_dir, "tree_branch_lengths.png"))

