library(ggtree)
library(ape)

# Filenames
TREE <- "/home/eliott.tempez/Documents/archaea_data/whole_tree.nwk"
FOCAL_SPECIES <- "GCA_001433455@Thermococcus_barophilus_CH5"

# Read the tree
tree <- read.tree(TREE)

# Plot the tree
ggtree(tree, layout = "circular", branch.length = "none") +
  geom_tippoint(aes(subset = (label == FOCAL_SPECIES)), color = "red", size = 3)
