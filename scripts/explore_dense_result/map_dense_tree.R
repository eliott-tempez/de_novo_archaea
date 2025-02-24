library(ggtree)
library(ape)

TRG_FILE <- "/home/eliott.tempez/Documents/archaea_data/dense/GCA_001433455@Thermococcus_barophilus_CH5/TRG_match_matrix.tsv"
GENE <- "TBCH5v1_RS00965"
TREE <- "/home/eliott.tempez/Documents/archaea_data/whole_tree.nwk"

gene_name <- paste0(GENE, ".mRNA.0_elongated")
# Read the TRG file
trg <- read.table(TRG_FILE, header = TRUE, sep = "\t", row.names = 1)
# Keep the row corresponding to the gene of interest
trg <- trg[gene_name, ]
# Get a different color for each cell type
colors <- sapply(trg, function(x) {
  if (x == "CDS") {
    "green"
  } else if (grepl("gS", x)) {
    "purple"
  } else if (grepl("gNS", x)) {
    "blue"
  } else {
    "black"
  }
})

# Read the tree
tree <- read.tree(TREE)
# Plot the tree
p <- ggtree(tree) + geom_tiplab() + theme_tree2() + geom_point2(aes(color = colors), size = 2)
print(p)
