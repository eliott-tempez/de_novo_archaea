library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)

trg_file <- "/home/eliott.tempez/Documents/archaea_data/dense/Thermococcus_barophilus_CH1/TRG_match_matrix.tsv"
gene_name <- "TBCH1v1_0272.mRNA.0_elongated"
tree <- "/home/eliott.tempez/Documents/archaea_data/whole_tree.nwk"
focal_species = "Thermococcus_barophilus_CH1"


# Read the TRG file
trg <- read.table(trg_file, header = TRUE, sep = "\t", row.names = 1)
# Keep the row corresponding to the gene of interest
trg <- trg[gene_name, ]
trg <- t(trg)
colnames(trg) <- "type"
rownames(trg) <- gsub("\\.", "@", rownames(trg))
for (i in seq_len(nrow(trg))) {
  if (trg[i, "type"] == "CDS") {
    trg[i, "type"] <- "CDS"
  } else if (substr(trg[i, "type"], 1, 2) == "gS") {
    trg[i, "type"] <- "gS"
  } else if (substr(trg[i, "type"], 1, 3) == "gNS") {
    trg[i, "type"] <- "gNS"
  } else {
    trg[i, "type"] <- "no match"
  }
}


# Read the tree
tree <- read.tree(tree)
# Plot the tree
p <- ggtree(tree, branch.length = "none") +
  geom_tippoint(aes(subset = (label == focal_species)), shape = 23,
                color = "black", fill = "yellow",
                size = 3, stroke = 1, alpha = 1)
p <- p + new_scale_fill()
p <- gheatmap(p, data = trg, offset = 0.1, width = 0.1)
p