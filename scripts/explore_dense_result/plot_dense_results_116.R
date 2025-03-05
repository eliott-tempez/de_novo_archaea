library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
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
genera_dir <- "/home/eliott.tempez/Documents/archaea_data/genera/out/"
dense_dir <- "/home/eliott.tempez/Documents/archaea_data/dense/"
cds_dir <- "/home/eliott.tempez/Documents/archaea_data/complete_122/CDS/"
genomes_list <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/genera_archaea/genomes_list.txt"
output_dir <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/explore_dense_results/"


# Read genomes list
genomes <- scan(genomes_list, what = "", sep = " ")[1:116]
# Get the info for each genome
data <- data.frame(matrix(NA, nrow = length(genomes), ncol = 3))
rownames(data) <- genomes
colnames(data) <- c("n_denovo", "n_trg", "n_cds")
for (g in genomes) {
  # Filenames
  dense_file <- paste0(dense_dir, g, "/denovogenes.tsv")
  genera_file <- Sys.glob(paste0(genera_dir, g, "/*gene_age_summary.tsv"))[1]
  cds_file <- paste0(cds_dir, g, "_CDS.faa")
  # Number of denovo
  denovo <- read.table(dense_file, header = TRUE, sep = "\t")
  n_denovo <- nrow(denovo)
  # Number of TRGs
  trgs <- read.table(genera_file, header = TRUE, sep = "\t", comment.char = "")
  n_trg <- sum(trgs[trgs$phylorank >= 7, 1])
  # Number of CDS
  lines_fasta <- readLines(cds_file)
  n_cds <- sum(grepl("^>", lines_fasta))
  # Populate the data
  data[g, "n_denovo"] <- n_denovo
  data[g, "n_trg"] <- n_trg
  data[g, "n_cds"] <- n_cds
}

# Normalised de novo & trg number
data$n_trg_norm <- data$n_trg / data$n_cds * 100


# Plot
# Read the tree
tree <- read.tree(tree_file)
# Plot the tree
p <- ggtree(tree, layout = "circular", branch.length = "none")
# Add the heatmaps
## De novo ##
p <- p + new_scale_fill()
p <- gheatmap(p, data[, "n_denovo", drop = FALSE],
              width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "black", name = "De novo (#)",
                      guide = guide_colorbar(order = 1))
## TRGs ##
p <- p + new_scale_fill()
p <- gheatmap(p, data[, "n_trg", drop = FALSE], offset = 2,
              width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "red", name = "TRGs (#)",
                      guide = guide_colorbar(order = 2), limits = c(0, max(data$n_trg)))
## TRGs normalised ##
p <- p + new_scale_fill()
p <- gheatmap(p, data[, "n_trg_norm", drop = FALSE], offset = 3,
              width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "#009E73", high = "#6b00b2", name = "TRGs\nnormalised (%)",
                      guide = guide_colorbar(order = 3), limits = c(0, 100))

# Add the title
p <- p + ggtitle("Number of de novo genes and TRGs for each genome") +
  theme(plot.title = element_text(hjust = 0.5, vjust = -10))
p

ggsave(paste0(output_dir, "denovo_trg_116.png"))