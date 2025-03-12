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
genomes_dir <- "/home/eliott.tempez/Documents/archaea_data/complete_122/fasta_renamed/"
cds_dir <- "/home/eliott.tempez/Documents/archaea_data/complete_122/reannotated_CDS/"
genomes_list <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/genera_archaea/genomes_list.txt"


# Read genomes list
genomes <- scan(genomes_list, what = "", sep = " ")[1:116]
# Get the info for each genome
data <- data.frame(matrix(NA, nrow = length(genomes), ncol = 3))
rownames(data) <- genomes
colnames(data) <- c("gc_percent", "n_cds", "genome_size")
for (g in genomes) {
    # Filenames
    fasta_file <- paste0(genomes_dir, g, ".fa")
    cds_file <- paste0(cds_dir, g, "_CDS.faa")
    # GC content
    fasta <- readLines(fasta_file)
    # get rid of > lines
    fasta <- fasta[!grepl(">", fasta)]
    # Get the number of nucleotides
    n_nucleotides <- sum(nchar(fasta))
    # Get the GC content
    n_gc <- sum(nchar(gsub("[^GCgc]", "", paste(fasta, collapse = ""))))
    gc <- n_gc / n_nucleotides * 100
    # Number of CDS
    lines_fasta <- readLines(cds_file)
    n_cds <- sum(grepl("^>", lines_fasta))
    # Add the data
    data[g, "gc_percent"] <- gc
    data[g, "n_cds"] <- n_cds
    data[g, "genome_size"] <- n_nucleotides / 1e6  # Convert to megabases
}


# Plot the tree
tree <- read.tree(tree_file)
# Plot the tree
p <- ggtree(tree, layout = "circular", branch.length = "none")
# Add the heatmaps
# Genome length
p <- p + new_scale_fill()
p <- gheatmap(p, data[, "genome_size", drop = FALSE], width = 0.05, colnames = FALSE) + 
    scale_fill_gradient(name = "Genome size (Mb)", guide = guide_colorbar(order = 1), low = "#56B4E9", high = "#E69F00")
# Number of CDS
p <- p + new_scale_fill()
p <- gheatmap(p, data[, "n_cds", drop = FALSE], offset = 2, width = 0.05, colnames = FALSE) + 
    scale_fill_gradient(name = "Number of CDS", guide = guide_colorbar(order = 2), low = "#F0E442", high = "#009E73")
# GC content
p <- p + new_scale_fill()
p <- gheatmap(p, data[, "gc_percent", drop = FALSE], offset = 4, width = 0.05, colnames = FALSE) + 
    scale_fill_gradient(name = "GC content (%)", guide = guide_colorbar(order = 3), low = "#0072B2", high = "#CC79A7")
p

ggsave("/home/eliott.tempez/Documents/M2_Stage_I2BC/results/initial_data_analysis_figures/initial_tree.png")
