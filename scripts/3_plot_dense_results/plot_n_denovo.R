library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
library(yaml)
library(ggpubr)
library(viridis)
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
genera_dir <- paths$genera_dir
dense_dir <- paths$dense_dir
cds_dir <- paths$cds_dir
genomes_list <- paths$genomes_list
fa_dir <- paths$fa_dir
species_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/data/different_species.tsv"
output_dir <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/3_plot_dense_results/"



# Read genomes list
genomes <- scan(genomes_list, what = "", sep = " ")[1:116]
# Get the info for each genome
data <- data.frame(matrix(NA, nrow = length(genomes), ncol = 4))
rownames(data) <- genomes
colnames(data) <- c("n_denovo", "n_trg", "n_cds", "gc_perc")
for (g in genomes) {
  # Filenames
  dense_file <- paste0(dense_dir, g, "/denovogenes.tsv")
  genera_file <- Sys.glob(paste0(genera_dir, g, "/*gene_age_summary.tsv"))[1]
  cds_file <- paste0(cds_dir, g, "_CDS.faa")
  ## Number of denovo
  denovo <- read.table(dense_file, header = TRUE, sep = "\t")
  denovo_names <- denovo$CDS
  n_denovo <- nrow(denovo)
  ## Number of TRGs
  trgs <- read.table(genera_file, header = TRUE, sep = "\t", comment.char = "")
  n_trg <- sum(trgs[trgs$phylorank >= 7, 1])
  ## Number of CDS
  lines_fasta <- readLines(cds_file)
  n_cds <- sum(grepl("^>", lines_fasta))
  ## GC content
  fa_file <- paste0(fa_dir, g, ".fa")
  fasta <- readLines(fa_file)
  # get rid of > lines
  fasta <- fasta[!grepl(">", fasta)]
  # Get the number of nucleotides
  n_nucleotides <- sum(nchar(fasta))
  # Get the GC content
  n_gc <- sum(nchar(gsub("[^GCgc]", "", paste(fasta, collapse = ""))))
  gc <- n_gc / n_nucleotides * 100

  # Populate the data
  data[g, "n_denovo"] <- n_denovo
  data[g, "n_trg"] <- n_trg
  data[g, "n_cds"] <- n_cds
  data[g, "gc_perc"] <- gc
  data[g, "n_good_candidates"] <- n_good_candidates
}

# Distinct species
species_data <- read.table(species_file, header = TRUE, sep = "\t")
colnames(species_data) <- c("genome", "genre", "species.ani", "mash", "aai",
"rename")
data$species <- species_data[match(rownames(data), species_data$genome), "species.ani"]



########### Tree plot ###########
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
p <- gheatmap(p, data[, "n_trg", drop = FALSE], offset = 1,
              width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "red", name = "TRGs (#)",
                      guide = guide_colorbar(order = 2), limits = c(0, max(data$n_trg)))
## GC% ##
p <- p + new_scale_fill()
p <- gheatmap(p, data[, "gc_perc", drop = FALSE], offset = 2,
              width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "#009E73", high = "#6b00b2", name = "GC %",
                      guide = guide_colorbar(order = 3))

print(p)
ggsave(paste0(output_dir, "denovo_trg_116.png"))
## Species ##
p <- p + new_scale_fill()
p <- gheatmap(p, data[, "species", drop = FALSE], offset = 4,
              width = .1, colnames = FALSE, legend_title = NULL) +
  scale_fill_manual(guide = "none",
                    values = viridis::viridis(length(unique(data$species)), option = "turbo"))
print(p)
ggsave(paste0(output_dir, "denovo_trg_116_with_species.png"))