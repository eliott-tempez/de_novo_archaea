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
genera_dir <- paths$genera_dir
dense_dir <- paths$dense_dir
cds_dir <- paths$cds_dir
genomes_list <- paths$genomes_list
fa_dir <- paths$fa_dir
output_dir <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/13_plot_dense_results/global_results"
intergenic_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/10_analyse_intergenic/intergenic_lengths.tsv"


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
}

# Normalised de novo & trg number
#data$n_trg_norm <- data$n_trg / data$n_cds * 100

# Get the species with the GCA identifier
species_GCA <- genomes[which(grepl("GCA", rownames(data)))]

# Intergenic mean lengths
intergenic <- read.table(intergenic_file, header = TRUE, sep = "\t", row.names = 1)



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
## TRGs normalised ##
p <- p + new_scale_fill()
p <- gheatmap(p, data[, "gc_perc", drop = FALSE], offset = 2,
              width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "#009E73", high = "#6b00b2", name = "GC %",
                      guide = guide_colorbar(order = 3))

# Add the title
p <- p + ggtitle("Number of de novo genes and TRGs for each genome") +
  theme(plot.title = element_text(hjust = 0.5, vjust = -10))
p

ggsave(paste0(output_dir, "denovo_trg_116.png"))



########### Tree plot with GCA ###########
# Read the tree
tree <- read.tree(tree_file)
# Plot the tree
p <- ggtree(tree, layout = "circular", branch.length = "none") +
  geom_tippoint(aes(subset = (label %in% species_GCA)), shape = 23,
                color = "black", fill = "yellow",
                size = 3, stroke = 1, alpha = 1)
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
                      guide = guide_colorbar(order = 3))

# Add the title
p <- p + ggtitle("Number of de novo genes and TRGs for each genome") +
  theme(plot.title = element_text(hjust = 0.5, vjust = -10))
p

ggsave(paste0(output_dir, "denovo_trg_116_GCA.png"))



########## Distribution of de novo genes ##########
ggplot(data, aes(x = n_denovo)) +
  geom_histogram(fill = "#8261dd", color = "black") +
  labs(title = "Distribution of the number of de novo genes\nfound for each genome",
       x = "Number of de novo genes",
       y = "Number of genomes") +
  theme(plot.title = element_text(hjust = 0.5, vjust = -10))
ggsave(paste0(output_dir, "denovo_distribution_116.png"))



########## Distribution of de novo genes with GCA ##########
data$status <- ifelse(rownames(data) %in% species_GCA, "GCA", "Other")
# Divide in 2 dataframes
data_GCA <- data[data$status == "GCA", ]
data_other <- data[data$status == "Other", ]
# Plot
ggplot() +
  geom_histogram(data = data_GCA, aes(x = n_denovo, fill = "GCA"), alpha = 0.6, color = "black") +
  geom_histogram(data = data_other, aes(x = n_denovo, fill = "Other"), alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("GCA" = "#5370f0", "Other" = "#e76ba5"), name = NULL) +
  labs(title = "Distribution of the number of de novo genes\nfound for each genome",
        x = "Number of de novo genes",
        y = "Number of genomes") +
  theme(plot.title = element_text(hjust = 0.5, vjust = -10))
ggsave(paste0(output_dir, "denovo_distribution_116_GCA.png"))




########## Distribution of intergenic lengths ##########
ggplot(intergenic, aes(x = mean_intergenic_length)) +
  geom_histogram(fill = "#8261dd", color = "black") +
  labs(title = "Distribution of the mean intergenic lengths",
       x = "Mean intergenic length (bp)",
       y = "Number of genomes") +
  theme(plot.title = element_text(hjust = 0.5, vjust = -10))
ggsave(paste0(output_dir, "intergenic_distribution_116.png"))








########## Distribution of intergenic lengths with GCA ##########
intergenic$status <- ifelse(rownames(intergenic) %in% species_GCA, "GCA", "Other")
# Divide in 2 dataframes
intergenic_GCA <- intergenic[intergenic$status == "GCA", ]
intergenic_other <- intergenic[intergenic$status == "Other", ]
# Plot
ggplot() +
  geom_histogram(data = intergenic_GCA, aes(x = mean_intergenic_length, fill = "GCA"), alpha = 0.6, color = "black") +
  geom_histogram(data = intergenic_other, aes(x = mean_intergenic_length, fill = "Other"), alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("GCA" = "#5370f0", "Other" = "#e76ba5"), name = NULL) +
  labs(title = "Distribution of the mean intergenic lengths",
       x = "Mean intergenic length (bp)",
       y = "Number of genomes") +
  theme(plot.title = element_text(hjust = 0.5, vjust = -10))
ggsave(paste0(output_dir, "intergenic_distribution_116_GCA.png"))