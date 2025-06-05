library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
library(yaml)
library(dplyr)
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



# Functions
# Get the lca node for a given cluster
get_lca_node <- function(tree, tips) {
  tips <- intersect(tips, tree$tip.label)
  # If only one genome in the cluster, return tip point
  if (length(tips) == 1) {
    return(which(tree$tip.label == tips))
  } else if (length(tips) > 1) {
    return(getMRCA(tree, tips))
  } else {
    return(NA)  # No valid tips
  }
}



# Filenames
tree_file <- paths$tree_file
dense_dir <- paths$dense_dir
genomes_list <- paths$genomes_list
fa_dir <- paths$fa_dir
output_dir <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/18_clustering/"
good_candidates <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/good_candidates.txt"
cluster_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/18_clustering/good_candidates_clustering.tsv"



# Read genomes list
genomes <- scan(genomes_list, what = "", sep = " ")[1:116]
# Read list of good denovo
good_candidates <- scan(good_candidates, what = "", sep = "\n")
# Get the info for each genome
data <- data.frame(matrix(NA, nrow = length(genomes), ncol = 3))
denovo_genomes <- data.frame(matrix(NA, nrow = 0, ncol = 2))
colnames(denovo_genomes) <- c("denovo", "genome")
rownames(data) <- genomes
colnames(data) <- c("n_denovo", "n_good_candidates", "gc_perc")
for (g in genomes) {
  # Filenames
  dense_file <- paste0(dense_dir, g, "/denovogenes.tsv")
  ## Number of denovo
  denovo <- read.table(dense_file, header = TRUE, sep = "\t")
  denovo_names <- denovo$CDS
  # Add the match denovo/genome to dtf
  for (i in seq_along(denovo_names)) {
    denovo_name <- denovo_names[i]
    denovo_genomes <- rbind(denovo_genomes,
                            data.frame(denovo = denovo_name, genome = g))
  }
  n_denovo <- nrow(denovo)
  # Good candidates
  n_good_candidates <- sum(denovo_names %in% good_candidates)
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
  data[g, "gc_perc"] <- gc
  data[g, "n_good_candidates"] <- n_good_candidates
}
data$n_bad_candidates <- data$n_denovo - data$n_good_candidates
n_good <- sum(data$n_good_candidates)
n_bad <- sum(data$n_bad_candidates)



# Extract the clusters
cluster_df <- read.table(cluster_file, header = FALSE, sep = "\t")
colnames(cluster_df) <- c("cluster", "denovo")
# Add the genomes
cluster_df <- left_join(cluster_df, denovo_genomes, by = "denovo")
# Get the lcas
cluster_lcas <- cluster_df %>%
  group_by(cluster) %>%
  summarise(
    genomes_in_cluster = list(unique(genome)),
    .groups = 'drop'
  ) %>%
  mutate(
    lca_node = sapply(genomes_in_cluster, function(tips) get_lca_node(tree, tips))
  ) %>%
  filter(!is.na(lca_node))



# Read the tree
tree <- read.tree(tree_file)
# Plot the tree
p <- ggtree(tree, layout = "circular", branch.length = "none")
tree_data <- p$data

# Match LCA nodes to tree coordinates
cluster_lcas_coords <- cluster_lcas %>%
  left_join(tree_data, by = c("lca_node" = "node"))
# Keep only unique lcas and count
cluster_lcas_summary <- cluster_lcas_coords %>%
  group_by(lca_node, x, y) %>%
  summarise(
    count = n(),
    label = first(label),
    .groups = "drop"
  )



# Add the heatmaps
## Good denovo ##
p <- p + new_scale_fill()
p <- gheatmap(p, data[, "n_good_candidates", drop = FALSE], offset = 0.5,
              width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "#179207",
                      name = paste0("No integrity\npreservation\n(n = ", n_good, ")"),
                      guide = guide_colorbar(order = 1))

## Bad denovo ##
p <- p + new_scale_fill()
p <- gheatmap(p, data[, "n_bad_candidates", drop = FALSE], offset = 1.5,
              width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "#4f535a",
                      name = paste0("Integrity\npreserved\n(n =", n_bad, ")"),
                      guide = guide_colorbar(order = 2),
                      limits = c(0, max(data$n_bad_candidates)))

# Add the title
p <- p + ggtitle("Integrity of the denovo NC matches\nand location of the cluster LCAs") +
  theme(plot.title = element_text(hjust = 0.5, vjust = -10))

# Add the points
p <- p + geom_point(data = cluster_lcas_summary,
                    aes(x = x, y = y),
                    inherit.aes = FALSE, fill = "#c51515", shape = 21,
                    color = "black", stroke = 0.5) +
  scale_size(range = c(1.5, 4),
             name = paste0("Cluster LCA position\n(n =", nrow(cluster_lcas_coords), ")")) +
  theme(legend.position = "right")
p

ggsave(paste0(output_dir, "denovo_good_bad.png"))
