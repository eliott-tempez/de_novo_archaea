library(yaml)
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


#config <- yaml.load_file("/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/my_functions/filepaths.yaml")
config <- yaml.load_file("/home/eliott/Documents/UNI/M2/Stage/M2_stage_I2BC/scripts/my_functions/filepaths.yaml")
#paths <- config$local_paths
paths <- config$home_paths
tree_file <- paths$tree_file
#input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/20_de_novo_conservation/de_novo_conservation.tsv"
#input_file_orf <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/20_de_novo_conservation/orf_conservation.tsv"
#output_folder <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/20_de_novo_conservation/figures"
input_file <- "/home/eliott/Documents/UNI/M2/Stage/M2_stage_I2BC/results/20_de_novo_conservation/de_novo_conservation.tsv"
input_file_orf <- "/home/eliott/Documents/UNI/M2/Stage/M2_stage_I2BC/results/20_de_novo_conservation/orf_conservation.tsv"
output_folder <- "/home/eliott/Documents/UNI/M2/Stage/M2_stage_I2BC/results/20_de_novo_conservation/figures"


# Read the de novo file
conservation_data <- read.table(input_file, header = TRUE, sep = "\t")
unique_denovo <- unique(conservation_data$denovo)
# Read the ORF conservation file
orf_conservation_data <- read.table(input_file_orf, header = TRUE, sep = "\t")
# Read the tree
tree <- read.tree(tree_file)


# For each de novo gene
for (denovo in unique_denovo) {
  # Filter the data for the current de novo gene
  denovo_data <- conservation_data[conservation_data$denovo == denovo, ]
  # Use the neighbour genome as rowname
  rownames(denovo_data) <- denovo_data$neighbour_genome
  # Get the focal genome
  focal_genome <- denovo_data$genome[1]
  # Add match as "cds" for the focal genome
  denovo_data[focal_genome, ] <- c(denovo, focal_genome, focal_genome, "cds")

  # Filter the orf conservation info
  orf_data <- orf_conservation_data[orf_conservation_data$genome == focal_genome, ]
  rownames(orf_data) <- orf_data$neighbour_genome
  orf_data[focal_genome, ] <- c(focal_genome, focal_genome, 0, 100, 0, 0)
  
  # Convert numeric columns to numeric (needed because c() converts everything to character)
  orf_data$n_cds <- as.numeric(orf_data$n_cds)
  orf_data$n_integral_nc <- as.numeric(orf_data$n_integral_nc)
  orf_data$n_fragmented_nc <- as.numeric(orf_data$n_fragmented_nc)
  orf_data$n_no_match <- as.numeric(orf_data$n_no_match)
  # Define the colors for the status
  status_colors <- c("cds" = "#e67b00",
                     "no match" = "#9e9e9e",
                     "integral nc" = "#00799e",
                     "fragmented nc" = "#8261dd")

  # Create the tree
  p <- ggtree(tree, layout = "circular", branch.length = "none")

  # Add the heatmap for the de novo data
  p <- p + new_scale_fill()
  p <- gheatmap(p, denovo_data[tree$tip.label, "status", drop = FALSE], offset = 1,
                width = .05, colnames = FALSE) +
    scale_fill_manual(
      values = status_colors,
      breaks = names(status_colors),
      guide = "none"
    )

  # Add the heatmaps for the ORF conservation data
  p <- p + new_scale_fill()
  p <- gheatmap(p, orf_data[tree$tip.label, "n_cds",
                            drop = FALSE], offset = 3,
                width = .05, colnames = FALSE) +
    scale_fill_gradient(low = "white", high = "#e67b00",
                        name = "% of CDS",
                        guide = guide_colorbar(order = 1),
                        limits = c(0, 100))

  p <- p + new_scale_fill()
  p <- gheatmap(p, orf_data[tree$tip.label, "n_integral_nc",
                            drop = FALSE], offset = 4,
                width = .05, colnames = FALSE) +
    scale_fill_gradient(low = "white", high = "#00799e",
                        name = "% of nc\n(integral)",
                        guide = guide_colorbar(order = 2),
                        limits = c(0, 100))

  p <- p + new_scale_fill()
  p <- gheatmap(p, orf_data[tree$tip.label, "n_fragmented_nc",
                            drop = FALSE], offset = 5,
                width = .05, colnames = FALSE) +
    scale_fill_gradient(low = "white", high = "#8261dd",
                        name = "% of nc\n(fragmented)",
                        guide = guide_colorbar(order = 3),
                        limits = c(0, 100))

  p <- p + new_scale_fill()
  p <- gheatmap(p, orf_data[tree$tip.label, "n_no_match",
                            drop = FALSE], offset = 6,
                width = .05, colnames = FALSE) +
    scale_fill_gradient(low = "white", high = "#9e9e9e",
                        name = "% of no match",
                        guide = guide_colorbar(order = 4),
                        limits = c(0, 100))


  # Add a diamond for the focal genome
  p <- p +  geom_tippoint(aes(subset = (label == focal_genome)),
                          shape = 23, color = "black",
                          fill = "yellow", size = 3,
                          stroke = 1, alpha = 1,
                          show.legend = FALSE)

  # Add a title
  denovo_name <- gsub("_gene_mRNA", "", denovo)
  p <- p + ggtitle(paste("Conservation of de novo gene", denovo_name,
                         "\n(internal circle) and 100 random ORFs\nin focal genome (external circles)")) +
    theme(plot.title = element_text(hjust = 0, margin = margin(l = 20), size = 16, face = "bold"))

  ggsave(paste0(output_folder, "/", denovo, ".png"))
}
