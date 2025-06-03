library(dplyr)
library(ggplot2)
# Save a plot
ggsave <- function(..., bg = "white",
                   width = 1000, height = 1000,
                   units = "px", dpi = 100) {
  ggplot2::ggsave(..., bg = bg,
                  width = width,
                  height = height,
                  units = units,
                  dpi = dpi)
}




descriptors_file <- input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/17_compare_denovo_sequences/sequence_features_good_candidates_all.csv"
focal_genome <- "GCA_001484685@Thermococcus_sp_2319x1"
good_candidates_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/good_candidates.txt"
out_folder <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/19_case_study/"


# Read the descriptors file
data <- read.csv(descriptors_file, sep = "\t", header = TRUE)
# Keep only the genomes that have denovo genes
genomes <- unique(data[data$type == "denovo", "genome"])
data <- data %>%
  filter(genome %in% genomes)
descriptors <- setdiff(colnames(data), c("genome", "cds", "type"))

# Read the good candidates file
denovo <- scan(good_candidates_file, what = "", sep = "\n")


# Make plot for each descriptor
for (descriptor in descriptors) {
  # Extract the data
  data_local <- data[, c("genome", "cds", descriptor)]
  # Remove the iorfs
  data_local <- data_local[!grepl("iorf", data_local$cds), ]
  # Add the genome
  data_local$genome <- ifelse(data_local$genome == focal_genome, "Focal genome", "Other genomes")
  # Focal genome median
  focal_median <- median(data_local[data_local$genome == "Focal genome", descriptor], na.rm = TRUE)


  # Create the plot
  p <- ggplot(data_local, aes(y = .data[[descriptor]])) +
    geom_boxplot(outliers = FALSE, width = 0.5, fill = "#8a8888") +
    xlim(c(-0.5, 0.5)) +
    geom_hline(aes(yintercept = focal_median, color = "Case study\ngenome median"), linetype = "dashed") +
    labs(color = "") +
    scale_color_manual(values = c("Case study\ngenome median" = "red")) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ggtitle("Feature distribution for all CDSs in the\ngenomes containing denovo genes")
    

  print(p)

  # Save the plot
  ggsave(paste0(out_folder, descriptor, ".png"))
}
