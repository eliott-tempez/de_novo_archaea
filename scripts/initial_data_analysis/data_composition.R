library(ggplot2)
library(dplyr)
library(tidyr)
ggsave <- function(..., bg = "white", width = 1000, height = 1000, units = "px", dpi = 100) {
  ggplot2::ggsave(..., bg = bg, width = width, height = height, units = units, dpi = dpi)
}

# Import the file
file <- "/home/eliott.tempez/Documents/archaea_data/Thermocomplete_table.csv"
data <- read.csv(file)

# Which columns are filled for at least 95% of all organisms
full_col <- c()
len_max <- dim(data)[1] * .95
for (col in colnames(data)) {
  print(col)
  print(data[, col])
  if (length(which(data[, col] != "")) > len_max) {
    full_col <- c(full_col, col)
  }
}
other_col <- setdiff(colnames(data), full_col)

# Delete non informative columns
data_high <- data[, setdiff(full_col, c("Genome", "Geographic.Region"))]
data_low <- cbind(data[,"Name"], data[, setdiff(other_col, c("Assembly.accession", "acces_EGM", "EGM_name"))])



## Start with columns for which we have a good amount of info
# Genome size
ggplot(data_high, aes(x = Genome_Size_Mb)) +
  geom_histogram(binwidth = 0.1, fill = "#A6C3AA", color = "black") +
  labs(x = "Genome size (Mb)", y = "Count", title = "Distribution of the genome sizes in the dataset") +
  theme_minimal()
# export figure
ggsave("/home/eliott.tempez/Documents/M2_Stage_I2BC/results/initial_data_analysis_figures/genome_size.png")


# GC content
ggplot(data_high, aes(x = GC_content.)) +
  geom_histogram(binwidth = 1, fill = "#C3A6AC", color = "black") +
  labs(x = "GC content (%)", y = "Count", title = "Distribution of the GC content in the dataset") +
  theme_minimal()
# export figure
ggsave("/home/eliott.tempez/Documents/M2_Stage_I2BC/results/initial_data_analysis_figures/GC_content.png")

# Average CDS and intergenic length
# Combine the data into a long format for ggplot
data_long <- data_high %>%
  select(Average_CDS_length, Average_intergenic_length) %>%
  pivot_longer(cols = everything(), names_to = "Type", values_to = "Length")

# Plot histograms with ggplot
ggplot(data_long, aes(x = Length, fill = Type)) +
  geom_histogram(binwidth = bin_width, position = "identity", alpha = 0.75, color = "black") +
  scale_fill_manual(values = c("Average_CDS_length" = "#C3A6C3", "Average_intergenic_length" = "#C3C3A6")) +
  labs(x = "Average length (nucl)", y = "Count", title = "Average CDS and intergenic region length in the dataset") +
  theme_minimal() +
  theme(legend.title = element_blank())
# export figure
ggsave("/home/eliott.tempez/Documents/M2_Stage_I2BC/results/initial_data_analysis_figures/dna_legth.png")