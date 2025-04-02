library(ggplot2)
library(tidyr)
ggsave <- function(..., bg = "white",
                   width = 1000, height = 1000,
                   units = "px", dpi = 100) {
  ggplot2::ggsave(..., bg = bg,
                  width = width,
                  height = height,
                  units = units,
                  dpi = dpi)
}

input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/13_plot_dense_results/sequences/sequence_features.csv"
out_folder <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/13_plot_dense_results/sequences"
data <- read.table(input_file, header = TRUE, sep = "\t")
n_cds <- nrow(data[data$type == "cds"])
n_trg <- nrow(data[data$type == "trg"])
n_denovo <- nrow(data[data$type == "denovo"])

# Pivot to longer
data$type <- data$origin
data <- pivot_longer(data, cols = c("gc_content", "aromaticity", "instability", "mean_flexibility", "hydropathy", "len_nu"), names_to = "feature", values_to = "value")


##### Sequence length #####
data_len <- data[data$feature == "gene_length", ]
data_len$origin <- factor(data_len$origin, levels = c("cds", "trg", "denovo"))
ggplot(data_len, aes(x = origin, y = value, fill = origin)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c") +
  labs(title = "Sequence length distribution",
       x = "Sequence type",
       y = "Length (bp)") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave(paste0(out_folder, "/sequence_length.png"))



##### GC content #####
data_gc <- data[data$feature == "gc", ]
data_gc$origin <- factor(data_gc$origin, levels = c("cds", "trg", "denovo"))
ggplot(data_gc, aes(x = origin, y = value, fill = origin)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c") +
  labs(title = "GC content distribution",
       x = "Sequence type",
       y = "GC content (%)") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave(paste0(out_folder, "/gc_content.png"))



##### Aromaticity #####
data_aro <- data[data$feature == "aromaticity", ]
data_aro$origin <- factor(data_aro$origin, levels = c("cds", "trg", "denovo"))
ggplot(data_aro, aes(x = origin, y = value, fill = origin)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c") +
  labs(title = "Aromaticity distribution",
       x = "Sequence type",
       y = "Aromaticity") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave(paste0(out_folder, "/aromaticity.png"))



##### Instability index #####
data_inst <- data[data$feature == "instability", ]
data_inst$origin <- factor(data_inst$origin, levels = c("cds", "trg", "denovo"))
ggplot(data_inst, aes(x = origin, y = value, fill = origin)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c") +
  labs(title = "Instability index distribution",
       x = "Sequence type",
       y = "Instability index") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave(paste0(out_folder, "/instability_index.png"))



##### Flexibility #####
data_flex <- data[data$feature == "flexibility", ]
data_flex$origin <- factor(data_flex$origin, levels = c("cds", "trg", "denovo"))
ggplot(data_flex, aes(x = origin, y = value, fill = origin)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c") +
  labs(title = "Mean flexibility distribution",
       x = "Sequence type",
       y = "Mean flexibility") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))



##### Hydrophobicity #####
data_hydro <- data[data$feature == "hydropathy", ]
data_hydro$origin <- factor(data_hydro$origin, levels = c("cds", "trg", "denovo"))
ggplot(data_hydro, aes(x = origin, y = value, fill = origin)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c") +
  labs(title = "Mean hydrophobicity distribution",
       x = "Sequence type",
       y = "Mean hydrophobicity") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave(paste0(out_folder, "/hydrophobicity.png"))