library(ggplot2)
library(tidyr)
library(ggpubr)
library(rstatix)
ggsave <- function(..., bg = "white",
                   width = 1000, height = 1000,
                   units = "px", dpi = 100) {
  ggplot2::ggsave(..., bg = bg,
                  width = width,
                  height = height,
                  units = units,
                  dpi = dpi)
}

#input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/13_plot_dense_results/sequences/sequence_features_good_candidates.csv"
#pval_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/13_plot_dense_results/sequences/pvalues.tsv"
input_file <- "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/results/13_plot_dense_results/sequences/sequence_features_good_candidates.csv"
pval_file <- "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/results/13_plot_dense_results/sequences/pvalues.tsv"
out_folder <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/13_plot_dense_results/sequences"
data <- read.table(input_file, header = TRUE, sep = "\t")
n_cds <- nrow(data[data$type == "cds", ])
n_trg <- nrow(data[data$type == "trg", ])
n_denovo <- nrow(data[data$type == "denovo", ])
# Import pvalues
pvals <- read.table(pval_file, header = TRUE, sep = "\t")
# Fix colnames
colnames(pvals) <- c("type_1", "type_2", "feature", "pval")

# Pivot to longer
data <- pivot_longer(data, cols = c("gc_content", "aromaticity", "instability", "mean_flexibility", "hydropathy", "len_nu"), names_to = "feature", values_to = "value")


##### Sequence length #####
data_len <- data[data$feature == "len_nu", ]
data_len$type <- factor(data_len$type, levels = c("cds", "trg", "denovo"))
pval_denovo_trg <- pvals[pvals$type_1 == "denovo" & pvals$type_2 == "trg" & pvals$feature == "length", "pval"]
pval_denovo_cds <- pvals[pvals$type_1 == "denovo" & pvals$type_2 == "cds" & pvals$feature == "length", "pval"]
pval_trg_cds <- pvals[pvals$type_1 == "trg" & pvals$type_2 == "cds" & pvals$feature == "length", "pval"]
ggplot(data_len, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "Sequence length distribution",
       x = "Sequence type",
       y = "Length (bp)") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12)) +
  coord_trans(y = "log10") +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")")))
#ggsave(paste0(out_folder, "/sequence_length.png"))



##### GC content #####
data_gc <- data[data$feature == "gc_content", ]
data_gc$type <- factor(data_gc$type, levels = c("cds", "trg", "denovo"))
ggplot(data_gc, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c") +
  labs(title = "GC ratio distribution",
       x = "Sequence type",
       y = "GC ratio: sequence GC % / genome GC %") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")")))
#ggsave(paste0(out_folder, "/gc_content.png"))



##### Aromaticity #####
data_aro <- data[data$feature == "aromaticity", ]
data_aro$type <- factor(data_aro$type, levels = c("cds", "trg", "denovo"))
ggplot(data_aro, aes(x = type, y = value, fill = type)) +
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
        axis.text.y = element_text(size = 12)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")")))
#ggsave(paste0(out_folder, "/aromaticity.png"))



##### Instability index #####
data_inst <- data[data$feature == "instability", ]
data_inst$type <- factor(data_inst$type, levels = c("cds", "trg", "denovo"))
ggplot(data_inst, aes(x = type, y = value, fill = type)) +
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
        axis.text.y = element_text(size = 12)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")")))
#ggsave(paste0(out_folder, "/instability_index.png"))



##### Flexibility #####
data_flex <- data[data$feature == "mean_flexibility", ]
data_flex$type <- factor(data_flex$type, levels = c("cds", "trg", "denovo"))
ggplot(data_flex, aes(x = type, y = value, fill = type)) +
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
        axis.text.y = element_text(size = 12)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")")))
#ggsave(paste0(out_folder, "/mean_flexibility.png"))



##### Hydrophobicity #####
data_hydro <- data[data$feature == "hydropathy", ]
data_hydro$type <- factor(data_hydro$type, levels = c("cds", "trg", "denovo"))
ggplot(data_hydro, aes(x = type, y = value, fill = type)) +
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
        axis.text.y = element_text(size = 12)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")")))
#ggsave(paste0(out_folder, "/hydrophobicity.png"))
