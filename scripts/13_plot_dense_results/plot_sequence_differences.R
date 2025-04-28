library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)
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


#input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/13_plot_dense_results/sequences/sequence_features_good_candidates_all.csv"
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
colnames(pvals) <- c("group1", "group2", "feature", "p")



# Get the y positions to print the p values
get_y_positions <- function(data, local_pvals, fact) {
  # Order rows
  local_pvals <- local_pvals %>%
    mutate(order = case_when(
      group1 == "trg" & group2 == "cds" ~ 1,
      group1 == "denovo" & group2 == "trg" ~ 2,
      group1 == "denovo" & group2 == "cds" ~ 3,
      TRUE ~ 4
    )) %>%
    arrange(order) %>%
    select(-order)
  # Get the top of the whisker and the p-val for each type
  types <- c("denovo", "trg", "cds")
  whisker_tops <- c()
  min_whisker_base <- max(data$value)
  for (type in types) {
    all_val <- as.numeric(as.vector(data[data$type == type, "value"])$value)
    whisker_top <- min(IQR(all_val) * 1.5 + quantile(all_val, 0.75), max(all_val))
    whisker_base <- max(quantile(all_val, 0.25) - IQR(all_val) * 1.5, min(all_val))
    whisker_tops <- c(whisker_tops, whisker_top)
    if (whisker_base < min_whisker_base) {
      min_whisker_base <- whisker_base
    }
  }
  reference_val <- max(whisker_tops) - min_whisker_base
  y_mat <- cbind("type" = types, "whisker_top" = whisker_tops)
  print(y_mat)
  # Add y position
  pos <- max(as.numeric(y_mat[, "whisker_top"])) + reference_val * fact
  for (i in 1:3) {
    if (local_pvals[i, "p"] <= 0.05) {
      local_pvals[i, "y.position"] <- pos
      pos <- pos + reference_val * fact
    } else {
      local_pvals[i, "y.position"] <- NA
    }
  }
  return(local_pvals)
}


# Get the p-values in a matrix to print on the graph
get_pvals <- function(desc, data, fact) {
  local_pvals <- pvals[pvals$feature == desc, c("group1", "group2", "p")]
  # y position
  local_pvals <- get_y_positions(data, local_pvals, fact)
  # significance level
  local_pvals$p.signif <- case_when(
    local_pvals$p <= 0.00001 ~ "****",
    local_pvals$p <= 0.0001 ~ "***",
    local_pvals$p <= 0.001 ~ "**",
    local_pvals$p <= 0.05 ~ "*",
    TRUE ~ ""
  )
  print(local_pvals)
  return(local_pvals)
}


# Pivot to longer
descriptors <- setdiff(colnames(data), c("genome", "cds", "type"))
data <- pivot_longer(data, cols = all_of(descriptors), names_to = "feature", values_to = "value")


##### Sequence length #####
data_len <- data[data$feature == "length", ]
data_len$type <- factor(data_len$type, levels = c("cds", "trg", "denovo"))
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
  scale_y_continuous(breaks = c(50, 100, 500, 1000, 2000),
                     labels = c("50", "100", "500", "1000", "2000")) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"))) +
  stat_pvalue_manual(get_pvals("length", data_len, 0.4), label = "p.signif", inherit.aes = FALSE, hide.ns = TRUE) +
  annotate("text", x = 3.3, y = max(data_len$value) * 1.5,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
#ggsave(paste0(out_folder, "/sequence_length.png"))



##### GC content #####
data_gc <- data[data$feature == "gc_rate", ]
data_gc$type <- factor(data_gc$type, levels = c("cds", "trg", "denovo"))
ggplot(data_gc, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
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
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"))) +
  stat_pvalue_manual(get_pvals("gc_rate", data_gc, 0.1), label = "p.signif", inherit.aes = FALSE, hide.ns = "p") +
  annotate("text", x = 3.3, y = max(data_gc$value) * 1.1,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
#ggsave(paste0(out_folder, "/gc_content.png"))



##### Aromaticity #####
data_aro <- data[data$feature == "aromaticity", ]
data_aro$type <- factor(data_aro$type, levels = c("cds", "trg", "denovo"))
ggplot(data_aro, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
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
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"))) +
  stat_pvalue_manual(get_pvals("aromaticity", data_aro, 0.1), label = "p.signif", inherit.aes = FALSE, hide.ns = "p") + 
  annotate("text", x = 3.3, y = max(data_aro$value) * 1.1,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
#ggsave(paste0(out_folder, "/aromaticity.png"))



##### Instability index #####
data_inst <- data[data$feature == "instability", ]
data_inst$type <- factor(data_inst$type, levels = c("cds", "trg", "denovo"))
ggplot(data_inst, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
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
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"))) +
  stat_pvalue_manual(get_pvals("instability", data_inst, 0.15), label = "p.signif", inherit.aes = FALSE, hide.ns = "p") +
  annotate("text", x = 3.3, y = max(data_inst$value) * 1.1,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
#ggsave(paste0(out_folder, "/instability_index.png"))



##### Flexibility #####
data_flex <- data[data$feature == "mean_flexibility", ]
data_flex$type <- factor(data_flex$type, levels = c("cds", "trg", "denovo"))
ggplot(data_flex, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
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
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"))) +
  stat_pvalue_manual(get_pvals("mean_flexibility", data_flex, 0.2), label = "p.signif", inherit.aes = FALSE, hide.ns = "p") +
  annotate("text", x = 3.3, y = max(data_flex$value) * 1.05,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
#ggsave(paste0(out_folder, "/mean_flexibility.png"))



##### Hydrophobicity #####
data_hydro <- data[data$feature == "hydropathy", ]
data_hydro$type <- factor(data_hydro$type, levels = c("cds", "trg", "denovo"))
ggplot(data_hydro, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
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
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"))) +
  stat_pvalue_manual(get_pvals("hydropathy", data_hydro, 0.2), label = "p.signif", inherit.aes = FALSE, hide.ns = "p") +
  annotate("text", x = 3.3, y = max(data_hydro$value) * 1.05,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
#ggsave(paste0(out_folder, "/hydrophobicity.png"))
