###################################################################
############################## IMPORT #############################
###################################################################

library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)
library(grid)

ggsave <- function(..., bg = "white",
                   width = 1000, height = 1000,
                   units = "px", dpi = 100) {
  ggplot2::ggsave(..., bg = bg,
                  width = width,
                  height = height,
                  units = units,
                  dpi = dpi)
}


input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/13_plot_dense_results/sequences/sequence_features_good_candidates_all.csv"
pval_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/13_plot_dense_results/sequences/2_bins/pvalues_2_bins.tsv"
bins_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/13_plot_dense_results/sequences/2_bins/bin_indexes_2.csv"
out_folder <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/13_plot_dense_results/sequences/2_bins"
#input_file <- "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/results/13_plot_dense_results/sequence_features_good_candidates_all.csv"
#pval_file <- "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/results/13_plot_dense_results/pvalues_2_bins.tsv"
#bins_file <- "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/results/13_plot_dense_results/sequences/2_bins/bin_indexes_2.csv"
#out_folder <- "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/results/13_plot_dense_results/sequences/2_bins"


data <- read.table(input_file, header = TRUE, sep = "\t")
pvals <- read.table(pval_file, header = TRUE, sep = "\t")
pvals$p <- pvals$pval
n_bins <- length(unique(pvals[c("bin1", "bin2")]))
# Pivot to longer format
descriptors <- setdiff(colnames(data), c("genome", "cds", "type"))
data <- pivot_longer(data, cols = all_of(descriptors), names_to = "feature", values_to = "value")
pvals <- pivot_longer(pvals, cols = all_of(descriptors), names_to = "feature", values_to = "pval")
# Add the bins to the data
bins <- read.table(bins_file, header = TRUE, sep = " ")
data <- data %>%
  left_join(bins, by = c("cds"))
data$bin <- as.character(data$bin)
# GC %
min_gc <- 39.4
max_gc <- 56.4
bin_limits <- seq(min_gc, max_gc, length.out = n_bins + 1)
bin_labels <- paste0(round(bin_limits[-length(bin_limits)], 1), " - ", round(bin_limits[-1], 1), " %")
signif_label <- "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05"




###################################################################
############################ FUNCTIONS ############################
###################################################################

# Get the y positions to print the p values
get_y_positions <- function(data, local_pvals, fact) {

# Get different possibilities
denovo_0 <- ((local_pvals$type1 == "denovo" & local_pvals$bin1 == "0") |
               (local_pvals$type2 == "denovo" & local_pvals$bin2 == "0"))
denovo_1 <- ((local_pvals$type1 == "denovo" & local_pvals$bin1 == "1") |
               (local_pvals$type2 == "denovo" & local_pvals$bin2 == "1"))
trg_0 <- ((local_pvals$type1 == "trg" & local_pvals$bin1 == "0") |
            (local_pvals$type2 == "trg" & local_pvals$bin2 == "0"))
trg_1 <- ((local_pvals$type1 == "trg" & local_pvals$bin1 == "1") |
            (local_pvals$type2 == "trg" & local_pvals$bin2 == "1"))
cds_0 <- ((local_pvals$type1 == "cds" & local_pvals$bin1 == "0") |
            (local_pvals$type2 == "cds" & local_pvals$bin2 == "0"))
cds_1 <- ((local_pvals$type1 == "cds" & local_pvals$bin1 == "1") |
            (local_pvals$type2 == "cds" & local_pvals$bin2 == "1"))
# Order rows
local_pvals <- local_pvals %>%
  mutate(order = case_when(
    cds_0 & denovo_1 ~ 1,
    trg_0 & denovo_1 ~ 2,
    trg_1 & cds_0 ~ 3,
    cds_0 & cds_1 ~ 4,
    trg_0 & trg_1 ~ 5,
    denovo_0 & denovo_1 ~ 6,
    trg_0 & cds_1 ~ 7,
    trg_1 & denovo_0 ~ 8,
    denovo_0 & trg_1 ~ 9,
    trg_0 & cds_0 ~ 10,
    cds_0 & denovo_0 ~ 11,
    denovo_0 & cds_1 ~ 12,
    cds_1 & denovo_1 ~ 13,
    cds_0 & trg_0 ~ 14,
    trg_0 & denovo_0 ~ 15,
    cds_1 & trg_1 ~ 16,
    trg_1 & denovo_1 ~ 17,
    TRUE ~ 18
  )) %>%
  arrange(order)

  # Get the top of the whisker and the p-val for each type
  types <- c("denovo", "trg", "cds")
  whisker_tops <- c()
  min_whisker_base <- max(data$value)
  y_mat <- matrix(NA, nrow = 6, ncol = 3)
  colnames(y_mat) <- c("type", "bin", "whisker_top")
  i <- 0
  for (type in types) {
    for (bin in c("0", "1")) {
      i <- i + 1
      all_val <- as.numeric(as.vector(data[data$type == type & data$bin == bin, "value"])$value)
      whisker_top <- min(IQR(all_val) * 1.5 + quantile(all_val, 0.75), max(all_val))
      whisker_base <- max(quantile(all_val, 0.25) - IQR(all_val) * 1.5, min(all_val))
      y_mat[i, ] <- c(type, bin, whisker_top)
      if (whisker_base < min_whisker_base) {
        min_whisker_base <- whisker_base
      }
    }
  }
  reference_val <- max(as.numeric(y_mat[, "whisker_top"])) - min_whisker_base

  # Add y position
  pos <- max(as.numeric(y_mat[, "whisker_top"])) + reference_val * fact
  for (i in 1:nrow(local_pvals)) {
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
  local_pvals <- pvals[pvals$feature == desc, c("type1", "type2", "bin1", "bin2", "p")]
  # Create combined x-axis groups
  local_pvals$group1 <- paste0(local_pvals$bin1, "_", local_pvals$type1)
  local_pvals$group2 <- paste0(local_pvals$bin2, "_", local_pvals$type2)
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




###################################################################
############################### MAIN ##############################
###################################################################

###### Sequence length ######
data_len <- data[data$feature == "length", ]
data_len$type <- factor(data_len$type, levels = c("cds", "trg", "denovo"))
data_len$group <- paste0(data_len$bin, "_", data_len$type)
data_len_summary <- data_len %>%
  group_by(bin, type) %>%
  summarise(n = n(), .groups = "drop")

ggplot(data_len, aes(x = bin, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  geom_text(data = data_len_summary,
            aes(x = bin,
                y = 65,
                label = paste0("n = ", n),
                group = type),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 4) +
  labs(title = "Sequence length distribution",
       x = "% GC (whole genome)",
       y = "Length (bp)") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  coord_trans(y = "log10") +
  scale_y_continuous(breaks = c(50, 100, 500, 1000, 2000),
                     labels = c("50", "100", "500", "1000", "2000")) +
  scale_x_discrete(labels = c("0" = bin_labels[1],
                              "1" = bin_labels[2])) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_blank())
#  stat_pvalue_manual(get_pvals("length", data_len, 0.5),
#                     label = "p.signif",
#                     inherit.aes = FALSE,
#                     hide.ns = TRUE) +
#  annotate("text", x = 3.3, y = max(data_len$value) * 1,
#           label = signif_label, hjust = 1, vjust = 1, 
#           size = 3, color = "black")
ggsave(paste0(out_folder, "/sequence_length.png"))



###### GC ratio ######
data_gc <- data[data$feature == "gc_rate", ]
data_gc$type <- factor(data_gc$type, levels = c("cds", "trg", "denovo"))
data_gc$group <- paste0(data_gc$bin, "_", data_gc$type)
data_gc_summary <- data_gc %>%
  group_by(bin, type) %>%
  summarise(n = n(), .groups = "drop")

ggplot(data_gc, aes(x = bin, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  geom_text(data = data_gc_summary,
            aes(x = bin,
                y = 0.5,
                label = paste0("n = ", n),
                group = type),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 4) +
  labs(title = "GC ratio distribution",
       x = "% GC (whole genome)",
       y = "GC ratio: sequence GC % / genome GC %") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  scale_x_discrete(labels = c("0" = bin_labels[1],
                              "1" = bin_labels[2])) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_blank())
#  stat_pvalue_manual(get_pvals("gc_rate", data_gc, 0.5),
#                     label = "p.signif",
#                     inherit.aes = FALSE,
#                     hide.ns = TRUE) +
#  annotate("text", x = 3.3, y = max(data_len$value) * 1,
#           label = signif_label, hjust = 1, vjust = 1,
#           size = 3, color = "black")
ggsave(paste0(out_folder, "/gc_content.png"))



###### GC ratio (intergenic) ######
data_gc_inter <- data[data$feature == "inter_gc_rate", ]
data_gc_inter$type <- factor(data_gc_inter$type, levels = c("cds", "trg", "denovo"))
data_gc_inter$group <- paste0(data_gc_inter$bin, "_", data_gc_inter$type)
data_gc_inter_summary <- data_gc_inter %>%
  group_by(bin, type) %>%
  summarise(n = n(), .groups = "drop")

ggplot(data_gc_inter, aes(x = bin, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  geom_text(data = data_gc_inter_summary,
            aes(x = bin,
                y = 0.6,
                label = paste0("n = ", n),
                group = type),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 4) +
  labs(title = "GC ratio distribution (intergenic)",
       x = "% GC (whole genome)",
       y = "GC ratio: sequence GC % / intergenic ORFs GC %") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  scale_x_discrete(labels = c("0" = bin_labels[1],
                              "1" = bin_labels[2])) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_blank())
ggsave(paste0(out_folder, "/gc_content_intergenic.png"))



###### Aromaticity ######
data_aro <- data[data$feature == "aromaticity", ]
data_aro$type <- factor(data_aro$type, levels = c("cds", "trg", "denovo"))
data_aro$group <- paste0(data_aro$bin, "_", data_aro$type)
data_aro_summary <- data_aro %>%
  group_by(bin, type) %>%
  summarise(n = n(), .groups = "drop")

ggplot(data_aro, aes(x = bin, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  geom_text(data = data_aro_summary,
            aes(x = bin,
                y = -0.02,
                label = paste0("n = ", n),
                group = type),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 4) +
  labs(title = "Aromaticity distribution",
       x = "% GC (whole genome)",
       y = "Aromaticity") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  scale_x_discrete(labels = c("0" = bin_labels[1],
                              "1" = bin_labels[2])) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_blank())
ggsave(paste0(out_folder, "/aromaticity.png"))



###### Instability index ######
data_inst <- data[data$feature == "instability", ]
data_inst$type <- factor(data_inst$type, levels = c("cds", "trg", "denovo"))
data_inst$group <- paste0(data_inst$bin, "_", data_inst$type)
data_inst_summary <- data_inst %>%
  group_by(bin, type) %>%
  summarise(n = n(), .groups = "drop")

ggplot(data_inst, aes(x = bin, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  geom_text(data = data_inst_summary,
            aes(x = bin,
                y = -13,
                label = paste0("n = ", n),
                group = type),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 4) +
  labs(title = "Instability index distribution",
       x = "% GC (whole genome)",
       y = "Instability index") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  scale_x_discrete(labels = c("0" = bin_labels[1],
                              "1" = bin_labels[2])) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_blank())
ggsave(paste0(out_folder, "/instability_index.png"))



###### Flexibility #######
data_flex <- data[data$feature == "mean_flexibility", ]
data_flex$type <- factor(data_flex$type, levels = c("cds", "trg", "denovo"))
data_flex$group <- paste0(data_flex$bin, "_", data_flex$type)
data_flex_summary <- data_flex %>%
  group_by(bin, type) %>%
  summarise(n = n(), .groups = "drop")

ggplot(data_flex, aes(x = bin, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  geom_text(data = data_flex_summary,
            aes(x = bin,
                y = 0.95,
                label = paste0("n = ", n),
                group = type),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 4) +
  labs(title = "Mean flexibility distribution",
       x = "% GC (whole genome)",
       y = "Mean flexibility") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  scale_x_discrete(labels = c("0" = bin_labels[1],
                              "1" = bin_labels[2])) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_blank())
ggsave(paste0(out_folder, "/mean_flexibility.png"))



###### Hydropathy ######
data_hydro <- data[data$feature == "hydropathy", ]
data_hydro$type <- factor(data_hydro$type, levels = c("cds", "trg", "denovo"))
data_hydro$group <- paste0(data_hydro$bin, "_", data_hydro$type)
data_hydro_summary <- data_hydro %>%
  group_by(bin, type) %>%
  summarise(n = n(), .groups = "drop")

ggplot(data_hydro, aes(x = bin, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  geom_text(data = data_hydro_summary,
            aes(x = bin,
                y = -1.9,
                label = paste0("n = ", n),
                group = type),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 4) +
  labs(title = "Hydropathy distribution",
       x = "% GC (whole genome)",
       y = "Hydropathy") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  scale_x_discrete(labels = c("0" = bin_labels[1],
                              "1" = bin_labels[2])) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_blank())
ggsave(paste0(out_folder, "/hydrophobicity.png"))


####### HCA ######
data_hca <- data[data$feature == "hca", ]
data_hca$type <- factor(data_hca$type, levels = c("cds", "trg", "denovo"))
data_hca$group <- paste0(data_hca$bin, "_", data_hca$type)
data_hca_summary <- data_hca %>%
  group_by(bin, type) %>%
  summarise(n = n(), .groups = "drop")

ggplot(data_hca, aes(x = bin, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  geom_text(data = data_hca_summary,
            aes(x = bin,
                y = -11,
                label = paste0("n = ", n),
                group = type),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 4) +
  labs(title = "HCA distribution",
       x = "% GC (whole genome)",
       y = "HCA") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  scale_x_discrete(labels = c("0" = bin_labels[1],
                              "1" = bin_labels[2])) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_blank())
ggsave(paste0(out_folder, "/hca.png"))



###### AA use ######
polar_aa <- c("S", "T", "N", "Q")
hydrophobic_aa <- c("V", "I", "L", "M", "F", "W", "Y")
positive_aa <- c("K", "R", "H")
negative_aa <- c("D", "E")
pg_aa <- c("G", "P")
a_aa <- c("A")
c_aa <- c("C")
aa_types <- list(polar_aa, hydrophobic_aa, positive_aa, negative_aa, pg_aa, a_aa, c_aa)
aa_types_names <- c("polar", "hydrophobic", "positive", "negative", "proline-glycine", "alanine", "cysteine")

## Plot ##
for (i in seq_along(aa_types)) {
  aa_type <- aa_types[[i]]
  aa_type_name <- aa_types_names[[i]]
  aa_plots <- c()
  for (aa in aa_type) {
    data_aa <- data[data$feature == paste0(aa, "_use"), ]
    data_aa$type <- factor(data_aa$type, levels = c("cds", "trg", "denovo"))
    data_aa$group <- paste0(data_aa$bin, "_", data_aa$type)
    data_aa_summary <- data_aa %>%
      group_by(bin, type) %>%
      summarise(n = n(), .groups = "drop")

    p <- ggplot(data_aa, aes(x = bin, y = value, fill = type)) +
      geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
      labs(title = aa,
           x = "",
           y = "") +
      scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
      theme_minimal() +
      scale_x_discrete(labels = c("0" = bin_labels[1],
                                  "1" = bin_labels[2])) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
            axis.text.y = element_text(size = 16), legend.position = "none")

    # Add to the list of plots
    aa_plots <- c(aa_plots, list(p))
  }
  fig <- ggarrange(plotlist = aa_plots)
  annotate_figure(fig, bottom = text_grob("% GC (whole genome)\n", size = 14),
                  left = text_grob("% use", rot = 90, size = 14),
                  top = text_grob(paste("Amino-acid distribution:",
                                        aa_type_name),
                                  size = 18))
  ggsave(paste0(out_folder, "/", aa_type_name, "_use.png"))
}




