###################################################################
############################## IMPORT #############################
###################################################################

# Librairies
library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)
library(grid)


# User-defined parameters
n_bins <- 3
plot_pvals <- FALSE


# Files
input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/13_plot_dense_results/sequences/sequence_features_good_candidates_all.csv"
out_folder <- paste0("/home/eliott.tempez/Documents/M2_Stage_I2BC/results/13_plot_dense_results/sequences/", n_bins, "_bins/")
pval_fils <- paste0(out_folder, "pvalues_", n_bins, "_bins.tsv")
bins_file <- paste0(out_folder, "bin_indexes_", n_bins, ".csv")
#input_file <- "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/results/13_plot_dense_results/sequence_features_good_candidates_all.csv"
#pval_file <- "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/results/13_plot_dense_results/pvalues_2_bins.tsv"
#bins_file <- "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/results/13_plot_dense_results/sequences/2_bins/bin_indexes_2.csv"
#out_folder <- "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/results/13_plot_dense_results/sequences/2_bins"


# Read data
data <- read.table(input_file, header = TRUE, sep = "\t")
if (plot_pvals) {
  pvals <- read.table(pval_file, header = TRUE, sep = "\t")
}
# Pivot to longer format
descriptors <- setdiff(colnames(data), c("genome", "cds", "type"))
data <- pivot_longer(data, cols = all_of(descriptors), names_to = "feature", values_to = "value")
if (plot_pvals) {
  pvals <- pivot_longer(pvals, cols = all_of(descriptors), names_to = "feature", values_to = "pval")
  pvals$p <- pvals$pval
}
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



add_dummy_rows <- function(data, feature, n_bins) {
  # Filter the data
  data_local <- data[data$feature == feature, ]
  # Add unique groups
  data_local$group <- paste0(data_local$bin, "_", data_local$type)
  # Create dummy rows for blank spaces on the x axis
  n_dummy <- (n_bins - 1)
  data_levels <- c("0_cds", "0_trg", "0_denovo")
  for (i in 1:n_dummy) {
    blank <- paste0(i - 1, "_blank")
    cds <- paste0(i, "_cds")
    trg <- paste0(i, "_trg")
    denovo <- paste0(i, "_denovo")
    data_levels <- c(data_levels, blank, cds, trg, denovo)
  }
  data_local$group <- as.character(data_local$group)
  data_local$group <- factor(data_local$group, levels = data_levels)
  for (i in 1:n_dummy) {
    dummy <- data_local[1, ]
    dummy$value <- NA
    dummy$group <- paste0(i - 1, "_blank")
    data_local <- rbind(data_local, dummy)
  }
  return(data_local)
}



# Get the number of CDS for each condition
get_ncds_conditions <- function(data, n_bins) {
  # Get the names of the dummy rows
  dummy_rows <- c()
  for (i in 1:(n_bins - 1)) {
    dummy_rows <- c(dummy_rows, paste0(i - 1, "_blank"))
  }
  # Get the data
  data_len_summary <- data %>%
    filter(!(group %in% dummy_rows)) %>%
    group_by(bin, type, group) %>%
    summarise(n = n(), .groups = "drop")
  return(data_len_summary)
}





# Get the y positions to print the p values
get_y_positions <- function(data, local_pvals, fact, non_signif) {
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
    if (!non_signif) {
      if (local_pvals[i, "p"] <= 0.05) {
        local_pvals[i, "y.position"] <- pos
        pos <- pos + reference_val * fact
      } else {
        local_pvals[i, "y.position"] <- NA
      }
    } else {
      if (local_pvals[i, "p"] > 0.05) {
        local_pvals[i, "y.position"] <- pos
        pos <- pos + reference_val * fact
      } else {
        local_pvals[i, "y.position"] <- NA
      }
    }
  }
  return(local_pvals)
}


# Get the p-values in a matrix to print on the graph
get_pvals <- function(desc, data, fact, non_signif = FALSE) {
  # Remove blank row
  data <- data[data$group != "blank", ]
  local_pvals <- pvals[pvals$feature == desc, c("type1", "type2", "bin1", "bin2", "p")]
  # Create combined x-axis groups
  local_pvals$group1 <- paste0(local_pvals$bin1, "_", local_pvals$type1)
  local_pvals$group2 <- paste0(local_pvals$bin2, "_", local_pvals$type2)
  # y position
  local_pvals <- get_y_positions(data, local_pvals, fact, non_signif)
  # significance level
  if (non_signif) {
    local_pvals$p.signif <- case_when(
      local_pvals$p >= 0.05 ~ "ns",
      TRUE ~ ""
    )
  } else {
    local_pvals$p.signif <- case_when(
    local_pvals$p <= 0.00001 ~ "****",
    local_pvals$p <= 0.0001 ~ "***",
    local_pvals$p <= 0.001 ~ "**",
    local_pvals$p <= 0.05 ~ "*",
    TRUE ~ "")
  }
  print(local_pvals)
  return(local_pvals)
}





###################################################################
############################### MAIN ##############################
###################################################################

###### Sequence length ######
data_len <- add_dummy_rows(data, "length", n_bins)
data_len_summary <- get_ncds_conditions(data_len, n_bins)

# Plot
ggplot(data_len, aes(x = group, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  geom_text(data = data_len_summary,
            aes(x = group, y = 60, label = paste0("n = ", n), group = type),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 4) +
  labs(title = "Sequence length distribution",
       x = "% GC (whole genome)",
       y = "Length (bp)") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  coord_trans(y = "log10") +
  scale_y_continuous(breaks = c(50, 100, 500, 1000, 2000),
                     labels = c("50", "100", "500", "1000", "2000")) +
  scale_x_discrete(labels = c("0_trg" = bin_labels[1],
                              "1_trg" = bin_labels[2],
                              "0_cds" = "", "0_denovo" = "",
                              "1_cds" = "", "1_denovo" = "",
                              "blank" = "")) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_blank()) 
  stat_pvalue_manual(get_pvals("length", data_len, .9, TRUE),
                     label = "p.signif",
                     inherit.aes = FALSE,
                     hide.ns = FALSE)
ggsave(paste0(out_folder, "/sequence_length.png"))



###### GC ratio ######
data_gc <- data[data$feature == "gc_rate", ]
data_gc$group <- paste0(data_gc$bin, "_", data_gc$type)
# Add a dummy level in the x-axis
data_gc$group <- as.character(data_gc$group)
data_gc$group <- factor(data_gc$group, levels = c("0_cds", "0_trg", "0_denovo", "blank", "1_cds", "1_trg", "1_denovo"))
# Create a dummy row for the blank space
dummy <- data_gc[1, ]
dummy$value <- NA
dummy$group <- "blank"
data_gc <- rbind(data_gc, dummy)

# Get summary data
data_gc_summary <- data_gc %>%
  filter(group != "blank") %>%  # skip the dummy in the summary
  group_by(bin, type, group) %>%
  summarise(n = n(), .groups = "drop")

# Plot
ggplot(data_gc, aes(x = group, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  geom_text(data = data_gc_summary,
            aes(x = group,
                y = 0.53,
                label = paste0("n = ", n),
                group = type),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 4) +
  labs(title = "GC ratio distribution",
       x = "% GC (whole genome)",
       y = "GC ratio: sequence GC % / genome GC %") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  scale_x_discrete(labels = c("0_trg" = bin_labels[1],
                              "1_trg" = bin_labels[2],
                              "0_cds" = "", "0_denovo" = "",
                              "1_cds" = "", "1_denovo" = "",
                              "blank" = "")) +
  scale_y_continuous(breaks = seq(0, 1.3, 0.2)) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_blank()) +
  stat_pvalue_manual(get_pvals("gc_rate", data_gc, 0.12, FALSE),
                     label = "p.signif",
                     inherit.aes = FALSE,
                     hide.ns = TRUE) +
  annotate("text", x = 4, y = 1.8,
           label = signif_label, hjust = 1, vjust = 1,
           size = 3, color = "black")
ggsave(paste0(out_folder, "/gc_content.png"))



###### GC ratio (intergenic) ######
data_gc_inter <- data[data$feature == "inter_gc_rate", ]
data_gc_inter$group <- paste0(data_gc_inter$bin, "_", data_gc_inter$type)
# Add a dummy level in the x-axis
data_gc_inter$group <- as.character(data_gc_inter$group)
data_gc_inter$group <- factor(data_gc_inter$group, levels = c("0_cds", "0_trg", "0_denovo", "blank", "1_cds", "1_trg", "1_denovo"))
# Create a dummy row for the blank space
dummy <- data_gc_inter[1, ]
dummy$value <- NA
dummy$group <- "blank"
data_gc_inter <- rbind(data_gc_inter, dummy)

# Get summary data
data_gc_inter_summary <- data_gc_inter %>%
  filter(group != "blank") %>%  # skip the dummy in the summary
  group_by(bin, type, group) %>%
  summarise(n = n(), .groups = "drop")

# Plot
ggplot(data_gc_inter, aes(x = group, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  geom_text(data = data_gc_inter_summary,
            aes(x = group,
                y = 0.62,
                label = paste0("n = ", n),
                group = type),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 4) +
  labs(title = "GC ratio distribution (intergenic)",
       x = "% GC (whole genome)",
       y = "GC ratio: sequence GC % / intergenic ORFs GC %") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  scale_x_discrete(labels = c("0_trg" = bin_labels[1],
                              "1_trg" = bin_labels[2],
                              "0_cds" = "", "0_denovo" = "",
                              "1_cds" = "", "1_denovo" = "",
                              "blank" = "")) +
  scale_y_continuous(breaks = seq(0, 1.6, 0.2)) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_blank()) +
  stat_pvalue_manual(get_pvals("inter_gc_rate", data_gc_inter, .1, FALSE),
                     label = "p.signif",
                     inherit.aes = FALSE,
                     hide.ns = TRUE) +
  annotate("text", x = 4, y = 2,
           label = signif_label, hjust = 1, vjust = 1,
           size = 3, color = "black")
ggsave(paste0(out_folder, "/gc_content_intergenic.png"))



###### Aromaticity ######
data_aro <- data[data$feature == "aromaticity", ]
data_aro$group <- paste0(data_aro$bin, "_", data_aro$type)
# Add a dummy level in the x-axis
data_aro$group <- as.character(data_aro$group)
data_aro$group <- factor(data_aro$group, levels = c("0_cds", "0_trg", "0_denovo", "blank", "1_cds", "1_trg", "1_denovo"))
# Create a dummy row for the blank space
dummy <- data_aro[1, ]
dummy$value <- NA
dummy$group <- "blank"
data_aro <- rbind(data_aro, dummy)

# Get summary data
data_aro_summary <- data_aro %>%
  filter(group != "blank") %>%  # skip the dummy in the summary
  group_by(bin, type, group) %>%
  summarise(n = n(), .groups = "drop")

# Compute p-values
pvals_aro <- get_pvals("aromaticity", data_aro, .1, FALSE)

# Plot
p <- ggplot(data_aro, aes(x = group, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  geom_text(data = data_aro_summary,
            aes(x = group,
                y = -0.02,
                label = paste0("n = ", n),
                group = type),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 4) +
  labs(title = "Aromaticity distribution",
       x = "% GC (whole genome)",
       y = "Aromaticity") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  scale_x_discrete(labels = c("0_trg" = bin_labels[1],
                              "1_trg" = bin_labels[2],
                              "0_cds" = "", "0_denovo" = "",
                              "1_cds" = "", "1_denovo" = "",
                              "blank" = "")) +
  scale_y_continuous(breaks = seq(0, 0.25, 0.05)) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_blank())

# Add stat_pvalue_manual and annotate only if y.position column is not entirely NA
if (!all(is.na(pvals_aro$y.position))) {
  p <- p +
    stat_pvalue_manual(pvals_aro,
                       label = "p.signif",
                       inherit.aes = FALSE,
                       hide.ns = TRUE) +
    annotate("text", x = 4, y = .27,
             label = signif_label, hjust = 1, vjust = 1,
             size = 3, color = "black")
}

# Print the plot
p
ggsave(paste0(out_folder, "/aromaticity.png"))



###### Instability index ######
data_inst <- data[data$feature == "instability", ]
data_inst$group <- paste0(data_inst$bin, "_", data_inst$type)
# Add a dummy level in the x-axis
data_inst$group <- as.character(data_inst$group)
data_inst$group <- factor(data_inst$group, levels = c("0_cds", "0_trg", "0_denovo", "blank", "1_cds", "1_trg", "1_denovo"))
# Create a dummy row for the blank space
dummy <- data_inst[1, ]
dummy$value <- NA
dummy$group <- "blank"
data_inst <- rbind(data_inst, dummy)

# Get summary data
data_inst_summary <- data_inst %>%
  filter(group != "blank") %>%  # skip the dummy in the summary
  group_by(bin, type, group) %>%
  summarise(n = n(), .groups = "drop")

# Plot
ggplot(data_inst, aes(x = group, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  geom_text(data = data_inst_summary,
            aes(x = group,
                y = -13,
                label = paste0("n = ", n),
                group = type),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 4) +
  labs(title = "Instability index distribution",
       x = "% GC (whole genome)",
       y = "Instability index") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  scale_x_discrete(labels = c("0_trg" = bin_labels[1],
                              "1_trg" = bin_labels[2],
                              "0_cds" = "", "0_denovo" = "",
                              "1_cds" = "", "1_denovo" = "",
                              "blank" = "")) +
  scale_y_continuous(breaks = seq(0, 75, 25)) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_blank()) +
  stat_pvalue_manual(get_pvals("instability", data_inst, .1, FALSE),
                     label = "p.signif",
                     inherit.aes = FALSE,
                     hide.ns = TRUE) +
  annotate("text", x = 4, y = 140,
           label = signif_label, hjust = 1, vjust = 1,
           size = 3, color = "black")
ggsave(paste0(out_folder, "/instability_index.png"))



###### Flexibility #######
data_flex <- data[data$feature == "mean_flexibility", ]
data_flex$group <- paste0(data_flex$bin, "_", data_flex$type)
# Add a dummy level in the x-axis
data_flex$group <- as.character(data_flex$group)
data_flex$group <- factor(data_flex$group, levels = c("0_cds", "0_trg", "0_denovo", "blank", "1_cds", "1_trg", "1_denovo"))
# Create a dummy row for the blank space
dummy <- data_flex[1, ]
dummy$value <- NA
dummy$group <- "blank"
data_flex <- rbind(data_flex, dummy)

# Get summary data
data_flex_summary <- data_flex %>%
  filter(group != "blank") %>%  # skip the dummy in the summary
  group_by(bin, type, group) %>%
  summarise(n = n(), .groups = "drop")

# Plot
ggplot(data_flex, aes(x = group, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  geom_text(data = data_flex_summary,
            aes(x = group,
                y = 0.95,
                label = paste0("n = ", n),
                group = type),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 4) +
  labs(title = "Mean flexibility distribution",
       x = "% GC (whole genome)",
       y = "Mean flexibility") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  scale_x_discrete(labels = c("0_trg" = bin_labels[1],
                              "1_trg" = bin_labels[2],
                              "0_cds" = "", "0_denovo" = "",
                              "1_cds" = "", "1_denovo" = "",
                              "blank" = "")) +
  scale_y_continuous(breaks = seq(0.96, 1.04, 0.02)) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_blank()) +
  stat_pvalue_manual(get_pvals("mean_flexibility", data_flex, .1, FALSE),
                     label = "p.signif",
                     inherit.aes = FALSE,
                     hide.ns = TRUE) +
  annotate("text", x = 4, y = 1.06,
           label = signif_label, hjust = 1, vjust = 1,
           size = 3, color = "black")
ggsave(paste0(out_folder, "/mean_flexibility.png"))



###### Hydropathy ######
data_hydro <- data[data$feature == "hydropathy", ]
data_hydro$group <- paste0(data_hydro$bin, "_", data_hydro$type)
# Add a dummy level in the x-axis
data_hydro$group <- as.character(data_hydro$group)
data_hydro$group <- factor(data_hydro$group, levels = c("0_cds", "0_trg", "0_denovo", "blank", "1_cds", "1_trg", "1_denovo"))
# Create a dummy row for the blank space
dummy <- data_hydro[1, ]
dummy$value <- NA
dummy$group <- "blank"
data_hydro <- rbind(data_hydro, dummy)

# Get summary data
data_hydro_summary <- data_hydro %>%
  filter(group != "blank") %>%  # skip the dummy in the summary
  group_by(bin, type, group) %>%
  summarise(n = n(), .groups = "drop")

# Plot
ggplot(data_hydro, aes(x = group, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  geom_text(data = data_hydro_summary,
            aes(x = group,
                y = -1.9,
                label = paste0("n = ", n),
                group = type),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 4) +
  labs(title = "Hydropathy distribution",
       x = "% GC (whole genome)",
       y = "Hydropathy") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  scale_x_discrete(labels = c("0_trg" = bin_labels[1],
                              "1_trg" = bin_labels[2],
                              "0_cds" = "", "0_denovo" = "",
                              "1_cds" = "", "1_denovo" = "",
                              "blank" = "")) +
  scale_y_continuous(breaks = seq(-2, 2, 1)) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_blank()) +
  stat_pvalue_manual(get_pvals("hydropathy", data_hydro, .1, FALSE),
                     label = "p.signif",
                     inherit.aes = FALSE,
                     hide.ns = TRUE) +
  annotate("text", x = 4, y = 3.5,
           label = signif_label, hjust = 1, vjust = 1,
           size = 3, color = "black")
ggsave(paste0(out_folder, "/hydrophobicity.png"))


####### HCA ######
data_hca <- data[data$feature == "hca", ]
data_hca$group <- paste0(data_hca$bin, "_", data_hca$type)
# Add a dummy level in the x-axis
data_hca$group <- as.character(data_hca$group)
data_hca$group <- factor(data_hca$group, levels = c("0_cds", "0_trg", "0_denovo", "blank", "1_cds", "1_trg", "1_denovo"))
# Create a dummy row for the blank space
dummy <- data_hca[1, ]
dummy$value <- NA
dummy$group <- "blank"
data_hca <- rbind(data_hca, dummy)

# Get summary data
data_hca_summary <- data_hca %>%
  filter(group != "blank") %>%  # skip the dummy in the summary
  group_by(bin, type, group) %>%
  summarise(n = n(), .groups = "drop")

# Plot
ggplot(data_hca, aes(x = group, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  geom_text(data = data_hca_summary,
            aes(x = group,
                y = -11,
                label = paste0("n = ", n),
                group = type),
            position = position_dodge(width = 0.75), vjust = -0.5, size = 4) +
  labs(title = "HCA distribution",
       x = "% GC (whole genome)",
       y = "HCA") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
  theme_minimal() +
  scale_x_discrete(labels = c("0_trg" = bin_labels[1],
                              "1_trg" = bin_labels[2],
                              "0_cds" = "", "0_denovo" = "",
                              "1_cds" = "", "1_denovo" = "",
                              "blank" = "")) +
  scale_y_continuous(breaks = seq(-10, 2, 4)) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 12),
        legend.title = element_blank()) +
  stat_pvalue_manual(get_pvals("hca", data_hca, .1, FALSE),
                     label = "p.signif",
                     inherit.aes = FALSE,
                     hide.ns = TRUE) +
  annotate("text", x = 4, y = 17,
           label = signif_label, hjust = 1, vjust = 1,
           size = 3, color = "black")
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
    data_aa$group <- paste0(data_aa$bin, "_", data_aa$type)
    # Add a dummy level in the x-axis
    data_aa$group <- as.character(data_aa$group)
    data_aa$group <- factor(data_aa$group, levels = c("0_cds", "0_trg", "0_denovo", "blank", "1_cds", "1_trg", "1_denovo"))
    # Create a dummy row for the blank space
    dummy <- data_aa[1, ]
    dummy$value <- NA
    dummy$group <- "blank"
    data_aa <- rbind(data_aa, dummy)

    # Get summary data
    data_aa_summary <- data_aa %>%
      filter(group != "blank") %>%  # skip the dummy in the summary
      group_by(bin, type, group) %>%
      summarise(n = n(), .groups = "drop")
    
    # Compute pvals 
    pval_aa <- get_pvals(paste0(aa, "_use"), data_aa, .3, FALSE)

    # Plot
    p <- ggplot(data_aa, aes(x = group, y = value, fill = type)) +
      geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
      labs(title = aa,
           x = "",
           y = "") +
      scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
      theme_minimal() +
      scale_x_discrete(labels = c("0_trg" = bin_labels[1],
                                  "1_trg" = bin_labels[2],
                                  "0_cds" = "", "0_denovo" = "",
                                  "1_cds" = "", "1_denovo" = "",
                                  "blank" = "")) +
      scale_y_continuous(breaks = seq(0, 26, 4)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
            axis.text.y = element_text(size = 16), legend.position = "none",
            axis.text.x = element_text(size = 14)) +
      theme(panel.border = element_rect(colour = "#ada9a9", fill = NA, size = 1))

      if (!all(is.na(pval_aa$y.position))) {
        p <- p + stat_pvalue_manual(get_pvals(paste0(aa, "_use"), data_aa, .1, FALSE),
                          label = "p.signif",
                          inherit.aes = FALSE,
                          hide.ns = TRUE)
      }

    # Add to the list of plots
    aa_plots <- c(aa_plots, list(p))
  }
  fig <- ggarrange(plotlist = aa_plots, common.legend = TRUE, legend = "right")
  annotate_figure(fig, bottom = text_grob("% GC (whole genome)\n", size = 14),
                  left = text_grob("% use", rot = 90, size = 14),
                  top = text_grob(paste("Amino-acid distribution:",
                                        aa_type_name, "\n"),
                                  size = 18),
                  right = text_grob(paste0("\n", signif_label), rot = 90, size = 10))
  ggsave(paste0(out_folder, "/aa_", aa_type_name, "_use.png"))
}




