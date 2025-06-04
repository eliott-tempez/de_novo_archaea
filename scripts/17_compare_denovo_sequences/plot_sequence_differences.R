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


home <- TRUE


if (home) {
  input_file <- "/home/eliott/Documents/UNI/M2/Stage/M2_stage_I2BC/results/17_compare_denovo_sequences/sequence_features_good_candidates_all.csv"
  pval_file <- "/home/eliott/Documents/UNI/M2/Stage/M2_stage_I2BC/results/17_compare_denovo_sequences/1_bins/pvalues_1_bins.tsv"
  out_folder <- "/home/eliott/Documents/UNI/M2/Stage/M2_stage_I2BC/results/17_compare_denovo_sequences/1_bins/"
} else {
  input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/17_compare_denovo_sequences/sequence_features_good_candidates_all.csv"
  pval_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/17_compare_denovo_sequences/1_bins/pvalues_1_bins.tsv"
  out_folder <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/17_compare_denovo_sequences/1_bins/"
}



only_denovo_genomes <- FALSE
if (only_denovo_genomes) {
  out_folder <- paste0(out_folder, "only_denovo_genomes/")
}


data <- read.table(input_file, header = TRUE, sep = "\t")
# Keep only the denovo genomes
if (only_denovo_genomes) {
  data <- data[data$genome %in% unique(data[data$type == "denovo", "genome"]), ]
}
n_cds <- nrow(data[data$type == "cds", ])
n_trg <- nrow(data[data$type == "trg", ])
n_denovo <- nrow(data[data$type == "denovo", ])
n_iorf <- nrow(data[data$type == "iorf", ])
# Import pvalues
pvals <- read.table(pval_file, header = TRUE, sep = "\t")



# Get the y positions to print the p values
get_y_positions <- function(data, local_pvals, fact) {
  # Order rows
  local_pvals <- local_pvals %>%
    mutate(order = case_when(
      (group1 == "trg" & group2 == "cds") | (group1 == "cds" & group2 == "trg") ~ 1,
      (group1 == "trg" & group2 == "denovo") | (group1 == "denovo" & group2 == "trg") ~ 2,
      (group1 == "denovo" & group2 == "iorf") | (group1 == "iorf" & group2 == "denovo") ~ 3,
      (group1 == "cds" & group2 == "denovo") | (group1 == "denovo" & group2 == "cds") ~ 4,
      (group1 == "trg" & group2 == "iorf") | (group1 == "iorf" & group2 == "trg") ~ 5,
      (group1 == "cds" & group2 == "iorf") | (group1 == "iorf" & group2 == "cds") ~ 6,
      TRUE ~ 7
    )) %>%
    arrange(order) %>%
    select(-order)
  # Get the top of the whisker and the p-val for each type
  types <- c("denovo", "trg", "cds", "iorf")
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
  # Add y position
  pos <- max(as.numeric(y_mat[, "whisker_top"])) + reference_val * fact
  for (i in 1:6) {
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
  # Remove NAs
  data <- data[!is.na(data$value), ]
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
# Add fake pvals if missing for some descriptors
missing_pvals <- setdiff(descriptors, setdiff(colnames(pvals), c("type1", "type2", "bin1", "bin2")))
if (length(missing_pvals) > 0) {
  print(paste0("There are ", length(missing_pvals), " missing p-values"))
  for (desc in missing_pvals) {
    for (desc in missing_pvals) {
        pvals[, desc] <- 1
    }
  }
}
data <- pivot_longer(data, cols = all_of(descriptors), names_to = "feature", values_to = "value")
pvals <- pivot_longer(pvals, cols = all_of(descriptors), names_to = "feature", values_to = "p")
colnames(pvals) <- c("group1", "group2", "bin1", "bin2", "feature", "p")


##### Sequence length #####
data_len <- data[data$feature == "length", ]
data_len$value <- data_len$value / 3
data_len$type <- factor(data_len$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_len, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "Sequence length distribution (aa)",
       x = "Sequence type",
       y = "Length (aa)") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"),
                              "iorf" = paste0("iorf\n(n = ", n_iorf, ")"))) +
  scale_y_continuous(breaks = seq(0, 800, 200)) +
  stat_pvalue_manual(get_pvals("length", data_len, 0.08), label = "p.signif", inherit.aes = FALSE, hide.ns = TRUE, tip.length = 0.001) +
  annotate("text", x = 2.5, y = 1150,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
ggsave(paste0(out_folder, "/sequence_length.png"))



##### GC rate #####
data_gc <- data[data$feature == "gc_rate", ]
data_gc$type <- factor(data_gc$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_gc, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "GC ratio distribution",
       x = "Sequence type",
       y = "GC ratio: sequence GC % / genome GC %") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"),
                              "iorf" = paste0("iorf\n(n = ", n_iorf, ")"))) +
  scale_y_continuous(breaks = seq(0, 1.25, 0.25)) +
  stat_pvalue_manual(get_pvals("gc_rate", data_gc, 0.1), label = "p.signif", inherit.aes = FALSE, hide.ns = "p", tip.length = 0.005) +
  annotate("text", x = 2.5, y = 1.65,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
ggsave(paste0(out_folder, "/gc_content.png"))



##### GC rate (intergenic) #####
data_gc_inter <- data[data$feature == "inter_gc_rate", ]
data_gc_inter$type <- factor(data_gc_inter$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_gc_inter, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "GC ratio distribution (intergenic)",
       x = "Sequence type",
       y = "GC ratio: sequence GC % / intergenic ORFs GC %") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"),
                              "iorf" = paste0("iorf\n(n = ", n_iorf, ")"))) +
  scale_y_continuous(breaks = seq(0, 1.4, 0.2)) +
  stat_pvalue_manual(get_pvals("gc_rate", data_gc_inter, 0.1), label = "p.signif", inherit.aes = FALSE, hide.ns = "p", tip.length = 0.005) +
  annotate("text", x = 2.5, y = 1.8,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
ggsave(paste0(out_folder, "/gc_content_intergenic.png"))



##### Aromaticity #####
data_aro <- data[data$feature == "aromaticity", ]
data_aro$type <- factor(data_aro$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_aro, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "Aromaticity distribution",
       x = "Sequence type",
       y = "Aromaticity") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"),
                              "iorf" = paste0("iorf\n(n = ", n_iorf, ")"))) +
  scale_y_continuous(breaks = seq(0, 0.25, 0.05)) +
  stat_pvalue_manual(get_pvals("aromaticity", data_aro, 0.01), label = "p.signif", inherit.aes = FALSE, hide.ns = "p", tip.length = 0.005) + 
  annotate("text", x = 2.5, y = 0.3,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
ggsave(paste0(out_folder, "/aromaticity.png"))



##### Instability index #####
data_inst <- data[data$feature == "instability", ]
data_inst$type <- factor(data_inst$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_inst, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "Instability index distribution",
       x = "Sequence type",
       y = "Instability index") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"),
                              "iorf" = paste0("iorf\n(n = ", n_iorf, ")"))) +
  stat_pvalue_manual(get_pvals("instability", data_inst, 0.07), label = "p.signif", inherit.aes = FALSE, hide.ns = "p", tip.length = 0.005) +
  annotate("text", x = 2.5, y = 180,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black") +
  scale_y_continuous(breaks = seq(0, 125, 25))
ggsave(paste0(out_folder, "/instability_index.png"))



##### Flexibility #####
data_flex <- data[data$feature == "mean_flexibility", ]
data_flex$type <- factor(data_flex$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_flex, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "Mean flexibility distribution",
       x = "Sequence type",
       y = "Mean flexibility") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"),
                              "iorf" = paste0("iorf\n(n = ", n_iorf, ")"))) +
  stat_pvalue_manual(get_pvals("mean_flexibility", data_flex, 0.1), label = "p.signif", inherit.aes = FALSE, hide.ns = "p", tip.length = 0.005) +
  scale_y_continuous(breaks = seq(0.940, 1.040, 0.02)) +
  annotate("text", x = 2.5, y = 1.05,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
ggsave(paste0(out_folder, "/mean_flexibility.png"))



##### Hydrophobicity #####
data_hydro <- data[data$feature == "hydropathy", ]
data_hydro$type <- factor(data_hydro$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_hydro, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "Mean hydrophobicity distribution",
       x = "Sequence type",
       y = "Mean hydrophobicity") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"),
                              "iorf" = paste0("iorf\n(n = ", n_iorf, ")"))) +
  stat_pvalue_manual(get_pvals("hydropathy", data_hydro, 0.1), label = "p.signif", inherit.aes = FALSE, hide.ns = "p") +
  annotate("text", x = 2.5, y = 2.2,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
ggsave(paste0(out_folder, "/hydrophobicity.png"))



##### HCA #####
data_hca <- data[data$feature == "hca", ]
data_hca$type <- factor(data_hca$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_hca, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "HCA distribution",
       x = "Sequence type",
       y = "HCA") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"),
                              "iorf" = paste0("iorf\n(n = ", n_iorf, ")"))) +
  scale_y_continuous(breaks = seq(-10, 10, 5)) +
  stat_pvalue_manual(get_pvals("hca", data_hca, 0.08), label = "p.signif", inherit.aes = FALSE, hide.ns = "p", tip.length = 0.005) +
  annotate("text", x = 2.5, y = 18,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
ggsave(paste0(out_folder, "/hca.png"))




###### Intrinsic disorder (IUpred) ######
data_disord <- data[data$feature == "disord", ]
data_disord$type <- factor(data_disord$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_disord, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "Intrinsic disorder distribution (IUPred)",
       x = "Sequence type",
       y = "Intrinsic disorder") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"),
                              "iorf" = paste0("iorf\n(n = ", n_iorf, ")"))) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  stat_pvalue_manual(get_pvals("disord", data_disord, 0.08), label = "p.signif", inherit.aes = FALSE, hide.ns = "p", tip.length = 0.005) +
  annotate("text", x = 2.5, y = 1.5,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
ggsave(paste0(out_folder, "/intrinsic_disorder.png"))




###### Aggregation (tango) ######
data_agg <- data[data$feature == "aggreg", ]
data_agg$type <- factor(data_agg$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_agg, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "Aggregation distribution (Tango)",
       x = "Sequence type",
       y = "Aggregation") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"),
                              "iorf" = paste0("iorf\n(n = ", n_iorf, ")"))) +
  scale_y_continuous(breaks = seq(0, 0.9, 0.3)) +
  stat_pvalue_manual(get_pvals("aggreg", data_agg, 0.08), label = "p.signif", inherit.aes = FALSE, hide.ns = "p", tip.length = 0.005) +
  annotate("text", x = 2.5, y = 1.1,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
ggsave(paste0(out_folder, "/aggregation.png"))




##### AA use #####
### Polar ###
data_polar <- data[data$feature == "polar_use", ]
data_polar$type <- factor(data_polar$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_polar, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "Polar AA use distribution (S, T, N, Q)",
       x = "Sequence type",
       y = "Polar AA use") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"),
                              "iorf" = paste0("iorf\n(n = ", n_iorf, ")"))) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.25)) +
  stat_pvalue_manual(get_pvals("polar_use", data_polar, 0.08), label = "p.signif", inherit.aes = FALSE, hide.ns = "p", tip.length = 0.005) +
  annotate("text", x = 2.5, y = 0.68,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
ggsave(paste0(out_folder, "/aa_polar_use.png"))

### hydrophobic ###
data_hydro_use <- data[data$feature == "hydrophobic_use", ]
data_hydro_use$type <- factor(data_hydro_use$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_hydro_use, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "Hydrophobic AA use distribution (M, Y, V, L, I, F, W)",
       x = "Sequence type",
       y = "Hydrophobic AA use") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"),
                              "iorf" = paste0("iorf\n(n = ", n_iorf, ")"))) +
  stat_pvalue_manual(get_pvals("hydrophobic_use", data_hydro_use, 0.08), label = "p.signif", inherit.aes = FALSE, hide.ns = "p", tip.length = 0.005) +
  scale_y_continuous(breaks = seq(0, 0.75, 0.25)) +
  annotate("text", x = 2.5, y = 0.95,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
ggsave(paste0(out_folder, "/aa_hydrophobic_use.png"))

### positive ###
data_pos <- data[data$feature == "positive_use", ]
data_pos$type <- factor(data_pos$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_pos, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "Positive AA use distribution (K, R, H)",
       x = "Sequence type",
       y = "Positive AA use") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"),
                              "iorf" = paste0("iorf\n(n = ", n_iorf, ")"))) +
  scale_y_continuous(breaks = seq(0, 0.4, 0.1)) +
  stat_pvalue_manual(get_pvals("positive_use", data_pos, 0.08), label = "p.signif", inherit.aes = FALSE, hide.ns = "p", tip.length = 0.005) +
  annotate("text", x = 2.5, y = 0.48,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
ggsave(paste0(out_folder, "/aa_positive_use.png"))

### negative ###
data_neg <- data[data$feature == "negative_use", ]
data_neg$type <- factor(data_neg$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_neg, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "Negative AA use distribution (D, E)",
       x = "Sequence type",
       y = "Negative AA use") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"),
                              "iorf" = paste0("iorf\n(n = ", n_iorf, ")"))) +
  scale_y_continuous(breaks = seq(0, 0.3, 0.1)) +
  stat_pvalue_manual(get_pvals("negative_use", data_neg, 0.09), label = "p.signif", inherit.aes = FALSE, hide.ns = "p", tip.length = 0.005) +
  annotate("text", x = 2.5, y = 0.47,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
ggsave(paste0(out_folder, "/aa_negative_use.png"))

### Proline glycine ###
data_pro_gly <- data[data$feature == "proline.glycine_use", ]
data_pro_gly$type <- factor(data_pro_gly$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_pro_gly, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "Proline and glycine use distribution",
       x = "Sequence type",
       y = "Proline and glycine AA use") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"),
                              "iorf" = paste0("iorf\n(n = ", n_iorf, ")"))) +
  scale_y_continuous(breaks = seq(0, 0.3, 0.1)) +
  stat_pvalue_manual(get_pvals("proline.glycine_use", data_pro_gly, 0.1), label = "p.signif", inherit.aes = FALSE, hide.ns = "p") +
  annotate("text", x = 2.5, y = 0.3,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
ggsave(paste0(out_folder, "/aa_proline_glycine_use.png"))

### Cysteine ###
data_cys <- data[data$feature == "cysteine_use", ]
data_cys$type <- factor(data_cys$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_cys, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "Cysteine use distribution",
       x = "Sequence type",
       y = "Cysteine AA use") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"),
                              "iorf" = paste0("iorf\n(n = ", n_iorf, ")"))) +
  stat_pvalue_manual(get_pvals("cysteine_use", data_cys, 0.1), label = "p.signif", inherit.aes = FALSE, hide.ns = "p") +
  annotate("text", x = 2.5, y = 0.1,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
ggsave(paste0(out_folder, "/aa_cysteine_use.png"))

### Alanine ###
data_ala <- data[data$feature == "alanine_use", ]
data_ala$type <- factor(data_ala$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_ala, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "Alanine use distribution",
       x = "Sequence type",
       y = "Alanine AA use") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n = ", n_cds, ")"),
                              "trg" = paste0("trg\n(n = ", n_trg, ")"),
                              "denovo" = paste0("denovo\n(n = ", n_denovo, ")"),
                              "iorf" = paste0("iorf\n(n = ", n_iorf, ")"))) +
  scale_y_continuous(breaks = seq(0, 0.2, 0.05)) +
  stat_pvalue_manual(get_pvals("alanine_use", data_ala, 0.08), label = "p.signif", inherit.aes = FALSE, hide.ns = "p", tip.length = 0.005) +
  annotate("text", x = 2.5, y = 0.21,
           label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05", hjust = 1, vjust = 1, size = 3, color = "black")
ggsave(paste0(out_folder, "/aa_alanine_use.png"))



##### GC % ######
# Get the summary of the number of genomes
n_gen_cds <- nrow(unique(data[data$type == "cds", "genome"]))
n_gen_trg <- nrow(unique(data[data$type == "trg", "genome"]))
n_gen_denovo <- nrow(unique(data[data$type == "denovo", "genome"]))
n_gen_iorf <- nrow(unique(data[data$type == "iorf", "genome"]))

### Species gc ###
data_gc_species <- data[data$feature == "gc_species", ]
# Keep only one row per unique couple genome and type
data_gc_species <- data_gc_species[!duplicated(data_gc_species[, c("genome", "type")]), ]
data_gc_species$type <- factor(data_gc_species$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_gc_species, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "GC % distribution (species)",
       x = "Sequence type",
       y = "GC %") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n genomes = ", n_gen_cds, ")"),
                              "trg" = paste0("trg\n(n genomes = ", n_gen_trg, ")"),
                              "denovo" = paste0("denovo\n(n genomes = ", n_gen_denovo, ")"),
                              "iorf" = paste0("iorf\n(n genomes = ", n_gen_iorf, ")")))
ggsave(paste0(out_folder, "/gc_species.png"))


### Intergenic gc ###
data_gc_inter <- data[data$feature == "inter_gc_species", ]
data_gc_inter <- data_gc_inter[!duplicated(data_gc_inter[, c("genome", "type")]), ]
data_gc_inter$type <- factor(data_gc_inter$type, levels = c("cds", "trg", "denovo", "iorf"))
ggplot(data_gc_inter, aes(x = type, y = value, fill = type)) +
  geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
  labs(title = "GC % distribution (intergenic)",
       x = "Sequence type",
       y = "GC %") +
  scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("cds" = paste0("cds\n(n genomes = ", n_gen_cds, ")"),
                              "trg" = paste0("trg\n(n genomes = ", n_gen_trg, ")"),
                              "denovo" = paste0("denovo\n(n genomes = ", n_gen_denovo, ")"),
                              "iorf" = paste0("iorf\n(n genomes = ", n_gen_iorf, ")")))
ggsave(paste0(out_folder, "/gc_species_intergenic.png"))