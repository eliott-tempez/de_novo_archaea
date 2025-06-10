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


home <- FALSE
use_violins <- FALSE
print_pval_labels <- FALSE


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

  pos <- max(as.numeric(y_mat[, "whisker_top"])) + reference_val * (fact + 0.05)
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


# Remove the values that would be outliers from data
remove_outliers <- function(data) {
  # Remove outliers
  data <- data %>%
    filter(!is.na(value)) %>%
    group_by(type) %>%
    mutate(value = ifelse(value > quantile(value, 0.75) + IQR(value) * 1.5, NA, value)) %>%
    filter(!is.na(value)) %>%
    mutate(value = ifelse(value < quantile(value, 0.25) - IQR(value) * 1.5, NA, value)) %>%
    filter(!is.na(value)) %>%
    ungroup()
  return(data)
}


# Plot
plot_data <- function(data_local,
                      title, y_axis,
                      y_scale,
                      factor, tip_length,
                      pval_ypos = NA) {

  boxplot_width <- ifelse(use_violins, 0.2, 0.8)
  boxplot_alpha <- ifelse(use_violins, 0.7, 1)
  pvals_local <- get_pvals(data_local$feature[1], data_local, factor)
  plot_pvals <- !all(is.na(pvals_local$y.position))

  # Remove outliers for violin plot
  data_sans_outliers <- remove_outliers(data_local)

  # Plot
  p <- ggplot(data_local, aes(x = type, y = value, fill = type))

  if (use_violins) {
    p <- p + geom_violin(data = data_sans_outliers, na.rm = TRUE, colour = "#2c2c2c", scale = "width", alpha = 0.7)
  }

  p <- p + geom_boxplot(data = data_local, na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE, width = boxplot_width, alpha = boxplot_alpha) +
    labs(title = title,
         x = "Sequence type",
         y = y_axis) +
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
    scale_y_continuous(breaks = y_scale)

  if (plot_pvals) {
    p <- p + stat_pvalue_manual(pvals_local, label = "p.signif", inherit.aes = FALSE, hide.ns = TRUE, tip.length = tip_length, size = 5)
  }

  if (print_pval_labels) {
    p <- p + annotate("text", x = 3, y = pval_ypos,
                      label = "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05",
                      hjust = 1, vjust = 1, size = 4, color = "black")
  }

return(p)
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











##################### PLOTS #####################

##### Sequence length #####
data_len <- data[data$feature == "length", ]
data_len$value <- data_len$value / 3
data_len$type <- factor(data_len$type, levels = c("cds", "trg", "denovo", "iorf"))

p <- plot_data(data_len,
               title = "Sequence length distribution (aa)",
               y_axis = "Length (aa)",
               y_scale = seq(0, 800, 200),
               factor = 0.05,
               tip_length = 0.002)
p

if (!use_violins) {
  ggsave(paste0(out_folder, "/sequence_length.png"))
} else {
  ggsave(paste0(out_folder, "/sequence_length_violin.png"))
}


##### GC rate #####
data_gc <- data[data$feature == "gc_rate", ]
data_gc$type <- factor(data_gc$type, levels = c("cds", "trg", "denovo", "iorf"))

p <- plot_data(data_gc,
          title = "GC ratio distribution",
          y_axis = "GC ratio: sequence GC % / species GC %",
          y_scale = seq(0, 1.4, 0.2),
          factor = 0.05,
          tip_length = 0.005,
          pval_ypos = 1.8)
p

if (!use_violins) {
  ggsave(paste0(out_folder, "/gc_content.png"))
} else {
  ggsave(paste0(out_folder, "/gc_content_violin.png"))
}



##### GC rate (intergenic) #####
data_gc_inter <- data[data$feature == "inter_gc_rate", ]
data_gc_inter$type <- factor(data_gc_inter$type, levels = c("cds", "trg", "denovo", "iorf"))

p <- plot_data(data_gc_inter,
          title = "GC ratio distribution (intergenic)",
          y_axis = "GC ratio: sequence GC % / intergenic ORFs GC %",
          y_scale = seq(0, 1.4, 0.2),
          factor = 0.05,
          tip_length = 0.005,
          pval_ypos = 1.8)
p

if (!use_violins) {
  ggsave(paste0(out_folder, "/gc_content_intergenic.png"))
} else {
  ggsave(paste0(out_folder, "/gc_content_intergenic_violin.png"))
}



##### Aromaticity #####
data_aro <- data[data$feature == "aromaticity", ]
data_aro$type <- factor(data_aro$type, levels = c("cds", "trg", "denovo", "iorf"))

p <- plot_data(data_aro,
          title = "Aromaticity distribution",
          y_axis = "Aromaticity",
          y_scale = seq(0, 0.25, 0.05),
          factor = 0.05,
          tip_length = 0.005,
          pval_ypos = 0.3)
p

if (!use_violins) {
  ggsave(paste0(out_folder, "/aromaticity.png"))
} else {
  ggsave(paste0(out_folder, "/aromaticity_violin.png"))
}



##### Instability index #####
data_inst <- data[data$feature == "instability", ]
data_inst$type <- factor(data_inst$type, levels = c("cds", "trg", "denovo", "iorf"))

p <- plot_data(data_inst,
          title = "Instability index distribution",
          y_axis = "Instability index",
          y_scale = seq(0, 125, 25),
          factor = 0.05,
          tip_length = 0.005,
          pval_ypos = 180)
p

if (!use_violins) {
  ggsave(paste0(out_folder, "/instability_index.png"))
} else {
  ggsave(paste0(out_folder, "/instability_index_violin.png"))
}



##### Flexibility #####
data_flex <- data[data$feature == "mean_flexibility", ]
data_flex$type <- factor(data_flex$type, levels = c("cds", "trg", "denovo", "iorf"))

p <- plot_data(data_flex,
          title = "Mean flexibility distribution",
          y_axis = "Mean flexibility",
          y_scale = seq(0.940, 1.040, 0.02),
          factor = 0.05,
          tip_length = 0.005,
          pval_ypos = 1.05)
p

if (!use_violins) {
  ggsave(paste0(out_folder, "/mean_flexibility.png"))
} else {
  ggsave(paste0(out_folder, "/mean_flexibility_violin.png"))
}



##### Hydrophobicity #####
data_hydro <- data[data$feature == "hydropathy", ]
data_hydro$type <- factor(data_hydro$type, levels = c("cds", "trg", "denovo", "iorf"))

p <- plot_data(data_hydro,
          title = "Mean hydrophobicity distribution",
          y_axis = "Mean hydrophobicity",
          y_scale = seq(-2, 2.5, 0.5),
          factor = 0.1,
          tip_length = 0.005,
          pval_ypos = 2.2)
p

if (!use_violins) {
  ggsave(paste0(out_folder, "/hydrophobicity.png"))
} else {
  ggsave(paste0(out_folder, "/hydrophobicity_violin.png"))
}



##### HCA #####
data_hca <- data[data$feature == "hca", ]
data_hca$type <- factor(data_hca$type, levels = c("cds", "trg", "denovo", "iorf"))

p <- plot_data(data_hca,
          title = "HCA distribution",
          y_axis = "HCA",
          y_scale = seq(-10, 20, 5),
          factor = 0.05,
          tip_length = 0.005,
          pval_ypos = 18)
p

if (!use_violins) {
  ggsave(paste0(out_folder, "/hca.png"))
} else {
  ggsave(paste0(out_folder, "/hca_violin.png"))
}




###### Intrinsic disorder (IUpred) ######
data_disord <- data[data$feature == "disord", ]
data_disord$type <- factor(data_disord$type, levels = c("cds", "trg", "denovo", "iorf"))

p <- plot_data(data_disord,
          title = "Intrinsic disorder distribution (IUPred)",
          y_axis = "Intrinsic disorder",
          y_scale = seq(0, 1, 0.25),
          factor = 0.05,
          tip_length = 0.005,
          pval_ypos = 1.5)
p

if (!use_violins) {
  ggsave(paste0(out_folder, "/intrinsic_disorder.png"))
} else {
  ggsave(paste0(out_folder, "/intrinsic_disorder_violin.png"))
}




###### Aggregation (tango) ######
data_agg <- data[data$feature == "aggreg", ]
data_agg$type <- factor(data_agg$type, levels = c("cds", "trg", "denovo", "iorf"))

p <- plot_data(data_agg,
          title = "Aggregation distribution (Tango)",
          y_axis = "Aggregation",
          y_scale = seq(0, 0.9, 0.3),
          factor = 0.05,
          tip_length = 0.005,
          pval_ypos = 1.1)
p

if (!use_violins) {
  ggsave(paste0(out_folder, "/aggregation.png"))
} else {
  ggsave(paste0(out_folder, "/aggregation_violin.png"))
}




##### AA use #####
### Polar ###
data_polar <- data[data$feature == "polar_use", ]
data_polar$type <- factor(data_polar$type, levels = c("cds", "trg", "denovo", "iorf"))

p <- plot_data(data_polar,
          title = "Polar AA use distribution (S, T, N, Q)",
          y_axis = "Polar AA use",
          y_scale = seq(0, 0.5, 0.25),
          factor = 0.05,
          tip_length = 0.005,
          pval_ypos = 0.68)
p

if (!use_violins) {
  ggsave(paste0(out_folder, "/aa_polar_use.png"))
} else {
  ggsave(paste0(out_folder, "/aa_polar_use_violin.png"))
}

### hydrophobic ###
data_hydro_use <- data[data$feature == "hydrophobic_use", ]
data_hydro_use$type <- factor(data_hydro_use$type, levels = c("cds", "trg", "denovo", "iorf"))

p <- plot_data(data_hydro_use,
  title = "Hydrophobic AA use distribution (M, Y, V, L, I, F, W)",
  y_axis = "Hydrophobic AA use",
  y_scale = seq(0, 0.75, 0.25),
  factor = 0.08,
  tip_length = 0.005,
  pval_ypos = 0.95
)
p

if (!use_violins) {
  ggsave(paste0(out_folder, "/aa_hydrophobic_use.png"))
} else {
  ggsave(paste0(out_folder, "/aa_hydrophobic_use_violin.png"))
}

### positive ###
data_pos <- data[data$feature == "positive_use", ]
data_pos$type <- factor(data_pos$type, levels = c("cds", "trg", "denovo", "iorf"))

p <- plot_data(data_pos,
  title = "Positive AA use distribution (K, R, H)",
  y_axis = "Positive AA use",
  y_scale = seq(0, 0.4, 0.1),
  factor = 0.05,
  tip_length = 0.005,
  pval_ypos = 0.48
)
p

if (!use_violins) {
  ggsave(paste0(out_folder, "/aa_positive_use.png"))
} else {
  ggsave(paste0(out_folder, "/aa_positive_use_violin.png"))
}

### negative ###
data_neg <- data[data$feature == "negative_use", ]
data_neg$type <- factor(data_neg$type, levels = c("cds", "trg", "denovo", "iorf"))

p <- plot_data(data_neg,
  title = "Negative AA use distribution (D, E)",
  y_axis = "Negative AA use",
  y_scale = seq(0, 0.3, 0.1),
  factor = 0.05,
  tip_length = 0.005,
  pval_ypos = 0.47
)
p

if (!use_violins) {
  ggsave(paste0(out_folder, "/aa_negative_use.png"))
} else {
  ggsave(paste0(out_folder, "/aa_negative_use_violin.png"))
}

### Proline glycine ###
data_pro_gly <- data[data$feature == "proline.glycine_use", ]
data_pro_gly$type <- factor(data_pro_gly$type, levels = c("cds", "trg", "denovo", "iorf"))

p <- plot_data(data_pro_gly,
  title = "Proline and glycine use distribution",
  y_axis = "Proline and glycine AA use",
  y_scale = seq(0, 0.3, 0.1),
  factor = 0.05,
  tip_length = 0.005,
  pval_ypos = 0.3
)
p

if (!use_violins) {
  ggsave(paste0(out_folder, "/aa_proline_glycine_use.png"))
} else {
  ggsave(paste0(out_folder, "/aa_proline_glycine_use_violin.png"))
}

### Cysteine ###
data_cys <- data[data$feature == "cysteine_use", ]
data_cys$type <- factor(data_cys$type, levels = c("cds", "trg", "denovo", "iorf"))

p <- plot_data(data_cys,
  title = "Cysteine use distribution",
  y_axis = "Cysteine AA use",
  y_scale = seq(0, 0.2, 0.05),
  factor = 0.05,
  tip_length = 0.005,
  pval_ypos = 0.21
)
p

if (!use_violins) {
  ggsave(paste0(out_folder, "/aa_cysteine_use.png"))
} else {
  ggsave(paste0(out_folder, "/aa_cysteine_use_violin.png"))
}

### Alanine ###
data_ala <- data[data$feature == "alanine_use", ]
data_ala$type <- factor(data_ala$type, levels = c("cds", "trg", "denovo", "iorf"))

p <- plot_data(data_ala,
  title = "Alanine use distribution",
  y_axis = "Alanine AA use",
  y_scale = seq(0, 0.2, 0.05),
  factor = 0.05,
  tip_length = 0.005,
  pval_ypos = 0.21
)
p

if (!use_violins) {
  ggsave(paste0(out_folder, "/aa_alanine_use.png"))
} else {
  ggsave(paste0(out_folder, "/aa_alanine_use_violin.png"))
}



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