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
n_bins <- 2
plot_pvals <- TRUE
save_plots <- TRUE


# Files
input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/17_compare_denovo_sequences/sequence_features_good_candidates_all.csv"
out_folder <- paste0("/home/eliott.tempez/Documents/M2_Stage_I2BC/results/17_compare_denovo_sequences/", n_bins, "_bins/")
pval_file <- paste0(out_folder, "pvalues_", n_bins, "_bins.tsv")
bins_file <- paste0(out_folder, "bin_indexes_", n_bins, ".csv")
#input_file <- "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/results/17_compare_denovo_sequences/features_good_candidates_all.csv"
#pval_file <- "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/results/17_compare_denovo_sequences/pvalues_2_bins.tsv"
#bins_file <- "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/results/17_compare_denovo_sequences/sequences/2_bins/bin_indexes_2.csv"
#out_folder <- "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/results/17_compare_denovo_sequences/sequences/2_bins"


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
if (plot_pvals) {
  signif_label <- "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05"} else {
  signif_label <- ""}





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
  min_whisker_base <- max(data$value, na.rm = TRUE)
  y_mat <- matrix(NA, nrow = 6, ncol = 3)
  colnames(y_mat) <- c("type", "bin", "whisker_top")
  i <- 0
  for (type in types) {
    for (bin in c("0", "1")) {
      i <- i + 1
      all_val <- as.numeric(as.vector(data[data$type == type & data$bin == bin, "value"])$value)
      # Remove NA values
      all_val <- all_val[!is.na(all_val)]
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


# Get the x-axis labels
get_x_labels <- function(data) {
  x_labels <- c()
  groups <- unique(data$group)
  for (i in seq_along(groups)) {
    group <- as.character(levels(groups)[i])
    if (grepl("trg", group)) {
      n <- as.numeric(strsplit(group, "_")[[1]][1])
      x_labels <- c(x_labels, bin_labels[as.numeric(strsplit(group, "_")[[1]][1]) + 1])
    } else {
      x_labels <- c(x_labels, "")
    }
  }
  return(x_labels)
}



# Get the plot
get_plot <- function(data, data_len_summary, feature, title, print_pval = c(NA), scale_y = c(NA), n_y_pos = NA) {
  labels_x_scale <- get_x_labels(data)

  p <- ggplot(data, aes(x = group, y = value, fill = type)) +
    geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE) +
    labs(title = "",
         x = "% GC (whole genome)",
         y = title) +
    scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c")) +
    theme_minimal() +
    scale_x_discrete(labels = labels_x_scale) +
    theme(axis.text.x = element_text(size = 16),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 12),
          legend.title = element_blank())

  if (grepl("_use", feature)) {
    aa <- strsplit(feature, "_")[[1]][1]
    p <- p + labs(title = aa)
    p <- p +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
            axis.text.y = element_text(size = 16),
            legend.position = "none",
            axis.text.x = element_text(size = 8),
            axis.title.x = element_text(size = 0))
    p <- p +
      theme(panel.border = element_rect(colour = "#ada9a9", 
                                        fill = NA, linewidth = 1))
  }


  if (!all(is.na(scale_y))) {
    p <- p + scale_y_continuous(breaks = scale_y)
  }

  if (!is.na(n_y_pos)) {
    p <- p +
      geom_text(data = data_len_summary,
                aes(x = group,
                    y = n_y_pos,
                    label = paste0("n = ", n),
                    group = type),
                position = position_dodge(width = 0.75), 
                vjust = -0.5, size = 4)
  }

  if (!all(is.na(print_pval))) {
    pval_factor <- print_pval[1]
    only_ns <- print_pval[2]
    y_annotation <- print_pval[3]

    pvals_local <- get_pvals(feature, data, pval_factor, only_ns)

    if (!all(is.na(pvals_local$y.position))) {
      p <- p +
        stat_pvalue_manual(pvals_local,
                           label = "p.signif",
                           inherit.aes = FALSE,
                           hide.ns = !only_ns,
                           tip.length = 0.01) +
        annotate("text", x = 4, y = y_annotation,
                 label = signif_label, hjust = 1, vjust = 1,
                 size = 3, color = "black")
    }
  }
  return(p)
}




###################################################################
############################### MAIN ##############################
###################################################################

###### Sequence length ######
data_len <- add_dummy_rows(data, "length", n_bins)
data_len_summary <- get_ncds_conditions(data_len, n_bins)

# Pvals
if (plot_pvals) {
  pval_factor <- 0.1
  only_ns <- TRUE
  y_annotation <- 3000
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}

# Plot
p <- get_plot(data_len,
              data_len_summary,
              "length",
              "Sequence length distribution",
              n_y_pos = -50,
              print_pval = pval_vect,
              scale_y = c(50, 100, 500, 1000, 2000))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/sequence_length.png"), plot = p)
}



###### GC ratio ######
data_gc <- add_dummy_rows(data, "gc_rate", n_bins)
data_gc_summary <- get_ncds_conditions(data_gc, n_bins)

# Pvals
if (plot_pvals) {
  pval_factor <- 0.08
  only_ns <- FALSE
  y_annotation <- 1.5
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}

# Plot
p <- get_plot(data_gc,
              data_gc_summary,
              "gc_rate",
              "GC ratio: sequence GC % / genome GC %",
              n_y_pos = 0.53,
              print_pval = pval_vect,
              scale_y = seq(0, 1.3, 0.2))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/gc_content.png"), plot = p)
}



###### GC ratio (intergenic) ######
data_gc_inter <- add_dummy_rows(data, "inter_gc_rate", n_bins)
data_gc_inter_summary <- get_ncds_conditions(data_gc_inter, n_bins)

# Pvals
if (plot_pvals) {
  pval_factor <- 0.06
  only_ns <- FALSE
  y_annotation <- 1.8
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}

# Plot
p <- get_plot(data_gc_inter,
              data_gc_inter_summary,
              "inter_gc_rate",
              "GC ratio: sequence GC % / iORF GC %",
              n_y_pos = 0.62,
              print_pval = pval_vect,
              scale_y = seq(0, 1.6, 0.2))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/gc_content_intergenic.png"), plot = p)
}



###### Aromaticity ######
data_aro <- add_dummy_rows(data, "aromaticity", n_bins)
data_aro_summary <- get_ncds_conditions(data_aro, n_bins)

# Pvals
if (plot_pvals) {
  pval_factor <- 0.1
  only_ns <- FALSE
  y_annotation <- 0.27
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}

# Plot
p <- get_plot(data_aro,
              data_aro_summary,
              "aromaticity",
              "Aromaticity",
              n_y_pos = -0.02,
              print_pval = pval_vect,
              scale_y = seq(0, 0.25, 0.05))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/aromaticity.png"), plot = p)
}




###### Instability index ######
data_inst <- add_dummy_rows(data, "instability", n_bins)
data_inst_summary <- get_ncds_conditions(data_inst, n_bins)

# Pvals
if (plot_pvals) {
  pval_factor <- 0.06
  only_ns <- FALSE
  y_annotation <- 120
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}

# Plot
p <- get_plot(data_inst,
              data_inst_summary,
              "instability",
              "Instability index",
              n_y_pos = -15,
              print_pval = pval_vect,
              scale_y = seq(0, 75, 25))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/instability_index.png"), plot = p)
}



###### Flexibility #######
data_flex <- add_dummy_rows(data, "mean_flexibility", n_bins)
data_flex_summary <- get_ncds_conditions(data_flex, n_bins)

# Pvals
if (plot_pvals) {
  pval_factor <- 0.05
  only_ns <- FALSE
  y_annotation <- 1.05
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}

# Plot
p <- get_plot(data_flex,
              data_flex_summary,
              "mean_flexibility",
              "Mean flexibility",
              n_y_pos = 0.95,
              print_pval = pval_vect,
              scale_y = seq(0.96, 1.04, 0.02))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/mean_flexibility.png"), plot = p)
}



###### Hydrophobicity ######
data_hyd <- add_dummy_rows(data, "hydropathy", n_bins)
data_hyd_summary <- get_ncds_conditions(data_hyd, n_bins)

# Pvals
if (plot_pvals) {
  pval_factor <- 0.05
  only_ns <- FALSE
  y_annotation <- 3
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}

# Plot
p <- get_plot(data_hyd,
              data_hyd_summary,
              "hydropathy",
              "Hydrophobicity",
              n_y_pos = -2,
              print_pval = pval_vect,
              scale_y = seq(-2, 2, 1))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/hydrophobicity.png"), plot = p)
}



###### HCA ######
data_hca <- add_dummy_rows(data, "hca", n_bins)
data_hca_summary <- get_ncds_conditions(data_hca, n_bins)

# Pvals
if (plot_pvals) {
  pval_factor <- 0.1
  only_ns <- TRUE
  y_annotation <- 8
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}

# Plot
p <- get_plot(data_hca,
              data_hca_summary,
              "hca",
              "HCA",
              n_y_pos = -11,
              print_pval = pval_vect,
              scale_y = seq(-10, 2, 4))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/hca.png"), plot = p)
}




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
    # Get the data
    data_aa <- add_dummy_rows(data, paste0(aa, "_use"), n_bins)
    data_aa_summary <- get_ncds_conditions(data_aa, n_bins)

    # Pvals
    if (plot_pvals) {
      pval_factor <- 0.2
      only_ns <- FALSE
      y_annotation <- NA
      pval_vect <- c(pval_factor, only_ns, y_annotation)
    } else {
      pval_vect <- c(NA, NA, NA)
    }

    # Plot
    p <- get_plot(data_aa,
                  data_aa_summary,
                  paste0(aa, "_use"),
                  "",
                  print_pval = pval_vect)


    aa_plots <- c(aa_plots, list(p))
  }

  fig <- ggarrange(plotlist = aa_plots, common.legend = TRUE, legend = "right")
  fig <- annotate_figure(fig, bottom = text_grob("% GC (whole genome)\n", size = 14),
                  left = text_grob("% use", rot = 90, size = 14),
                  top = text_grob(paste("Amino-acid distribution:",
                                        aa_type_name, "\n"),
                                  size = 18),
                  right = text_grob(paste0("\n", signif_label), rot = 90, size = 10))
  print(fig)

  # Save the plot
  if (save_plots) {
    ggsave(paste0(out_folder, "/aa_", aa_type_name, "_use.png"), plot = fig)
  }
}


