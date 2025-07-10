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
only_denovo_genomes <- FALSE
use_violins <- TRUE
ONE_SIDED <- TRUE


# Files

input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/7_characterise_denovo/sequence_features_good_candidates_all.csv"
good_candidates_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/4_integrity_analysis/good_candidates.txt"
out_folder <- paste0("/home/eliott.tempez/Documents/M2_Stage_I2BC/results/8_plot_denovo_characterisation/", n_bins, "_bins/")
pval_file <- paste0(out_folder, "pvalues_", n_bins, "_bins.tsv")
if (ONE_SIDED) {
  pval_file <- paste0(out_folder, "pvalues_", n_bins, "_bins_one_sided.csv")
}
bins_file <- paste0(out_folder, "bin_indexes_", n_bins, ".csv")



# Read data
data <- read.table(input_file, header = TRUE, sep = "\t")
if (plot_pvals) {
  pvals <- read.table(pval_file, header = TRUE, sep = "\t")
}
if (only_denovo_genomes) {
  data <- data[data$genome %in% unique(data[data$type == "denovo", "genome"]), ]
}
# Sample 5000 sequences for cds, trg and iorf
n_cds <- nrow(data[data$type == "cds", ])
n_trg <- nrow(data[data$type == "trg", ])
n_denovo <- nrow(data[data$type == "denovo", ])
n_iorf <- nrow(data[data$type == "iorf", ])
if (n_cds > 5000) {
  data <- data[data$type != "cds" |
                 (data$type == "cds" &
                    seq_len(nrow(data[data$type == "cds", ])) %in% sample(n_cds, 5000)), ]
}
if (n_trg > 5000) {
  data <- data[data$type != "trg" |
                 (data$type == "trg" &
                    seq_len(nrow(data[data$type == "trg", ])) %in% sample(n_trg, 5000)), ]
}
if (n_iorf > 5000) {
  data <- data[data$type != "iorf" |
                 (data$type == "iorf" &
                    seq_len(nrow(data[data$type == "iorf", ])) %in% sample(n_iorf, 5000)), ]
}
# Pivot to longer format
descriptors <- setdiff(colnames(data), c("genome", "cds", "type"))
not_present <- setdiff(descriptors, colnames(pvals))
descriptors <- setdiff(descriptors, not_present)
data <- pivot_longer(data, cols = all_of(descriptors), names_to = "feature", values_to = "value")
if (plot_pvals) {
  pvals <- pivot_longer(pvals, cols = all_of(setdiff(descriptors, c("gc_species", "inter_gc_species"))), names_to = "feature", values_to = "pval")
  pvals$p <- pvals$pval
}

## Add the bins to the data
bins <- read.table(bins_file, header = TRUE, sep = "")
bins <- bins %>%
  left_join((data[c("cds", "genome")]), by = "cds")
if (length(unique(na.omit(bins$genome))) != 116 && !only_denovo_genomes) {
  stop("The bins file should contain all 116 genomes.")
}
# Get the bins per genome
bins <- bins %>%
  select(c("genome", "bin")) %>%
  distinct()
data <- data %>%
  left_join(bins, by = c("genome"))
data$bin <- as.character(data$bin)
# GC %
min_gc <- 39.4
max_gc <- 56.4
bin_limits <- seq(min_gc, max_gc, length.out = n_bins + 1)
bin_labels <- paste0(round(bin_limits[-length(bin_limits)], 1), " - ", round(bin_limits[-1], 1), " %")
if (plot_pvals) {
  signif_label <- "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05"} else {
  signif_label <- ""}
# Change the de novo labels to good candidates only
good_candidates <- scan(good_candidates_file, what = "", sep = "\n")
data[data$type == "denovo", "type"] <- "trg"
data[data$cds %in% good_candidates, "type"] <- "denovo"




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
  data_levels <- c("0_cds", "0_trg", "0_denovo", "0_iorf")
  if (n_dummy > 0) {
    for (i in 1:n_dummy) {
      blank <- paste0(i - 1, "_blank")
      cds <- paste0(i, "_cds")
      trg <- paste0(i, "_trg")
      denovo <- paste0(i, "_denovo")
      iorf <- paste0(i, "_iorf")
      data_levels <- c(data_levels, blank, cds, trg, denovo, iorf)
    }
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
iorf_0 <- ((local_pvals$type1 == "iorf" & local_pvals$bin1 == "0") |
             (local_pvals$type2 == "iorf" & local_pvals$bin2 == "0"))
iorf_1 <- ((local_pvals$type1 == "iorf" & local_pvals$bin1 == "1") |
             (local_pvals$type2 == "iorf" & local_pvals$bin2 == "1"))
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
    iorf_0 & trg_1 ~ 18,
    iorf_1 & trg_0 ~ 19,
    iorf_0 & denovo_1 ~ 20,
    iorf_1 & denovo_0 ~ 21,
    iorf_0 & cds_1 ~ 22,
    iorf_1 & cds_0 ~ 23,
    iorf_0 & trg_0 ~ 24,
    iorf_1 & trg_1 ~ 25,
    iorf_0 & cds_0 ~ 26,
    iorf_1 & cds_1 ~ 27,
    iorf_0 & denovo_0 ~ 28,
    iorf_1 & denovo_1 ~ 29,
    iorf_0 & iorf_1 ~ 30,
    TRUE ~ 31
  )) %>%
  arrange(order)

  # Get the top of the whisker and the p-val for each type
  types <- c("denovo", "trg", "cds", "iorf")
  whisker_tops <- c()
  print(head(data))
  min_whisker_base <- max(data$value, na.rm = TRUE)
  y_mat <- matrix(NA, nrow = 8, ncol = 3)
  colnames(y_mat) <- c("type", "bin", "whisker_top")
  i <- 0
  for (type in types) {
    for (bin in as.character(unique(data$bin))) {
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


# Get the plot
get_plot <- function(data_local, data_len_summary, feature, title, print_pval = c(NA), scale_y = c(NA), n_y_pos = NA) {
  labels_x_scale <- get_x_labels(data_local)

  factor <- 0.05
  boxplot_width <- ifelse(use_violins, 0.2, 0.8)
  boxplot_alpha <- ifelse(use_violins, 0.8, 1)
  pvals_local <- get_pvals(data_local$feature[1], data_local, factor)
  plot_pvals <- !all(is.na(pvals_local$y.position))

  # Remove outliers for violin plot
  data_sans_outliers <- add_dummy_rows(remove_outliers(data_local), feature, n_bins)

  p <- ggplot(data_local, aes(x = group, y = value, fill = type))

  if (use_violins) {
    p <- p + geom_violin(data = data_sans_outliers, na.rm = TRUE, colour = "#2c2c2c", scale = "width", alpha = 0.8)
  }


  p <- p + geom_boxplot(data = data_local, na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE, width = boxplot_width, alpha = boxplot_alpha) +
    labs(title = "",
         x = "% GC (whole genome)",
         y = title) +
    scale_fill_manual(values = c("#cc7f0a", "#ad4646", "#4d4c4c", "#693fb6"),
                      breaks = c("cds", "trg", "denovo", "iorf")) +
    theme_minimal() +
    scale_x_discrete(labels = labels_x_scale) +
    theme(axis.text.x = element_text(size = 22),
          axis.title.x = element_text(size = 24),
          axis.title.y = element_text(size = 24),
          axis.text.y = element_text(size = 22),
          legend.text = element_text(size = 20),)


  if (!all(is.na(scale_y))) {
    p <- p + scale_y_continuous(breaks = scale_y)
  }

  if (!all(is.na(print_pval))) {
    pval_factor <- print_pval[1]
    only_ns <- print_pval[2]
    y_annotation <- print_pval[3]

    pvals_local <- get_pvals(feature, data_local, pval_factor, only_ns)

    if (!all(is.na(pvals_local$y.position))) {
      p <- p +
        stat_pvalue_manual(pvals_local,
                           label = "p.signif",
                           inherit.aes = FALSE,
                           hide.ns = !only_ns,
                           tip.length = 0,
                           size = 10)

    }
  }
  return(p)
}




###################################################################
############################### MAIN ##############################
###################################################################

###### Sequence length ######
data_len <- add_dummy_rows(data, "length", n_bins)
data_len$type <- factor(data_len$type, levels = c("cds", "trg", "denovo", "iorf"))
# Use aa length
data_len$value <- data_len$value / 3
data_len_summary <- get_ncds_conditions(data_len, n_bins)

# Pvals
if (plot_pvals) {
  pval_factor <- 0.08
  only_ns <- FALSE
  y_annotation <- 3000
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}

# Plot
p <- get_plot(data_len,
              data_len_summary,
              "length",
              "Sequence length (aa)",
              n_y_pos = -50,
              print_pval = pval_vect,
              scale_y = seq(0, 800, 200))
p

# Save the plot
file_name <- paste0(out_folder, "sequence_length")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)



###### GC ratio ######
data_gc <- add_dummy_rows(data, "gc_rate", n_bins)
data_gc$type <- factor(data_gc$type, levels = c("cds", "trg", "denovo", "iorf"))
data_gc_summary <- get_ncds_conditions(data_gc, n_bins)

# Pvals
if (plot_pvals) {
  pval_factor <- 0.07
  only_ns <- FALSE
  y_annotation <- 1.8
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}

# Plot
p <- get_plot(data_gc,
              data_gc_summary,
              "gc_rate",
              "GC ratio: sequence GC % / genome GC %",
              n_y_pos = 0.4,
              print_pval = pval_vect,
              scale_y = seq(0, 1.4, 0.2))
p

# Save the plot
file_name <- paste0(out_folder, "gc_content")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)



###### GC ratio (intergenic) ######
data_gc_inter <- add_dummy_rows(data, "inter_gc_rate", n_bins)
data_gc_inter$type <- factor(data_gc_inter$type, levels = c("cds", "trg", "denovo", "iorf"))
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
              "GC ratio: sequence GC % / intergenic GC %",
              n_y_pos = 0.45,
              print_pval = pval_vect,
              scale_y = seq(0, 1.4, 0.2))
p

# Save the plot
file_name <- paste0(out_folder, "gc_content_intergenic")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)


###### GC ratio (iorfs) ######
data_gc_iorfs <- add_dummy_rows(data, "iorfs_gc_rate", n_bins)
data_gc_iorfs$type <- factor(data_gc_iorfs$type, levels = c("cds", "trg", "denovo", "iorf"))
data_gc_iorfs_summary <- get_ncds_conditions(data_gc_iorfs, n_bins)

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
p <- get_plot(data_gc_iorfs,
              data_gc_iorfs_summary,
              "iorfs_gc_rate",
              "GC ratio: sequence GC % / iorfs GC %",
              n_y_pos = 0.45,
              print_pval = pval_vect,
              scale_y = seq(0, 1.4, 0.2))
p

# Save the plot
file_name <- paste0(out_folder, "gc_content_iorfs")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)



###### Aromaticity ######
data_aro <- add_dummy_rows(data, "aromaticity", n_bins)
data_aro$type <- factor(data_aro$type, levels = c("cds", "trg", "denovo", "iorf"))
data_aro_summary <- get_ncds_conditions(data_aro, n_bins)

# Pvals
if (plot_pvals) {
  pval_factor <- 0.06
  only_ns <- FALSE
  y_annotation <- 0.4
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
              scale_y = seq(0, 0.3, 0.05))
p

# Save the plot
file_name <- paste0(out_folder, "aromaticity")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)




###### Instability index ######
data_inst <- add_dummy_rows(data, "instability", n_bins)
data_inst$type <- factor(data_inst$type, levels = c("cds", "trg", "denovo", "iorf"))
data_inst_summary <- get_ncds_conditions(data_inst, n_bins)

# Pvals
if (plot_pvals) {
  pval_factor <- 0.06
  only_ns <- FALSE
  y_annotation <- 180
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}

# Plot
p <- get_plot(data_inst,
              data_inst_summary,
              "instability",
              "Instability index",
              n_y_pos = -28,
              print_pval = pval_vect,
              scale_y = seq(0, 125, 25))
p

# Save the plot
file_name <- paste0(out_folder, "instability_index")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)



###### Flexibility #######
data_flex <- add_dummy_rows(data, "mean_flexibility", n_bins)
data_flex$type <- factor(data_flex$type, levels = c("cds", "trg", "denovo", "iorf"))
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

file_name <- paste0(out_folder, "mean_flexibility")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)



###### Hydrophobicity ######
data_hyd <- add_dummy_rows(data, "hydropathy", n_bins)
data_hyd$type <- factor(data_hyd$type, levels = c("cds", "trg", "denovo", "iorf"))
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
              n_y_pos = -2.2,
              print_pval = pval_vect,
              scale_y = seq(-2, 2, 1))
p

file_name <- paste0(out_folder, "hydrophobicity")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)



###### HCA ######
data_hca <- add_dummy_rows(data, "hca", n_bins)
data_hca$type <- factor(data_hca$type, levels = c("cds", "trg", "denovo", "iorf"))
data_hca_summary <- get_ncds_conditions(data_hca, n_bins)

# Pvals
if (plot_pvals) {
  pval_factor <- 0.05
  only_ns <- FALSE
  y_annotation <- 14
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
              scale_y = seq(-10, 10, 4))
p

file_name <- paste0(out_folder, "hca")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)



###### Intrinsic disorder (IUpred) ######
data_disord <- add_dummy_rows(data, "disord", n_bins)
data_disord$type <- factor(data_disord$type, levels = c("cds", "trg", "denovo", "iorf"))
data_disord_summary <- get_ncds_conditions(data_disord, n_bins)

# Pvals
if (plot_pvals) {
  pval_factor <- 0.05
  only_ns <- FALSE
  y_annotation <- 1.4
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}

# Plot
p <- get_plot(data_disord,
              data_disord_summary,
              "disord",
              "Intrinsic disorder (IUPred)",
              n_y_pos = -0.1,
              print_pval = pval_vect,
              scale_y = seq(0, 1, 0.2))
p

# Save the plot
file_name <- paste0(out_folder, "intrinsic_disorder")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)



###### Aggregation (tango) ######
data_agg <- add_dummy_rows(data, "aggreg", n_bins)
data_agg$type <- factor(data_agg$type, levels = c("cds", "trg", "denovo", "iorf"))
data_agg_summary <- get_ncds_conditions(data_agg, n_bins)

# Pvals
if (plot_pvals) {
  pval_factor <- 0.05
  only_ns <- FALSE
  y_annotation <- 1.2
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}

# Plot
p <- get_plot(data_agg,
              data_agg_summary,
              "aggreg",
              "Aggregation (Tango)",
              n_y_pos = -0.1,
              print_pval = pval_vect,
              scale_y = seq(0, .8, 0.2))
p

# Save the plot
file_name <- paste0(out_folder, "aggregation")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)


###### AA use ######
### polar ###
data_polar <- add_dummy_rows(data, "polar_use", n_bins)
data_polar$type <- factor(data_polar$type, levels = c("cds", "trg", "denovo", "iorf"))
data_polar_summary <- get_ncds_conditions(data_polar, n_bins)
# Pvals
if (plot_pvals) {
  pval_factor <- 0.05
  only_ns <- FALSE
  y_annotation <- 0.4
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}
# Plot
p <- get_plot(data_polar,
              data_polar_summary,
              "polar_use",
              "Polar AA distribution (S, T, N, Q)",
              n_y_pos = -0.05,
              print_pval = pval_vect,
              scale_y = seq(0, 0.5, 0.1))
p
# Save the plot
file_name <- paste0(out_folder, "aa_polar_use")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)


### hydrophobic ###
data_hydrophobic <- add_dummy_rows(data, "hydrophobic_use", n_bins)
data_hydrophobic$type <- factor(data_hydrophobic$type, levels = c("cds", "trg", "denovo", "iorf"))
data_hydrophobic_summary <- get_ncds_conditions(data_hydrophobic, n_bins)
# Pvals
if (plot_pvals) {
  pval_factor <- 0.05
  only_ns <- FALSE
  y_annotation <- 0.9
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}
# Plot
p <- get_plot(data_hydrophobic,
              data_hydrophobic_summary,
              "hydrophobic_use",
              "Hydrophobic AA distribution (M, Y, V, L, I, F, W)",
              n_y_pos = -0.05,
              print_pval = pval_vect,
              scale_y = seq(0, 0.8, 0.2))
p
# Save the plot
file_name <- paste0(out_folder, "aa_hydrophobic_use")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)


### positive ###
data_pos <- add_dummy_rows(data, "positive_use", n_bins)
data_pos$type <- factor(data_pos$type, levels = c("cds", "trg", "denovo", "iorf"))
data_pos_summary <- get_ncds_conditions(data_pos, n_bins)
# Pvals
if (plot_pvals) {
  pval_factor <- 0.05
  only_ns <- FALSE
  y_annotation <- 0.4
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}
# Plot
p <- get_plot(data_pos,
              data_pos_summary,
              "positive_use",
              "Positive AA distribution (K, R, H)",
              n_y_pos = -0.02,
              print_pval = pval_vect,
              scale_y = seq(0, 0.3, 0.1))
p
# Save the plot
file_name <- paste0(out_folder, "aa_positive_use")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)


### negative ###
data_neg <- add_dummy_rows(data, "negative_use", n_bins)
data_neg$type <- factor(data_neg$type, levels = c("cds", "trg", "denovo", "iorf"))
data_neg_summary <- get_ncds_conditions(data_neg, n_bins)
# Pvals
if (plot_pvals) {
  pval_factor <- 0.05
  only_ns <- FALSE
  y_annotation <- 0.5
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}
# Plot
p <- get_plot(data_neg,
              data_neg_summary,
              "negative_use",
              "Negative AA distribution (D, E)",
              n_y_pos = -0.02,
              print_pval = pval_vect,
              scale_y = seq(0, 0.3, 0.1))
p
# Save the plot
file_name <- paste0(out_folder, "aa_negative_use")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)


### proline glycine ### 
data_pro_gly <- add_dummy_rows(data, "proline.glycine_use", n_bins)
data_pro_gly$type <- factor(data_pro_gly$type, levels = c("cds", "trg", "denovo", "iorf"))
data_pro_gly_summary <- get_ncds_conditions(data_pro_gly, n_bins)
# Pvals
if (plot_pvals) {
  pval_factor <- 0.05
  only_ns <- FALSE
  y_annotation <- 0.5
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}
# Plot
p <- get_plot(data_pro_gly,
              data_pro_gly_summary,
              "proline.glycine_use",
              "Proline and glycine distribution",
              n_y_pos = -0.02,
              print_pval = pval_vect,
              scale_y = seq(0, 0.3, 0.1))
p
# Save the plot
file_name <- paste0(out_folder, "aa_proline_glycine_use")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)


### cysteine ###
data_cys <- add_dummy_rows(data, "cysteine_use", n_bins)
data_cys$type <- factor(data_cys$type, levels = c("cds", "trg", "denovo", "iorf"))
data_cys_summary <- get_ncds_conditions(data_cys, n_bins)
# Pvals
if (plot_pvals) {
  pval_factor <- 0.05
  only_ns <- FALSE
  y_annotation <- 0.2
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}
# Plot
p <- get_plot(data_cys,
              data_cys_summary,
              "cysteine_use",
              "Cysteine distribution",
              n_y_pos = -0.01,
              print_pval = pval_vect,
              scale_y = seq(0, 0.2, 0.05))
p
# Save the plot
file_name <- paste0(out_folder, "aa_cysteine_use")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)

### alanine ###
data_ala <- add_dummy_rows(data, "alanine_use", n_bins)
data_ala$type <- factor(data_ala$type, levels = c("cds", "trg", "denovo", "iorf"))
data_ala_summary <- get_ncds_conditions(data_ala, n_bins)
# Pvals
if (plot_pvals) {
  pval_factor <- 0.05
  only_ns <- FALSE
  y_annotation <- 0.23
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}
# Plot
p <- get_plot(data_ala,
              data_ala_summary,
              "alanine_use",
              "Alanine distribution",
              n_y_pos = -0.02,
              print_pval = pval_vect,
              scale_y = seq(0, 0.15, 0.05))
p
# Save the plot
file_name <- paste0(out_folder, "aa_alanine_use")
if (use_violins) file_name <- paste0(file_name, "_violin")
if (ONE_SIDED) file_name <- paste0(file_name, "_one_sided")
file_name <- paste0(file_name, ".png")
ggsave(file_name)




####### genome GC % (control) ######
plot_pvals <- FALSE

### whole genome ###
data_gc_genome <- add_dummy_rows(data, "gc_species", n_bins)
data_gc_genome$type <- factor(data_gc_genome$type, levels = c("cds", "trg", "denovo", "iorf"))
data_gc_genome$value <- data_gc_genome$value / 100 # Convert to fraction
data_gc_genome_summary <- data %>%
  group_by(bin, type) %>%
  summarize(n = n_distinct(genome), .groups = "drop") %>%
  mutate(group = paste0(bin, "_", type))
# Pvals
if (plot_pvals) {
  pval_factor <- 0.05
  only_ns <- FALSE
  y_annotation <- 0.37
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}
# Plot
p <- get_plot(data_gc_genome,
              data_gc_genome_summary,
              "gc_species",
              "GC distribution (whole genome)",
              n_y_pos = 0.37,
              print_pval = pval_vect,
              scale_y = seq(0.4, 0.6, 0.05))
p
# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/gc_species.png"), plot = p)
}


### intergenic GC ###
data_gc_inter_genome <- add_dummy_rows(data, "inter_gc_species", n_bins)
data_gc_inter_genome$type <- factor(data_gc_inter_genome$type, levels = c("cds", "trg", "denovo", "iorf"))
data_gc_inter_genome$value <- data_gc_inter_genome$value / 100 # Convert to fraction
data_gc_inter_genome_summary <- data %>%
  group_by(bin, type) %>%
  summarize(n = n_distinct(genome), .groups = "drop") %>%
  mutate(group = paste0(bin, "_", type))
# Pvals
if (plot_pvals) {
  pval_factor <- 0.05
  only_ns <- FALSE
  y_annotation <- 0.5
  pval_vect <- c(pval_factor, only_ns, y_annotation)
} else {
  pval_vect <- c(NA, NA, NA)
}
# Plot
p <- get_plot(data_gc_inter_genome,
              data_gc_inter_genome_summary,
              "inter_gc_species",
              "GC distribution (intergenic regions)",
              n_y_pos = 0.36,
              print_pval = pval_vect,
              scale_y = seq(0.35, 0.6, 0.05))
p
# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/gc_species_intergenic.png"), plot = p)
}
