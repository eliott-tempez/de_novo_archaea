###################################################################
############################## IMPORT #############################
###################################################################
# Librairies
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(ggpubr)


save_plots <- TRUE
home <- TRUE


# Files
if (!home){
  input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/17_compare_denovo_sequences/sequence_features_good_candidates_all.csv"
  good_candidates_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/good_candidates.txt"
  all_denovo_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/all_denovo.txt"
  pvals_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/17_compare_denovo_sequences/good_bad/pvalues_good_bad_denovo_one_sided.csv"
  out_folder <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/17_compare_denovo_sequences/good_bad/"
} else {
  input_file <- "/home/eliott/Documents/UNI/M2/Stage/M2_stage_I2BC/results/17_compare_denovo_sequences/sequence_features_good_candidates_all.csv"
  good_candidates_file <- "/home/eliott/Documents/UNI/M2/Stage/M2_stage_I2BC/results/14_get_noncoding_match/good_candidates.txt"
  all_denovo_file <- "/home/eliott/Documents/UNI/M2/Stage/M2_stage_I2BC/results/14_get_noncoding_match/all_denovo.txt"
  pvals_file <- "/home/eliott/Documents/UNI/M2/Stage/M2_stage_I2BC/results/17_compare_denovo_sequences/good_bad/pvalues_good_bad_denovo_one_sided.csv"
  out_folder <- "/home/eliott/Documents/UNI/M2/Stage/M2_stage_I2BC/results/17_compare_denovo_sequences/good_bad/"
}


# Read data
data <- read.table(input_file, header = TRUE, sep = "\t")
pvals <- read.table(pvals_file, header = TRUE, sep = "\t")
# Pivot to longer format
descriptors <- setdiff(colnames(data), c("genome", "cds", "type"))
not_present <- setdiff(descriptors, colnames(pvals))
descriptors <- setdiff(descriptors, not_present)
data <- pivot_longer(data, cols = all_of(descriptors), names_to = "feature", values_to = "value")
pvals <- pivot_longer(pvals, cols = all_of(descriptors), names_to = "feature", values_to = "pval")
pvals$p <- pvals$pval
signif_label <- "****: p <= 1e-5    ***: p <= 1e-4    **: p <= 1e-3    *: p <= 0.05"
# Change the de novo labels
all_denovo <- scan(all_denovo_file, what = "", sep = "\n")
good_candidates <- scan(good_candidates_file, what = "", sep = "\n")
bad_candidates <- setdiff(all_denovo, good_candidates)
data[data$cds %in% good_candidates, "type"] <- "good"
data[data$cds %in% bad_candidates, "type"] <- "bad"
# Keep only relevant data
data <- data[data$type %in% c("good", "bad"), ]



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


get_ncds_conditions <- function(data) {
  # Get the data
  data_len_summary <- data %>%
    group_by(type) %>%
    summarise(n = n(), .groups = "drop")
  return(data_len_summary)
}



# Get the p-values in a matrix to print on the graph
get_pvals <- function(desc, data, fact, non_signif = FALSE) {
  local_pvals <- pvals[pvals$feature == desc, "p"]
  local_pvals$group1 <- "good"
  local_pvals$group2 <- "bad"
  # y position
  local_pvals$y.position <- fact
  
  # significance level
  if (non_signif) {
    local_pvals$p.signif <- case_when(
      local_pvals$p >= 0.05 ~ "ns",
      TRUE ~ ""
    )
  } else {
    ast <- case_when(
      local_pvals$p <= 1e-5 ~ "****",
      local_pvals$p <= 1e-4 ~ "***",
      local_pvals$p <= 1e-3 ~ "**",
      local_pvals$p <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
    local_pvals$p.signif <- paste0(ast, "\np = ", round(local_pvals$p, 5))
  }
  print(local_pvals)
  return(local_pvals)
}


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
get_plot <- function(data, data_summary, feature, title, print_pval = c(NA), scale_y = c(NA), n_y_pos = NA) {

  data_sans_outliers <- remove_outliers(data)

  p <- ggplot(data, aes(x = type, y = value, fill = type)) +
    geom_violin(data = data_sans_outliers, na.rm = TRUE, colour = "#2c2c2c", scale = "width", alpha = 0.8) +
    geom_boxplot(na.rm = TRUE, colour = "#2c2c2c", outliers = FALSE, width = 0.2, alpha = 0.8) +
    labs(title = "",
         x = "",
         y = title) +
    scale_fill_manual(values = c("#4f535a", "#179207")) +
    theme_minimal() +
    scale_x_discrete(labels = c("Integrity", "No integrity")) +
    theme(axis.text.x = element_text(size = 22),
          axis.title.y = element_text(size = 24),
          axis.text.y = element_text(size = 22),
          legend.position = "none")


  if (!all(is.na(scale_y))) {
    p <- p + scale_y_continuous(breaks = scale_y)
  }

  if (!is.na(n_y_pos)) {
    p <- p +
      geom_text(data = data_summary,
                aes(x = type,
                    y = n_y_pos,
                    label = paste0("n = ", n),
                    group = type),
                position = position_dodge(width = 0.75), 
                vjust = -0.5, size = 8)
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
                           hide.ns = FALSE,
                           tip.length = 0.01,
                           label.size = 10)
    }
  }
  return(p)
}







###################################################################
############################### MAIN ##############################
###################################################################


###### Sequence length ######
data_len <- data[data$feature == "length", ]
# Use aa length
data_len$value <- data_len$value / 3
data_summary <- get_ncds_conditions(data_len)
data_len$type <- factor(data_len$type, levels = c("bad", "good"))


# Pvals
pval_pos <- 210
only_ns <- FALSE
y_annotation <- 400
pval_vect <- c(pval_pos, only_ns, y_annotation)

# Plot
p <- get_plot(data_len,
              data_summary,
              "length",
              "Sequence length (residues)",
              n_y_pos = 5,
              print_pval = pval_vect,
              scale_y = seq(0, 300, 50))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/sequence_length.png"), plot = p)
}



###### GC ratio ######
data_gc <- data[data$feature == "gc_rate", ]
data_gc$type <- factor(data_gc$type, levels = c("bad", "good"))

# Pvals
pval_pos <- 1.25
only_ns <- TRUE
y_annotation <- 0.55
pval_vect <- c(pval_pos, only_ns, y_annotation)

# Plot
p <- get_plot(data_gc,
              data_summary,
              "gc_rate",
              "GC ratio: sequence GC % / genome GC %",
              n_y_pos = 0.5,
              print_pval = pval_vect,
              scale_y = seq(0, 1.5, 0.2))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/gc_rate.png"), plot = p)
}



###### GC ratio (intergenic) ######
data_gc_inter <- data[data$feature == "iorfs_gc_rate", ]
data_gc_inter$type <- factor(data_gc_inter$type, levels = c("bad", "good"))

# Pvals
pval_pos <- 1.3
only_ns <- TRUE
y_annotation <- 0.55
pval_vect <- c(pval_pos, only_ns, y_annotation)

# Plot
p <- get_plot(data_gc_inter,
              data_summary,
              "inter_gc_rate",
              "GC ratio: sequence GC % / iORF GC %",
              n_y_pos = 0.58,
              print_pval = pval_vect,
              scale_y = seq(0, 1.5, 0.1))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/gc_rate_intergenic.png"), plot = p)
}



##### Aromaticity ######
data_arom <- data[data$feature == "aromaticity", ]
data_arom$type <- factor(data_arom$type, levels = c("bad", "good"))

# Pvals
pval_pos <- 0.26
only_ns <- TRUE
y_annotation <- 0.55
pval_vect <- c(pval_pos, only_ns, y_annotation)

# Plot
p <- get_plot(data_arom,
              data_summary,
              "aromaticity",
              "Aromaticity",
              n_y_pos = -0.02,
              print_pval = pval_vect,
              scale_y = seq(0, 0.2, 0.1))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/aromaticity.png"), plot = p)
}




###### Instability index ######
data_instab <- data[data$feature == "instability", ]
data_instab$type <- factor(data_instab$type, levels = c("bad", "good"))

# Pvals
pval_pos <- 95
only_ns <- TRUE
y_annotation <- 0.55
pval_vect <- c(pval_pos, only_ns, y_annotation)

# Plot
p <- get_plot(data_instab,
              data_summary,
              "instability",
              "Instability index",
              n_y_pos = -5,
              print_pval = pval_vect,
              scale_y = seq(0, 100, 20))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/instability.png"), plot = p)
}




###### Flexibility #######
data_flex <- data[data$feature == "mean_flexibility", ]
data_flex$type <- factor(data_flex$type, levels = c("bad", "good"))

# Pvals
pval_pos <- 1.04
only_ns <- TRUE
y_annotation <- 1.05
pval_vect <- c(pval_pos, only_ns, y_annotation)

# Plot
p <- get_plot(data_flex,
              data_summary,
              "mean_flexibility",
              "Flexibility index",
              n_y_pos = 0.96,
              print_pval = pval_vect,
              scale_y = seq(0.96, 1.04, 0.02))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/flexibility.png"), plot = p)
}




###### Hydrophobicity ######
data_hydrop <- data[data$feature == "hydropathy", ]
data_hydrop$type <- factor(data_hydrop$type, levels = c("bad", "good"))

# Pvals
pval_pos <- 1.7
only_ns <- TRUE
y_annotation <- 3
pval_vect <- c(pval_pos, only_ns, y_annotation)

# Plot
p <- get_plot(data_hydrop,
              data_summary,
              "hydropathy",
              "Hydrophobicity index",
              n_y_pos = -1.8,
              print_pval = pval_vect,
              scale_y = seq(-2, 2, 1))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/hydropathy.png"), plot = p)
}




###### HCA ######
data_hca <- data[data$feature == "hca", ]
data_hca$type <- factor(data_hca$type, levels = c("bad", "good"))

# Pvals
pval_pos <- 11
only_ns <- TRUE
y_annotation <- 8.5
pval_vect <- c(pval_pos, only_ns, y_annotation)

# Plot
p <- get_plot(data_hca,
              data_summary,
              "hca",
              "HCA index",
              n_y_pos = -11,
              print_pval = pval_vect,
              scale_y = seq(-10, 10, 4))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/hca.png"), plot = p)
}




###### Intrinsic disorder (IUpred) ######
data_iupred <- data[data$feature == "disord", ]
data_iupred$type <- factor(data_iupred$type, levels = c("bad", "good"))

# Pvals
pval_pos <- 0.9
only_ns <- TRUE
y_annotation <- 0.1
pval_vect <- c(pval_pos, only_ns, y_annotation)

# Plot
p <- get_plot(data_iupred,
              data_summary,
              "disord",
              "Intrinsic disorder (IUPred)",
              n_y_pos = -0.05,
              print_pval = pval_vect,
              scale_y = seq(0, 1, 0.2))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/iupred.png"), plot = p)
}




###### Aggregation (tango) ######
data_tango <- data[data$feature == "aggreg", ]
data_tango$type <- factor(data_tango$type, levels = c("bad", "good"))

# Pvals
pval_pos <- 0.7
only_ns <- TRUE
y_annotation <- 0.55
pval_vect <- c(pval_pos, only_ns, y_annotation)

# Plot
p <- get_plot(data_tango,
              data_summary,
              "aggreg",
              "Aggregation (Tango)",
              n_y_pos = -0.05,
              print_pval = pval_vect,
              scale_y = seq(0, 1, 0.2))
p

# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/tango.png"), plot = p)
}







###### AA use ######
### polar ###
data_polar <- data[data$feature == "polar_use", ]
data_polar$type <- factor(data_polar$type, levels = c("bad", "good"))
# Pvals
pval_pos <- 0.37
only_ns <- TRUE
y_annotation <- 0.55
pval_vect <- c(pval_pos, only_ns, y_annotation)
# Plot
p <- get_plot(data_polar,
              data_summary,
              "polar_use",
              "Polar AA distribution (S, T, N, Q)",
              n_y_pos = 0.07,
              print_pval = pval_vect,
              scale_y = seq(0, 1, 0.1))
p
# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/aa_polar_use.png"), plot = p)
}

### hydrophobic ###
data_hydrophobic <- data[data$feature == "hydrophobic_use", ]
data_hydrophobic$type <- factor(data_hydrophobic$type, levels = c("bad", "good"))
# Pvals
pval_pos <- 0.63
only_ns <- TRUE
y_annotation <- 0.55
pval_vect <- c(pval_pos, only_ns, y_annotation)
# Plot
p <- get_plot(data_hydrophobic,
              data_summary,
              "hydrophobic_use",
              "Hydrophobic AA distribution (M, Y, V, L, I, F, W)",
              n_y_pos = 0.07,
              print_pval = pval_vect,
              scale_y = seq(0, 1, 0.1))
p
# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/aa_hydrophobic_use.png"), plot = p)
}

### positive ###
data_positive <- data[data$feature == "positive_use", ]
data_positive$type <- factor(data_positive$type, levels = c("bad", "good"))
# Pvals
pval_pos <- 0.36
only_ns <- TRUE
y_annotation <- 0.55
pval_vect <- c(pval_pos, only_ns, y_annotation)
# Plot
p <- get_plot(data_positive,
              data_summary,
              "positive_use",
              "Positive AA distribution (K, R, H)",
              n_y_pos = -0.03,
              print_pval = pval_vect,
              scale_y = seq(0, 1, 0.1))
p
# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/aa_positive_use.png"), plot = p)
}

### negative ###
data_negative <- data[data$feature == "negative_use", ]
data_negative$type <- factor(data_negative$type, levels = c("bad", "good"))
# Pvals
pval_pos <- 0.28
only_ns <- TRUE
y_annotation <- 0.55
pval_vect <- c(pval_pos, only_ns, y_annotation)
# Plot
p <- get_plot(data_negative,
              data_summary,
              "negative_use",
              "Negative AA distribution (D, E)",
              n_y_pos = -0.03,
              print_pval = pval_vect,
              scale_y = seq(0, 1, 0.1))
p
# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/aa_negative_use.png"), plot = p)
}

### proline glycine ### 
data_proline_glycine <- data[data$feature == "proline.glycine_use", ]
data_proline_glycine$type <- factor(data_proline_glycine$type, levels = c("bad", "good"))
# Pvals
pval_pos <- 0.26
only_ns <- TRUE
y_annotation <- 0.55
pval_vect <- c(pval_pos, only_ns, y_annotation)
# Plot
p <- get_plot(data_proline_glycine,
              data_summary,
              "proline.glycine_use",
              "Proline and glycine AA distribution",
              n_y_pos = -0.03,
              print_pval = pval_vect,
              scale_y = seq(0, 1, 0.1))
p
# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/aa_proline-glycine_use.png"), plot = p)
}

### cysteine ###
data_cysteine <- data[data$feature == "cysteine_use", ]
data_cysteine$type <- factor(data_cysteine$type, levels = c("bad", "good"))
# Pvals
pval_pos <- 0.045
only_ns <- TRUE
y_annotation <- 0.55
pval_vect <- c(pval_pos, only_ns, y_annotation)
# Plot
p <- get_plot(data_cysteine,
              data_summary,
              "cysteine_use",
              "Cysteine AA distribution",
              n_y_pos = -0.005,
              print_pval = pval_vect,
              scale_y = seq(0, 1, 0.01))
p
# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/aa_cysteine_use.png"), plot = p)
}

### alanine ###
data_alanine <- data[data$feature == "alanine_use", ]
data_alanine$type <- factor(data_alanine$type, levels = c("bad", "good"))
# Pvals
pval_pos <- 0.17
only_ns <- TRUE
y_annotation <- 0.55
pval_vect <- c(pval_pos, only_ns, y_annotation)
# Plot
p <- get_plot(data_alanine,
              data_summary,
              "alanine_use",
              "Alanine AA distribution",
              n_y_pos = -0.01,
              print_pval = pval_vect,
              scale_y = seq(0, 1, 0.05))
p
# Save the plot
if (save_plots) {
  ggsave(paste0(out_folder, "/aa_alanine_use.png"), plot = p)
}