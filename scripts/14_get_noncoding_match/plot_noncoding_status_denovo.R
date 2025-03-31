library(ggVennDiagram)
library(ggplot2)
ggsave <- function(..., bg = "white",
                   width = 1000, height = 1000,
                   units = "px", dpi = 100) {
  ggplot2::ggsave(..., bg = bg,
                  width = width,
                  height = height,
                  units = units,
                  dpi = dpi)
}


input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/denovo_noncoding_status.tsv"


data <- read.table(input_file, header = TRUE, sep = "\t")
data[, c("intergenic", "f.0", "f.1", "f.2")] <- lapply(data[, c("intergenic", "f.0", "f.1", "f.2")], function(x) ifelse(x < 3, 0, x))



##### ORIGIN FOR EACH DE NOVO GENE #####
# Get the number of genes for each origin
origins <- data[, c("intergenic", "f.0", "f.1", "f.2")]
index_only_intergenic <- rownames(origins[apply(origins, 1, function(x) sum(x == 0) == 3) & origins$intergenic > 0, ])
index_only_f0 <- rownames(origins[apply(origins, 1, function(x) sum(x == 0) == 3) & origins$f.0 > 0, ])
index_only_f1 <- rownames(origins[apply(origins, 1, function(x) sum(x == 0) == 3) & origins$f.1 > 0, ])
index_only_f2 <- rownames(origins[apply(origins, 1, function(x) sum(x == 0) == 3) & origins$f.2 > 0, ])
index_several_origins <- rownames(origins[apply(origins, 1, function(x) sum(x == 0) < 3), ])
n_only_intergenic <- length(index_only_intergenic)
n_only_f0 <- length(index_only_f0)
n_only_f1 <- length(index_only_f1)
n_only_f2 <- length(index_only_f2)
n_several_origins <- length(index_several_origins)
# pie chart
# Combine the values and labels into a data frame for sorting
legend_data <- data.frame(
  value = c(n_only_intergenic, n_only_f0, n_only_f1, n_only_f2, n_several_origins),
  label = c("intergenic", "f+0", "f+1", "f+2", "several origins"),
  color = c("#009E73", "#E69F00", "#CC79A7", "#0072B2", "#D55E00")
)
# Update colors for 0 values to gray
legend_data$color[legend_data$value == 0] <- "gray"
# Sort the data frame by values in descending order
legend_data <- legend_data[order(-legend_data$value), ]
# Extract sorted values, labels, and colors
sorted_values <- legend_data$value
sorted_labels <- legend_data$label
sorted_colors <- legend_data$color
# Create the pie chart
pie(sorted_values,
    labels = NA,
    col = sorted_colors,
    main = "Origin of de novo genes in NC match in outgroup")
# Add the legend with sorted values and labels
legend("topright",
       legend = paste(sorted_labels, " (", sorted_values, ")", sep = ""),
       fill = sorted_colors)



##### VENN DIAGRAM FOR ALL ORIGINS #####
intergenic_genes <- data[data$intergenic > 0, "denovo_gene"]
f0_genes <- data[data$f.0 > 0, "denovo_gene"]
f1_genes <- data[data$f.1 > 0, "denovo_gene"]
f2_genes <- data[data$f.2 > 0, "denovo_gene"]
venn <- list(
  intergenic = intergenic_genes,
  f0 = f0_genes,
  f1 = f1_genes,
  f2 = f2_genes
)
# Venn diagram
ggVennDiagram(venn, label_alpha = 0) +
    ggtitle("Venn diagram of the origin of de novo genes in NC match in outgroup") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_gradientn(
    colors = c("white", "#ffbec8", "#ff7070"),
    values = scales::rescale(c(0, 4, 138)))



##### LENGTH OF THE MATCHES #####
len_only_intergenic <- data[index_only_intergenic, "noncoding_match_end"] - data[index_only_intergenic, "noncoding_match_start"]
len_only_f0 <- data[index_only_f0, "noncoding_match_end"] - data[index_only_f0, "noncoding_match_start"]
len_only_f1 <- data[index_only_f1, "noncoding_match_end"] - data[index_only_f1, "noncoding_match_start"]
len_only_f2 <- data[index_only_f2, "noncoding_match_end"] - data[index_only_f2, "noncoding_match_start"]
len_several_origins <- data[index_several_origins, "noncoding_match_end"] - data[index_several_origins, "noncoding_match_start"]
# Boxplot
boxplot(len_only_intergenic, len_only_f1, len_several_origins,
        names = c("intergenic", "f+1", "several origins"),
        col = c("#009E73", "#CC79A7", "#D55E00"),
        xlab = "Origin in the NC region of the outgroup",
        ylab = "Length of the match in the de novo gene (nt)",
        main = "Length distribution of the match between the de novo\ngene and the noncoding region of the outgroup")
# Add text annotations for the number of values
text(x = 1, y = 35,
     labels = paste("n =", length(len_only_intergenic)), pos = 3, cex = 0.8)
text(x = 2, y = 35,
     labels = paste("n =", length(len_only_f1)), pos = 3, cex = 0.8)
text(x = 3, y = 35,
     labels = paste("n =", length(len_several_origins)), pos = 3, cex = 0.8)



##### LENGTH OF THE MATCHES FOR GENES COMING FROM SEVERAL ORIGINS #####
len_intergenic_several_origins <- data[index_several_origins, "intergenic"]
len_intergenic_several_origins <- len_intergenic_several_origins[len_intergenic_several_origins > 0]
len_f0_several_origins <- data[index_several_origins, "f.0"]
len_f0_several_origins <- len_f0_several_origins[len_f0_several_origins > 0]
len_f1_several_origins <- data[index_several_origins, "f.1"]
len_f1_several_origins <- len_f1_several_origins[len_f1_several_origins > 0]
len_f2_several_origins <- data[index_several_origins, "f.2"]
len_f2_several_origins <- len_f2_several_origins[len_f2_several_origins > 0]
# Boxplot
boxplot(len_intergenic_several_origins, len_f0_several_origins, len_f1_several_origins, len_f2_several_origins,
        names = c("intergenic", "f+0", "f+1", "f+2"),
        col = c("#009E73", "#E69F00", "#CC79A7", "#0072B2"),
        xlab = "Frame of origin",
        ylab = "Length of the match in the de novo gene (nt)",
        main = "Length distribution of the match between the de novo\ngene and the noncoding region of the outgroup\nfor genes coming from several origins",
        ylim = c(-10, max(c(len_intergenic_several_origins, len_f0_several_origins, len_f1_several_origins, len_f2_several_origins))))
# Add text annotations for the number of values
text(x = 1, y = -10, labels = paste("n =", length(len_intergenic_several_origins)), pos = 3, cex = 0.8)
text(x = 2, y = -10, labels = paste("n =", length(len_f0_several_origins)), pos = 3, cex = 0.8)
text(x = 3, y = -10, labels = paste("n =", length(len_f1_several_origins)), pos = 3, cex = 0.8)
text(x = 4, y = -10, labels = paste("n =", length(len_f2_several_origins)), pos = 3, cex = 0.8)
