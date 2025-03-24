library(ggplot2)

# Files
tree_file <- "/home/eliott.tempez/Documents/archaea_data/whole_tree.nwk"
genomes_list <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/genera_archaea/genomes_list.txt"


# Get all genome names
genomes <- scan(genomes_list, what = "", sep = " ")[1:116]
genomes <- c("Thermococcus_sp_IRI_33c")


for (gen in genomes) {
  trg_file <- Sys.glob(paste0(genera_dir, gen, "/*gene_ages.tsv"))[1]
  data <- read.table(trg_file, header = TRUE, sep = "\t", comment.char = "")
  print(head(data))
  # Keep only the TRGs
  trgs <- data[data$rank >= 7, ]
  print(head(trgs))
  # Get the number of hits for each rank
  count_phyloranks <- data.frame(table(trgs$rank))
  count_phyloranks <- count_phyloranks[order(count_phyloranks$Var1, decreasing = TRUE), ]
  count_phyloranks$Var1 <- factor(count_phyloranks$Var1, levels = count_phyloranks$Var1[order(-count_phyloranks$Var1)])
  #dotplot
  ggplot(count_phyloranks, aes(x = Var1, y = Freq, size = Freq)) +
    geom_point(color = "black", fill = "#da928b", shape = 21) +
    scale_size_continuous(range = c(3, 10)) +
    labs(x = "Phylorank", y = "Number of hits", title = paste0("Number of TRGs for each phylogenetic age (N = ", nrow(trgs), ")")) +
    theme_minimal() +
    theme(legend.position = "none")

}
