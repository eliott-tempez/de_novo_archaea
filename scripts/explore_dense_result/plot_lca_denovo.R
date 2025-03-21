library(ape)
library(phytools)

input_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/explore_dense_results/denovo_noncoding_status.tsv"
tree_file <- "/home/eliott.tempez/Documents/archaea_data/whole_tree.nwk"

data <- read.table(input_file, header = TRUE, sep = "\t")
tree <- read.tree(tree_file)


# Function to get LCA node and age
get_lca_age <- function(tree, genome, outgroup) {
  lca_node <- findMRCA(tree, tips = c(genome, outgroup), type = "node")
  node_heights <- node.depth.edgelength(tree)  # Compute node heights
  lca_age <- max(node_heights) - node_heights[lca_node]  # Get relative age
  return(data.frame(lca_node = lca_node, lca_age = lca_age))
}


lca_results <- do.call(rbind, apply(data, 1, function(row) get_lca_age(tree, row["genome"], row["outgroup"])))

# Merge results with original dataframe
data <- cbind(data, lca_results)

# Print results
print(data)
