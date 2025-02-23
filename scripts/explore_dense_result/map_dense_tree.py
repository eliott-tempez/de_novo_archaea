import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo

TRG_FILE = "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/data/dense/archive_tree_v1/GCA_001433455@Thermococcus_barophilus_CH5/TRG_match_matrix.tsv"
GENE = "TBCH5v1_RS02920"
TREE = "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/data/whole_tree.nwk"


def get_trg_line():
    """Get the line corresponding to the gene"""
    trg_df = pd.read_csv(TRG_FILE, sep="\t", index_col=0)
    # get the indexes that have the gene name in their names
    gene_indexes = [index for index in trg_df.index if GENE in index]
    if len(gene_indexes) == 0:
        print("No gene found")
        exit()
    elif len(gene_indexes) > 1:
        print("Multiple genes found")
        exit()
    return trg_df.loc[gene_indexes[0]]


def print_tree(tree_file):
    """Print the tree"""
    tree = Phylo.read(tree_file, "newick")
    Phylo.draw(tree)




if __name__ == "__main__":
    # get info on TRG
    trg_line = get_trg_line()
    # print tree
    print_tree(TREE)
    