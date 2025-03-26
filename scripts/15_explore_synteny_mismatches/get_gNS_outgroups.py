import re
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
pd.set_option("future.no_silent_downcasting", True)


DENSE_DIR = "/home/eliott.tempez/Documents/archaea_data/dense/"
GENOMES_LIST = "/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/genera_archaea/genomes_list.txt"


def get_TRG_profile(genome, gene_name):
    """Get the TRG profile of a gene in a genome"""
    trg_file = os.path.join(DENSE_DIR, genome, "TRG_match_matrix.tsv")
    trg_df = pd.read_csv(trg_file, sep="\t", header=0)
    gene_row = trg_df[trg_df["CDS"] == f"{gene_name}_elongated"]
    # Replace rows with numbers
    gene_row = gene_row.replace(to_replace="noM", value=0)
    gene_row = gene_row.replace(to_replace=r".*gS.*", value=1, regex=True)
    gene_row = gene_row.replace(to_replace=r".*gNS.*", value=2, regex=True)
    gene_row = gene_row.replace(to_replace="CDS", value=3)
    profile_lst = list(gene_row.iloc[0])
    return(profile_lst[1:])

    



if __name__ == "__main__":
    # Read list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]

    # Read all gNS outgroups
    profiles = []
    genes_dict = {}
    gene_names = []
    for genome in genomes:
        gns_list = []
        trg_file = os.path.join(DENSE_DIR, genome, "TRG_match_matrix.tsv")
        trg_df = pd.read_csv(trg_file, sep="\t", header=0)
        # Iterate over the row in reverse until we hit the outgroup
        for row in trg_df.iterrows():
            i = len(row[1]) - 1
            cell = "noM"
            while cell == "noM":
                cell = row[1].iloc[i]
                i -= 1
            if "gNS" in str(cell):
                gene_name = "_".join(row[1].iloc[0].split("_")[:-1])
                gns_list.append(gene_name)
        
        # Get the de novo genes
        denovo_file = os.path.join(DENSE_DIR, genome, "denovogenes.tsv")
        denovo_df = pd.read_csv(denovo_file, sep="\t", header=0)
        if denovo_df.empty:
            continue
        denovo_names = denovo_df["CDS"].tolist()

        # Get rid of the gNS genes that are also de novo
        ourgroup_gns = [g for g in gns_list if g not in denovo_names]
        # Plot the TRG profile of the outgroup gNS genes
        for gene in ourgroup_gns:
            gene_names.append(gene)
            genes_dict[gene] = {"genome": genome, "colnames": trg_df.columns[1:]}
            profiles.append(get_TRG_profile(genome, gene))


    profiles_df = pd.DataFrame(profiles, index=gene_names)
    print(profiles_df.shape)
    # Heatmap in sets of 60 TRGs
    nrows = profiles_df.shape[0]
    for i in range(0, nrows, 60):
        profiles_df_subset = profiles_df.iloc[i:i+60]

        sns.set_theme(font_scale=0.6)
        # Set the colors
        colors = ["#A9A9A9", "#009e28", "#00799e", "#e67b00"]
        cmap = LinearSegmentedColormap.from_list('Custom', colors, len(colors))
        # Figure
        ax = sns.heatmap(profiles_df_subset, cmap=cmap, linewidths=.5, linecolor='lightgray')
        # Manually specify colorbar labelling after it's been generated
        colorbar = ax.collections[0].colorbar
        colorbar.set_ticks([0.375, 1.125, 1.875, 2.625])
        colorbar.set_ticklabels(['noM', 'gS', 'gNS', 'CDS'])
        # Labels
        ax.set_ylabel('')
        ax.set_xlabel('')
        ax.set_xticks([])
        _, labels = plt.yticks()
        plt.setp(labels, rotation=0)
        plt.show()

