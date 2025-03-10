import os
import matplotlib.pyplot as plt
plt.style.use('petroff10')

OLD_GFF = "/home/eliott.tempez/Documents/archaea_data/complete_122/gff3_no_fasta/"
NEW_GFF = "/home/eliott.tempez/Documents/archaea_data/complete_122/reannotated_gff/"
OUT_DIR = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/re_annotate_genomes/"


def get_features(gff_file):
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line_split = line.strip().split("\t")
            contig = line_split[0]
            feature_type = line_split[2]
            start = int(line_split[3])
            end = int(line_split[4])
            strand = line_split[6]
            attributes = line_split[8]
            yield contig, feature_type, start, end, strand, attributes


def get_gene_number(gff_file):
    gene_number = 0
    for _, feature, _, _, _, _ in get_features(gff_file):
        if feature == "gene":
            gene_number += 1
    return gene_number


def get_cds_number(gff_file):
    cds_number = 0
    for _, feature, _, _, _, _ in get_features(gff_file):
        if feature == "CDS":
            cds_number += 1
    return cds_number





if __name__ == "__main__":
    
    data_dict = {}

    for genome in os.listdir(OLD_GFF):
        if not genome.endswith(".gff3"):
            continue

        old_gff = os.path.join(OLD_GFF, genome)
        new_gff = os.path.join(NEW_GFF, genome)

        data_dict[genome] = {}
        data_dict[genome]["is_GCA"] = True if genome.startswith("GCA") else False
        data_dict[genome]["old_gene_nb"] = get_gene_number(old_gff)
        data_dict[genome]["new_gene_nb"] = get_gene_number(new_gff)
        data_dict[genome]["old_cds_nb"] = get_cds_number(old_gff)
        data_dict[genome]["new_cds_nb"] = get_cds_number(new_gff)
        

    # Plot results
    old_genes_nb_GCA = [data_dict[genome]["old_gene_nb"] for genome in data_dict if data_dict[genome]["is_GCA"]]
    new_genes_nb_GCA = [data_dict[genome]["new_gene_nb"] for genome in data_dict if data_dict[genome]["is_GCA"]]
    old_genes_nb_not_GCA = [data_dict[genome]["old_gene_nb"] for genome in data_dict if not data_dict[genome]["is_GCA"]]
    new_genes_nb_not_GCA = [data_dict[genome]["new_gene_nb"] for genome in data_dict if not data_dict[genome]["is_GCA"]]
    old_cds_nb_GCA = [data_dict[genome]["old_cds_nb"] for genome in data_dict if data_dict[genome]["is_GCA"]]
    new_cds_nb_GCA = [data_dict[genome]["new_cds_nb"] for genome in data_dict if data_dict[genome]["is_GCA"]]
    old_cds_nb_not_GCA = [data_dict[genome]["old_cds_nb"] for genome in data_dict if not data_dict[genome]["is_GCA"]]
    new_cds_nb_not_GCA = [data_dict[genome]["new_cds_nb"] for genome in data_dict if not data_dict[genome]["is_GCA"]]

    # Create subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    fig.suptitle("Comparison of gene and CDS number between\nold and new annotation (N = 116)")

    # Create boxplots and store the returned dictionaries
    bplot1 = ax1.boxplot([old_genes_nb_GCA, new_genes_nb_GCA], tick_labels=["old", "new"], patch_artist=True)
    bplot2 = ax2.boxplot([old_genes_nb_not_GCA, new_genes_nb_not_GCA], tick_labels=["old", "new"], patch_artist=True)
    bplot3 = ax3.boxplot([old_cds_nb_GCA, new_cds_nb_GCA], tick_labels=["old", "new"], patch_artist=True)
    bplot4 = ax4.boxplot([old_cds_nb_not_GCA, new_cds_nb_not_GCA], tick_labels=["old", "new"], patch_artist=True)

    # Set axis labels
    ax3.set_xlabel("GCA genomes")
    ax4.set_xlabel("Non GCA genomes")  
    ax1.set_ylabel("Number of genes")
    ax3.set_ylabel("Number of CDS")

    # Set the same y-scale for all subplots
    ymin = min([min(old_genes_nb_GCA), min(new_genes_nb_GCA), min(old_genes_nb_not_GCA), min(new_genes_nb_not_GCA), 
                min(old_cds_nb_GCA), min(new_cds_nb_GCA), min(old_cds_nb_not_GCA), min(new_cds_nb_not_GCA)]) - 100
    ymax = max([max(old_genes_nb_GCA), max(new_genes_nb_GCA), max(old_genes_nb_not_GCA), max(new_genes_nb_not_GCA), 
                max(old_cds_nb_GCA), max(new_cds_nb_GCA), max(old_cds_nb_not_GCA), max(new_cds_nb_not_GCA)]) + 100

    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_ylim(ymin, ymax)

    # Hide y labels for the right subplots
    ax2.set_yticklabels([])
    ax4.set_yticklabels([])

    # Set boxplot colors
    for bplot in [bplot1, bplot2, bplot3, bplot4]:
        for patch, color in zip(bplot['boxes'], ['#d9fac4', '#fac4d3']):
            patch.set_facecolor(color)
        for median in bplot['medians']:
            median.set_color('black')


    # Save figure
    fig.savefig(os.path.join(OUT_DIR, "gene_cds_number_comparison.png"))


