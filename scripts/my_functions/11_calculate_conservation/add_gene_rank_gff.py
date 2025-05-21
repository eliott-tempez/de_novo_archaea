import re
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors



GENOMES_LIST = "/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/genera_archaea/genomes_list.txt"
GENERA_DIR = "/home/eliott.tempez/Documents/archaea_data/genera/out/"
GFF_DIR = "/home/eliott.tempez/Documents/archaea_data/complete_122/gff3_no_fasta/"
OUTPUT_DIR = "/home/eliott.tempez/Documents/archaea_data/complete_122/gff3_coloured/"



if __name__ == "__main__":
    # Read the list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]
    # Create the output directory if doesn't exist
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # Create a color palette from black to red with 10 colors
    custom_cmap = mcolors.LinearSegmentedColormap.from_list("black_red", ["black", "red"])
    colors = [custom_cmap(i / 9) for i in range(10)]
    hex_colors = [mcolors.rgb2hex(color) for color in colors]

    # For each genome
    for genome in genomes:
        # Check if we have the gene age data
        gene_age_pattern = os.path.join(GENERA_DIR, genome, "*gene_ages.tsv")
        gene_age_files = glob.glob(gene_age_pattern)
        if not gene_age_files:
            continue
        gene_age_file = gene_age_files[0]
        gene_age_df = pd.read_csv(gene_age_file, sep="\t", header=0)
        # Check we have the gff file
        gff_file = os.path.join(GFF_DIR, genome + ".gff3")
        if not os.path.exists(gff_file):
            raise FileNotFoundError(f"Missing gff file for {genome}")
        
        # Open the input and output file
        output_file = os.path.join(OUTPUT_DIR, genome + ".gff3")
        with open(output_file, "w") as f_out:
            with open(gff_file, "r") as f_in:
                # Re-write the gff file with the color corresponding to the rank
                for line in f_in:
                    if line.startswith("#"):
                        f_out.write(line)
                        continue
                    line_split = line.split("\t")
                    if line_split[2] == "mRNA":
                        gene_name = re.search("ID=(\S+);", line_split[8]).group(1)
                        # Get the gene age
                        gene_rank = gene_age_df.loc[gene_age_df["#gene"] == gene_name, "rank"]
                        gene_rank = int(gene_rank.values[0]) if not gene_rank.empty and not pd.isna(gene_rank.values[0]) else None
                        gene_color = hex_colors[gene_rank - 1] if not gene_rank is None else "#6a6b6a"
                        # Add to the gff3 line
                        new_line = f'{line.strip()};color="{gene_color}"\n'
                        f_out.write(new_line)

            
        

