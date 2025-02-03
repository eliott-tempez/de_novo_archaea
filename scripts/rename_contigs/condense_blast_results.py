import os
import pandas as pd


BLAST_RESULTS = "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/results/rename_contigs/blast_result.txt"
OUTPUT_FILE = "/home/eltem/Documents/Cours/M2/Stage/M2_stage_I2BC/results/rename_contigs/blast_result_condensed.txt"


### Change the name of the genomes with 2+ contigs ###
dict_contigs = {}

# Read the blast results
data =  pd.read_csv(BLAST_RESULTS, sep="\t", header=None, names=["name", "gbk_contig", "fasta_contig", "evalue", "match_length"])
# Order by genome, contig, then e-value
data = data.sort_values(by=["name", "gbk_contig", "evalue"])
# Delete all lines with an e-value > 0.05
data = data[data["evalue"] <= 0.05]
unique_names = data["name"].unique()
# Extract data for each genome and contig
for name in unique_names:
    dict_contigs[name] = {}
    unique_contigs = data[data["name"] == name]["gbk_contig"].unique()
    for contig in unique_contigs:
        sub_data = data[(data["name"] == name) & (data["gbk_contig"] == contig)]
        dict_contigs[name][contig] = {}
        # Count each contig occurence
        count = sub_data["fasta_contig"].value_counts()
        # Add to dict
        for fasta_contig, nb in count.items():
            dict_contigs[name][contig][fasta_contig] = {}
            dict_contigs[name][contig][fasta_contig]["count"] = nb
            # Get the mean e-value
            dict_contigs[name][contig][fasta_contig]["evalue"] = sub_data[sub_data["fasta_contig"] == fasta_contig]["evalue"].mean()


with open(OUTPUT_FILE, "a") as f_out:
    for name in dict_contigs:
        f_out.write(f"{name}\n")
        for contig in dict_contigs[name]:
            f_out.write(f"{contig}\n")
            for fasta_contig in dict_contigs[name][contig]:
                f_out.write(f"{fasta_contig}\t{dict_contigs[name][contig][fasta_contig]['count']}\t{dict_contigs[name][contig][fasta_contig]['evalue']}\n")
        f_out.write("\n")