"""Get the LCA of the genes with high proportion of Archaea matches"""

import pandas as pd
import os
import subprocess
import tempfile

GENE_NAMES = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/explore_genera_results/archaea_prop_99_genes.txt"
BLAST_RESULT = "/home/eliott.tempez/Documents/archaea_data/re_blast_1000/reblast_output_1000.tsv"
TAXIDS = "/home/eliott.tempez/Documents/archaea_data/re_blast_1000/ncbi_lineages_2025-02-17.csv"

# Load taxonomy
taxid_to_rank = {}
with open(TAXIDS, "r") as f:
    for line in f:
        line_split = line.strip().split(",")
        taxid_to_rank[line_split[0]] = line_split[1]

# Get gene names
gene_names = []
with open(GENE_NAMES, "r") as f:
    for line in f:
        gene_names.append(line.strip())

# For each gene
for gene in gene_names:
    # Create a temporary file with the results for the gene
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as output_file:
        output_file_path = output_file.name
    subprocess.run(f'grep "^{gene}" {BLAST_RESULT} > {output_file_path}', shell=True, check=True)
    # Read the results
    blast_results = pd.read_csv(output_file_path, sep="\t", header=None)
    blast_results.columns = ["query", "subject", "taxids", "evalue", "qcov"]
    taxids = [str(t).split(";")[0] for t in blast_results["taxids"]]
    # Get the superkingdom for each taxid
    not_archaea = 0
    archaea_taxids = []
    for taxid in taxids:
        rank = taxid_to_rank.get(taxid, "Unknown")
        if rank != "Archaea":
            not_archaea += 1
        else:
            archaea_taxids.append(taxid)
    
    # Get the LCA
    unique_archaea_taxids = set(archaea_taxids)
    cmd = f"echo {" ".join(unique_archaea_taxids)} | taxonkit lca"
    result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE)
    lca = result.stdout.decode().strip().split("\t")[-1]

    # Print the results
    print(f"Gene: {gene}")
    print(f"Number of Archaea matches: {len(archaea_taxids)}")
    print(f"Number of non-Archaea matches: {not_archaea}")
    print(f"Archaea LCA: {lca}")
    print("\n")

    # Remove the temporary file
    os.remove(output_file_path)



