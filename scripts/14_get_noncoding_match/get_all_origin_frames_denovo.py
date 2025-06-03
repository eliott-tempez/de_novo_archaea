"""
This script extracts all origins for all the denovo 
good candidates, and all NC matches in synteny.
"""


##### Imports #####
import os
import re
import subprocess
from io import StringIO
import pandas as pd

# Files
GOOD_CANDIDATES_FILE = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/good_candidates.txt"
from my_functions.paths import GENOMES_LIST, DENSE_DIR
OUT_DIR = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/"



##### Functions #####
def get_nc_synteny_matches(denovo, genome):
    nc_genomes = []

    # Get the right row in the trg match matrix
    trg_file = os.path.join(DENSE_DIR, genome, "TRG_match_matrix.tsv")
    trg_df = pd.read_csv(trg_file, sep="\t", header=0)
    try:
        match = trg_df[trg_df["CDS"] == f"{denovo}_elongated"]
        if match.empty:
            raise ValueError(f"No corresponding row found for {denovo}_elongated in {genome}")
    except ValueError as e:
        print(e)
    # Iterate over the row in reverse and get all NC matches
    i = len(match.columns)
    while i > 0:
        i -= 1
        cell = match.iloc[0, i]
        if "gS" in cell:
            outgroup_nb = (cell[2:])
            if outgroup_nb.strip() == "":
                continue
            outgroup_genome = trg_df.columns[i]
            nc_genomes.append((outgroup_genome, int(outgroup_nb)))

    return nc_genomes


def run_get_origin_frames(denovo, genome, nc_genome):
    cmd = [
        "python",
        "14_get_noncoding_match/get_noncoding_match_status.py",
        "-d", denovo, "-f", genome, "-n", nc_genome
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    df = pd.read_csv(StringIO(result.stdout), sep="\t", header=0)
    return(df)




##### Main #####
if __name__ == "__main__":
    merged_df = pd.DataFrame()

    # Read list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]

    # Read good candidates
    good_candidates = []
    with open(GOOD_CANDIDATES_FILE, "r") as f:
        for line in f:
            good_candidates.append(line.strip())

    i = 0
    # Get all de novo good candidates information
    for genome in genomes:
        denovo_file = os.path.join(DENSE_DIR, genome, "denovogenes.tsv")
        # Skip if no denovo genes found
        denovo_df = pd.read_csv(denovo_file, sep="\t", header=0)
        if denovo_df.empty:
            continue
        for denovo_gene in denovo_df["CDS"]:
            if not denovo_gene in good_candidates:
                continue
            i += 1

            # Get all NC synteny matches
            nc_genomes = get_nc_synteny_matches(denovo_gene, genome)
            
            # Extract all origin frames
            for nc_genome, outgroup_nb in nc_genomes:
                result_df = run_get_origin_frames(denovo_gene, genome, nc_genome)
                result_df["outgroup_nb"] = outgroup_nb
                # Append result_df to merged_df
                merged_df = pd.concat([merged_df, result_df], ignore_index=True)
            print(f"Processed {i}/{len(good_candidates)} denovo genes...")

    # Save to file
    merged_df.to_csv(os.path.join(OUT_DIR, "denovo_noncoding_status.tsv"), sep="\t", index=False)





        
