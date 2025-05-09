import os
import sys
import re
import pandas as pd
from multiprocessing import Pool, cpu_count


sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from my_functions.paths import DENSE_DIR, GENOMES_LIST



def extract_denovo_names(focal_species, use_good_candidates=False):
    """Extract the names of the de novo genes"""
    # Make sure the file exists and read it
    if not use_good_candidates:
        denovo_file = os.path.join(DENSE_DIR, focal_species, "denovogenes.tsv")
        if os.path.exists(denovo_file):
            denovo_df = pd.read_csv(denovo_file, sep="\t", header=0)
        else:
            raise FileNotFoundError(f"No file {denovo_file}")
        # Extract the names
        denovo_names = denovo_df["CDS"].tolist()
    else:
        denovo_names = []
        denovo_file = "good_candidates.txt"
        with open(denovo_file, "r") as f:
            for line in f:
                denovo_names.append(line.strip())
    return denovo_names



if __name__ == "__main__":
    # Get the list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]
    denovo_names = []

    # Get the name of all denovo genes
    for genome in genomes:
        # Get the denovo names
        denovo_names += extract_denovo_names(genome, use_good_candidates=False)
    # Get the names of good candidates
    good_candidates = extract_denovo_names(genome, use_good_candidates=True)
    # Get the names of bad candidates
    bad_candidates = list(set(denovo_names) - set(good_candidates))

    # Read descriptors
    descriptors_file = "sequence_features_good_candidates_all.csv"
    descriptors_df = pd.read_csv(descriptors_file, sep="\t", header=0)
    descriptors = list(descriptors_df)
    descriptors.remove("genome")
    descriptors.remove("type")
    descriptors.remove("cds")
    
    
    # Get the indexes for each type of denovo
    good_indexes = []
    bad_indexes = []
    for cds in denovo_names:
        df_line = descriptors_df[descriptors_df["cds"] == cds]
        if df_line.empty:
            raise ValueError(f"No line for {cds} in {descriptors_file}")
        if cds in good_candidates:
            good_indexes.append(df_line.index[0])
        elif cds in bad_candidates:
            bad_indexes.append(df_line.index[0])
        else:
            raise ValueError(f"{cds} not in good or bad candidates")
    
    
    # Number of iterations
    n = 100000
    n_to_sample = min([len(good_indexes), len(bad_indexes)])
    signif_columns = ["good", "bad"] + descriptors
    signif_results = pd.DataFrame(columns=signif_columns)
    
    # Prepare args for each process
    iteration_args = [
        (i, n_to_sample, good_indexes, bad_indexes, descriptors_df, descriptors)
        for i in range(n)
    ]
    
    print("Running in parallel...")
    
    # Initialize progress tracking
    with Pool(cpu_count()) as pool:
        results = []
        