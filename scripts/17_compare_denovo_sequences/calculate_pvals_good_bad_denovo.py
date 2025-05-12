import os
import sys
import re
import random
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


def get_median_diff(df1, df2):
    # Get the median differences for each descriptor
    median_diff = {}
    for descriptor in df1.columns:
        if descriptor not in ["genome", "cds", "type"]:
            median_diff[descriptor] = abs(df1[descriptor].median() - df2[descriptor].median())
    return median_diff


def pool_cdss(df1, df2):
    # Pool the samples and get 2 random subsets
    pooled_df = pd.concat([df1, df2]).reset_index(drop=True)
    n = len(pooled_df)
    sampled_indexes_1 = random.sample(range(n), int(n/2))
    sampled_indexes_2 = [i for i in range(n) if i not in sampled_indexes_1]
    pooled_df_1 = pooled_df.iloc[sampled_indexes_1]
    pooled_df_2 = pooled_df.iloc[sampled_indexes_2]
    return pooled_df_1, pooled_df_2


def compare_medians(median_diff, random_diff):
    # Compare the median differences
    results = {}
    for descriptor in median_diff:
        if median_diff[descriptor] > random_diff[descriptor]:
            results[descriptor] = 1
        else:
            results[descriptor] = 0
    # Create a dataframe with the results
    results_lst = list(results.values())
    results_df = pd.DataFrame([results_lst], columns=list(results.keys()))
    return results_df


def calculate_pvalues(signif_results, descriptors):
    results = []
    # Get the number of 1s for each column
    signif = signif_results[descriptors].sum()
    # Get the number of samples
    n_samples = len(signif_results)
    # Get the pvalues
    pvalues = [1 - (x / n_samples) for x in signif]
    results.append(list(pvalues))
    # Create a dataframe with the results
    results_df = pd.DataFrame(results, columns=descriptors)
    # Save the results
    results_df.to_csv("pvalues_good_bad.csv", sep="\t", index=False)
    return results_df



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
    n = 1
    n_to_sample = min([len(good_indexes), len(bad_indexes)])
    signif_columns = descriptors
    signif_results = pd.DataFrame(columns=signif_columns)


    for i in range(n):
        # Sample the data
        good_sample = random.sample(good_indexes, n_to_sample)
        bad_sample = random.sample(bad_indexes, n_to_sample)
        good_df = descriptors_df.iloc[good_sample]
        bad_df = descriptors_df.iloc[bad_sample]
        print(good_df, bad_df)
        # Get the median difference
        median_diff = get_median_diff(good_df, bad_df)
        print(median_diff)

        # Pool the data and get the median diff
        pool1, pool2 = pool_cdss(good_df, bad_df)
        print(pool1, pool2)
        random_diff = get_median_diff(pool1, pool2)
        print(random_diff)

        # Compare the medians
        result = compare_medians(median_diff, random_diff)
        print(result)
        signif_results = pd.concat([signif_results, result], ignore_index=True)

    # Calculate p-values and save to file
    pvalues = calculate_pvalues(signif_results, descriptors)
    print("\nDone!")
    print(signif_results)
    print(pvalues)