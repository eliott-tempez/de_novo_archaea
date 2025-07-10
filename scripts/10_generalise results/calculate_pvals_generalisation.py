import os
import sys
import re
import random
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from multiprocessing import Pool, cpu_count


sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from my_functions.paths import DENSE_DIR, GENOMES_LIST
INCLUDE = "all" # Descriptors to run the test on
TWO_SIDED = False # If True, the test is two-sided, otherwise one-sided

def get_all_medians(descriptors_df, good_indexes, bad_indexes,descriptors):
    medians = {}
    # Calculate all the medians
    for descriptor in descriptors:
        good_df = descriptors_df.iloc[good_indexes]
        bad_df = descriptors_df.iloc[bad_indexes]
        median_good = good_df[descriptor].median()
        median_bad = bad_df[descriptor].median()
        medians[descriptor] = {
            "good": median_good,
            "bad": median_bad
        }
    return medians


def get_which_way(medians):
    # Determine which way to compare based on the medians if one_sided
    ways = {}
    for descriptor in medians:
        median_good = medians[descriptor]["good"]
        median_bad = medians[descriptor]["bad"]
        if median_good >= median_bad:
            ways[descriptor] = "good>bad"
        else:
            ways[descriptor] = "good<bad"
    return ways


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


def get_median_diff(good_df, bad_df, ways):
    # Get the median differences for each descriptor
    median_diff = {}
    if TWO_SIDED:
        for descriptor in good_df.columns:
            if descriptor not in ["genome", "cds", "type"]:
                median_diff[descriptor] = abs(good_df[descriptor].median() - bad_df[descriptor].median())
    else:
        for descriptor in good_df.columns:
            if descriptor not in ["genome", "cds", "type"]:
                way = ways[descriptor]
                if way == "good>bad":
                    median_diff[descriptor] = good_df[descriptor].median() - bad_df[descriptor].median()
                elif way == "good<bad":
                    median_diff[descriptor] = bad_df[descriptor].median() - good_df[descriptor].median()
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
    file_name = "pvalues_good_bad_denovo.csv"
    if not TWO_SIDED:
        file_name = "pvalues_good_bad_denovo_one_sided.csv"
    results_df.to_csv(file_name, sep="\t", index=False)
    return results_df


def process_iteration(good_indexes, bad_indexes, descriptors_df, descriptors, medians):
    # Sample the data
    n_to_sample = min([len(good_indexes), len(bad_indexes)])
    good_sample = random.sample(good_indexes, n_to_sample)
    bad_sample = random.sample(bad_indexes, n_to_sample)
    good_df = descriptors_df.iloc[good_sample]
    bad_df = descriptors_df.iloc[bad_sample]

    # Get which way to do the test
    ways = None
    if not TWO_SIDED:
        ways = get_which_way(medians)
    
    # Get the median difference
    median_diff = get_median_diff(good_df, bad_df, ways)

    # Pool the data and get the median diff
    pool1, pool2 = pool_cdss(good_df, bad_df)
    random_diff = get_median_diff(pool1, pool2, ways)

    # Compare the medians
    result = compare_medians(median_diff, random_diff)
    return result



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
    if INCLUDE == "all":
        descriptors = list(descriptors_df)
        descriptors.remove("genome")
        descriptors.remove("type")
        descriptors.remove("cds")
    else:
        descriptors = INCLUDE
        # Check if all descriptors are in the dataframe
        for descriptor in descriptors:
            if descriptor not in descriptors_df.columns:
                raise ValueError(f"Descriptor {descriptor} not found in {descriptors_file}")
    
    
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
    signif_results = pd.DataFrame(columns=descriptors)

    # Calculate the global medians to know which way to compare
    medians = None
    if not TWO_SIDED:
        medians = get_all_medians(descriptors_df, good_indexes, bad_indexes, descriptors)

    # Use ProcessPoolExecutor for parallel processing
    print("Starting parallel processing...")
    with ProcessPoolExecutor() as executor:
        # Submit tasks to the executor
        futures = [executor.submit(process_iteration, good_indexes, bad_indexes, descriptors_df, descriptors, medians) for _ in range(n)]

        # Track progress
        completed = 0
        progress_step = n // 10  # 10% of total iterations
        
        # Collect results as they complete
        for future in futures:
            result = future.result()
            signif_results = pd.concat([signif_results, result], ignore_index=True)
            completed += 1

            # Print progress every 10%
            if completed % progress_step == 0:
                print(f"Progress: {completed}/{n} iterations completed ({(completed / n) * 100:.0f}%)")

    # Calculate p-values and save to file
    pvalues = calculate_pvalues(signif_results, descriptors)
    print("\nDone!")