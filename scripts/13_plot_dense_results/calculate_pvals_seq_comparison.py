import os
import sys
import re
import random
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction as GC
from itertools import combinations

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from my_functions.paths import GENOMES_LIST, FA_DIR


NB_GC_BINS = 2


def get_species_gc_content(genome):
    fa_file = os.path.join(FA_DIR, genome + ".fa")
    if not os.path.exists(fa_file):
        raise FileNotFoundError(f"No file {fa_file}")
    seq = ""
    for record in SeqIO.parse(fa_file, "fasta"):
        seq += str(record.seq)
    return GC(seq)


def get_bin_indexes(descriptors_df, gc_dict, bin_limits):
    binned_indexes = []
    for i in range(len(bin_limits) - 1):
        bin_range = [bin_limits[i], bin_limits[i + 1]]
        genomes_in_bin = [g for g in gc_dict if gc_dict[g] >= bin_range[0] and gc_dict[g] < bin_range[1]]
        bin_indexes = descriptors_df[descriptors_df["genome"].isin(genomes_in_bin)].index.tolist()
        binned_indexes.append(bin_indexes)
    return binned_indexes


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


def compare_medians(median_diff, random_diff, type1, type2, bin1, bin2):
    # Compare the median differences
    results = {}
    for descriptor in median_diff:
        if median_diff[descriptor] > random_diff[descriptor]:
            results[descriptor] = 1
        else:
            results[descriptor] = 0
    # Create a dataframe with the results
    results_lst = [type1, type2, bin1, bin2] + list(results.values())
    results_df = pd.DataFrame([results_lst], columns=["type1", "type2", "bin1", "bin2"] + list(results.keys()))
    return results_df


def calculate_pvalues(signif_results, descriptors):
    print(signif_results[["type1", "type2", "bin1", "bin2"]].value_counts().reset_index(name='count'))
    



if __name__ == "__main__":
    # Get the list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]
    
    # Extract gc content
    gc_dict = {}
    min_gc = 1
    max_gc = 0
    for genome in genomes:
        # Get the species gc content
        genome_gc = get_species_gc_content(genome)
        gc_dict[genome] = genome_gc
        if genome_gc < min_gc:
            min_gc = genome_gc
        if genome_gc > max_gc:
            max_gc = genome_gc
    
    # Read descriptors 
    descriptors_file = "sequence_features_good_candidates_all.csv"
    descriptors_df = pd.read_csv(descriptors_file, sep="\t", header=0)
    descriptors = list(descriptors_df)
    descriptors.remove("genome")
    descriptors.remove("type")
    descriptors.remove("cds")
    
    # Get the indexes for each type of cds
    indexes = {}
    for type in ["denovo", "trg", "cds"]:
        indexes[type] = list(descriptors_df[descriptors_df["type"] == type].index)

    # Get the bins ranges
    bin_size = (max_gc - min_gc) / NB_GC_BINS
    bin_limits = [0] + [i * bin_size + min_gc for i in range(NB_GC_BINS + 1)][1:-1] + [1]
    
    # Get the different bin indexes
    bin_indexes_lst = get_bin_indexes(descriptors_df, gc_dict, bin_limits)
    # Get the bin indexes for each type of cds
    bin_indexes = {}
    for type in indexes:
        bin_indexes[type] = [list(set(indexes[type]) & set(bin)) for bin in bin_indexes_lst]
    
    # Number of iterations
    n = 10
    n_bins = len(bin_indexes["denovo"])
    n_to_sample = min([len(bin_indexes["denovo"][bin]) for bin in range(n_bins)])
    signif_columns = ["type1", "type2", "bin1", "bin2"] + descriptors
    signif_results = pd.DataFrame(columns=signif_columns)

    for i in range(n):
        # Sample each type for each bin
        samples = {}
        for type in bin_indexes:
            for bin in range(n_bins):
                key = (type, bin)
                draw = random.sample(bin_indexes[type][bin], n_to_sample)
                samples[key] = draw
                
        # Iterate on each possible combination
        for conf1, conf2 in combinations(samples.keys(), 2):
            type1, bin1 = conf1
            type2, bin2 = conf2
            
            # Get values
            sampled_type1_df = descriptors_df.iloc[samples[conf1]]
            sampled_type2_df = descriptors_df.iloc[samples[conf2]]
            # Get the median differences
            median_diff = get_median_diff(sampled_type1_df, sampled_type2_df)
            # Pool samples and get the random median differences
            pooled_df_1, pooled_df_2 = pool_cdss(sampled_type1_df, sampled_type2_df)
            random_diff = get_median_diff(pooled_df_1, pooled_df_2)
            # Compare medians
            signif = compare_medians(median_diff, random_diff, type1, type2, bin1, bin2)
            # Add the results to the dataframe
            signif_results = pd.merge(signif_results, signif, how="outer")
            
    pvalues = calculate_pvalues(signif_results, descriptors)
                

        

        
    
        
        