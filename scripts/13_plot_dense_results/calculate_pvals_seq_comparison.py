import os
import sys
import re
import random
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction as GC

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
    sampled_indexes_1 = random.sample(range(n), n/2)
    sampled_indexes_2 = [i for i in range(n) if i not in sampled_indexes_1]
    pooled_df_1 = pooled_df.iloc[sampled_indexes_1]
    pooled_df_2 = pooled_df.iloc[sampled_indexes_2]
    return pooled_df_1, pooled_df_2


def compare_medians(median_diff, random_diff, type1, type2, n1, n2, bin1, bin2):
    # Compare the median differences
    results = {}
    for descriptor in median_diff:
        if median_diff[descriptor] > random_diff[descriptor]:
            median_diff[descriptor] = 1
        else:
            median_diff[descriptor] = 0
        results[descriptor] = median_diff[descriptor]
    # Create a dataframe with the results
    results_lst = [type1, type2, n1, n2, bin1, bin2] + list(results.values())
    results_df = pd.DataFrame(results_lst, columns=["type1", "type2", "n1", "n2", "bin1", "bin2"] + list(results.keys()))
    return results_df





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
    descriptors.remove("genome", "cds", "type")
    
    # Get the indexes for each type of cds
    denovo_indexes = list(descriptors_df[descriptors_df["type"] == "denovo"].index)
    trg_indexes = list(descriptors_df[descriptors_df["type"] == "trg"].index)
    cds_indexes = list(descriptors_df[descriptors_df["type"] == "cds"].index)
    
    # Get the bins ranges
    bin_size = (max_gc - min_gc) / NB_GC_BINS
    bin_limits = [0] + [i * bin_size + min_gc for i in range(NB_GC_BINS + 1)][1:-1] + [1]
    
    # Get the different bin indexes
    bin_indexes = get_bin_indexes(descriptors_df, gc_dict, bin_limits)
    print([len(bin) for bin in bin_indexes])
    # Get the bin indexes for each type of cds
    denovo_bin_indexes = [list(set(denovo_indexes) & set(bin)) for bin in bin_indexes]
    trg_bin_indexes = [list(set(trg_indexes) & set(bin)) for bin in bin_indexes]
    cds_bin_indexes = [list(set(cds_indexes) & set(bin)) for bin in bin_indexes]
    
    # Number of iterations
    n = 1
    n_bins = len(bin_indexes)
    signif_columns = ["type1", "type2", "n1", "n2", "bin1", "bin2"] + descriptors
    denovo_trg_signif_results = pd.DataFrame(columns=signif_columns)
    denovo_cds_signif_results = pd.DataFrame(columns=signif_columns)
    trg_cds_signif_results = pd.DataFrame(columns=signif_columns)
    for i in range(n):
        # For each bin combination
        for i_bin_denovo in range(n_bins):
            sampled_denovo_indexes = denovo_bin_indexes[i_bin_denovo]
            n_to_sample = len(sampled_denovo_indexes)
            # Sample the trgs and de novo
            for i_bin_trg in range(n_bins):
                sampled_trg_indexes = random.sample(trg_bin_indexes[i_bin_trg], n_to_sample)
                for i_bin_cds in range(n_bins):
                    sampled_cds_indexes = random.sample(cds_bin_indexes[i_bin_cds], n_to_sample)
                    # Get the dfs
                    sampled_denovo_df = descriptors_df.iloc[sampled_denovo_indexes]
                    sampled_trg_df = descriptors_df.iloc[sampled_trg_indexes]
                    sampled_cds_df = descriptors_df.iloc[sampled_cds_indexes]
                    
                    # Get the median differences
                    denovo_trg_diff = get_median_diff(sampled_denovo_df, sampled_trg_df)
                    denovo_cds_diff = get_median_diff(sampled_denovo_df, sampled_cds_df)
                    trg_cds_diff = get_median_diff(sampled_trg_df, sampled_cds_df)
                    
                    # Pool samples and get the random median differences
                    denovo_trg_pool_1, denovo_trg_pool_2 = pool_cdss(sampled_denovo_df, sampled_trg_df)
                    denovo_cds_pool_1, denovo_cds_pool_2 = pool_cdss(sampled_denovo_df, sampled_cds_df)
                    trg_cds_pool_1, trg_cds_pool_2 = pool_cdss(sampled_trg_df, sampled_cds_df)
                    denovo_trg_random_diff = get_median_diff(denovo_trg_pool_1, denovo_trg_pool_2)
                    denovo_cds_random_diff = get_median_diff(denovo_cds_pool_1, denovo_cds_pool_2)
                    trg_cds_random_diff = get_median_diff(trg_cds_pool_1, trg_cds_pool_2)
                    
                    # Compare medians
                    denovo_trg_signif = compare_medians(denovo_trg_diff, denovo_trg_random_diff, "denovo", "trg", len(sampled_denovo_df), len(sampled_trg_df), i_bin_denovo, i_bin_trg)
                    denovo_cds_signif = compare_medians(denovo_cds_diff, denovo_cds_random_diff, "denovo", "cds", len(sampled_denovo_df), len(sampled_cds_df), i_bin_denovo, i_bin_cds)
                    trg_cds_signif = compare_medians(trg_cds_diff, trg_cds_random_diff, "trg", "cds", len(sampled_trg_df), len(sampled_cds_df), i_bin_trg, i_bin_cds)
                    # Add the results to the dataframes
                    denovo_trg_signif_results = pd.merge(denovo_trg_signif_results, denovo_trg_signif, how="outer")
                    print(denovo_trg_signif_results, "\n")
                    denovo_cds_signif_results = pd.merge(denovo_cds_signif_results, denovo_cds_signif, how="outer")
                    print(denovo_cds_signif_results, "\n")
                    trg_cds_signif_results = pd.merge(trg_cds_signif_results, trg_cds_signif, how="outer")
                    print(trg_cds_signif_results, "\n")
                    
                
            
            
            
    
            
            
        
        
    
        
        