import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction as GC


from my_functions.paths import GENOMES_LIST, FA_DIR


NB_GC_BINS = 3


def get_species_gc_content(genome):
    fa_file = os.path.join(FA_DIR, genome + ".fa")
    if not os.path.exists(fa_file):
        raise FileNotFoundError(f"No file {fa_file}")
    seq = ""
    for record in SeqIO.parse(fa_file, "fasta"):
        seq += str(record.seq)
    return GC(seq)


def binned_indexes(indexes, descriptors_df, gc_dict, bin_limits):
    binned_indexes = []
    df = descriptors_df.iloc[indexes]
    print(df)
    for i in range(len(bin_limits) - 1):
        bin_range = [bin_limits[i], bin_limits[i + 1]]
        genomes_in_bin = [g for g in gc_dict if gc_dict[g] >= bin_range[0] and gc_dict[g] < bin_range[1]]
        print(bin_range)
        print([(g, gc_dict[g]) for g in genomes_in_bin])
        bin_indexes = df[df["genome"].isin(genomes_in_bin)].index.tolist()
        print(df[df["genome"].isin(genomes_in_bin)])
        binned_indexes.append(bin_indexes)
    
    return binned_indexes





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
    descriptors_file = "sequence_features_good_candidates_all.tsv"
    descriptors_df = pd.read_csv(descriptors_file, sep="\t", header=0)
    
    # Get the indexes for each type of cds
    denovo_indexes = list(descriptors_df[descriptors_df["type"] == "denovo"].index)
    trg_indexes = list(descriptors_df[descriptors_df["type"] == "trg"].index)
    cds_indexes = list(descriptors_df[descriptors_df["type"] == "cds"].index)
    
    # Get the bins ranges
    bin_size = (max_gc - min_gc) / NB_GC_BINS
    bin_limits = [0] + [i * bin_size + min_gc for i in range(NB_GC_BINS + 1)][1:-1] + [1]
    
    # Get the de novo bins
    denovo_bins_indexes = binned_indexes(denovo_indexes, gc_dict, bin_limits)
        
        