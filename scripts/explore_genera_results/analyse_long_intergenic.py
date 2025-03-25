import os
import glob
from Bio import SeqIO
import re
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import metapredict as meta
import random


GENOMES_LIST = "/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/genera_archaea/genomes_list.txt"
DATA_DIR = "/home/eliott.tempez/Documents/archaea_data/complete_122"
GENERA_DIR = "/home/eliott.tempez/Documents/archaea_data/genera/out/"
OUT_DIR = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/explore_genera_results/intergenic_orfs_analysis"


def extract_intergenic(species):
    """Extract all intergenic sequences"""
    # Make sure the file exists and read it
    gff_file = os.path.join(DATA_DIR, "reannotated_gff_75/", species + ".gff3")
    gff = pd.read_csv(gff_file, sep="\t", header=None, comment="#")
    gff.columns = ["contig", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

    # Get length of genome
    with open(gff_file, "r") as f:
        line = f.readline()
        while "##sequence-region" not in line:
            line = next(f)
        genome_length = int(line.split()[3])

    # Extract data for each contig
    intergenic_segments = {}
    for record in gff["contig"].unique():
        record_str = str(record)
        intergenic_segments[record_str] = []
        # Get all nucl positions in CDSs from the gff file
        coding_loci = []
        for i, row in gff[gff["contig"] == record].iterrows():
            if row["type"] == "CDS":
                coding_loci += list(range(row["start"], row["end"]))
        # Get all the integers that are not in the list of coding loci
        full_set = set(range(min(coding_loci), max(coding_loci) + 1))
        noncoding_positions = list(full_set - set(coding_loci))
        # Convert to consecutive ranges
        noncoding_positions.sort()
        start = noncoding_positions[0]
        # Add first segment if intergenic
        if 0 not in set(coding_loci):
            intergenic_segments[record_str].append((0, min(coding_loci) - 1))
        # Add the rest of the segments
        for i in range(1, len(noncoding_positions)):
            if noncoding_positions[i] != noncoding_positions[i - 1] + 1:
                intergenic_segments[record_str].append((start, noncoding_positions[i - 1]))
                start = noncoding_positions[i]
        intergenic_segments[record_str].append((start, noncoding_positions[-1]))
        # Add last segment if intergenic
        if genome_length - 1 not in set(coding_loci):
            intergenic_segments[record_str].append((max(coding_loci) + 1, genome_length))


    # Extract the corresponding sequences from fasta file
    fa_file = os.path.join(DATA_DIR, "fasta_renamed/" + species + ".fa")
    fa_content = list(SeqIO.parse(fa_file, "fasta"))
    intergenic_dict = {}
    i = 0
    for contig in fa_content:
        if contig.name in intergenic_segments:
            for j in range(len(intergenic_segments[contig.name])):
                i += 1
                start, end = intergenic_segments[contig.name][j]
                intergenic_dict[f"{contig.name}_{i}"] = {"seq": contig.seq[start:end], "pos": (start, end), "contig": contig.name}
    return intergenic_dict


def extract_orfs(intergenic_dict):
    """Extract all ORFs from the intergenic sequences"""
    orf_dict = {}
    i = 0
    for inter in intergenic_dict:
        seq = intergenic_dict[inter]["seq"]
        start_inter = intergenic_dict[inter]["pos"][0]
        # Get all ORFs
        for frame in range(3):
            for i in range(frame, len(seq), 3):
                if i + 3 < len(seq):
                    codon = seq[i:i+3]
                    if codon == "ATG":
                        start = i
                        for j in range(i, len(seq), 3):
                            if j + 3 < len(seq):
                                codon = seq[j:j+3]
                                if codon in ["TAA", "TAG", "TGA"]:
                                    i += 1
                                    end = j + 3
                                    orf_dict[f"orf_{i}"] = {"inter": inter, "start": start + start_inter, "end": end + start_inter, "seq": seq[start:end]}
                                    break
    return orf_dict


def get_cds_seqs(genome):
    """Get the sequences of all CDSs in the genome"""
    seqs = []
    gff_file = os.path.join(DATA_DIR, "reannotated_gff_75/", genome + ".gff3")
    fasta_file = os.path.join(DATA_DIR, "fasta_renamed", genome + ".fa")
    gff = pd.read_csv(gff_file, sep="\t", header=None, comment="#")
    gff.columns = ["contig", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    fa_content = list(SeqIO.parse(fasta_file, "fasta"))
    for contig in gff["contig"].unique():
        for row in gff[gff["contig"] == contig].iterrows():
            if row[1]["type"] == "CDS":
                start = row[1]["start"]
                end = row[1]["end"]
                for record in fa_content:
                    if record.name == str(contig):
                        seq = record.seq[start:end]
                        seqs.append(seq)
                        break
    return seqs


def get_entropy(seq):
    """Calculate the entropy of a sequence"""
    seq = str(seq)
    freqs = {}
    for base in seq:
        if base in freqs:
            freqs[base] += 1
        else:
            freqs[base] = 1
    entropy = 0
    for base in freqs:
        p = freqs[base] / len(seq)
        entropy -= p * np.log2(p)
    return entropy


def get_hexamer_frequencies(sequences):
    """Calculate hexamer frequencies in given sequences."""
    hexamer_counts = defaultdict(int)
    total = 0

    for seq in sequences:
        seq = seq.upper()
        for i in range(len(seq) - 5):
            hexamer = seq[i:i+6]
            hexamer_counts[hexamer] += 1
            total += 1

    # Normalize frequencies
    hexamer_freq = {k: v / total for k, v in hexamer_counts.items()}
    return hexamer_freq


def calculate_log_odds_ratio(coding_freq, noncoding_freq):
    """Compute log-likelihood ratios for hexamers."""
    hexamers = set(coding_freq.keys()).union(noncoding_freq.keys())
    log_ratios = {}

    for hexamer in hexamers:
        fc = coding_freq.get(hexamer, 1e-6)  # Avoid zero frequency
        fn = noncoding_freq.get(hexamer, 1e-6)
        log_ratios[hexamer] = np.log2(fc / fn)

    return log_ratios


def score_sequence(seq, log_ratios):
    """Calculate hexamer score for a given sequence."""
    score = 0
    seq = seq.upper()
    for i in range(len(seq) - 5):
        hexamer = seq[i:i+6]
        if hexamer in log_ratios:
            score += log_ratios[hexamer]
    return score


def get_disordered_score(seq):
    """Get the mean disordered scores for a given sequence"""
    seq = str(seq)
    seq = re.sub(r"\*", "", seq)
    scores = meta.predict_disorder(seq)
    return np.mean(scores)




if __name__ == "__main__":
    hypotheticals = {}

    # Read the list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]


    # Chose a random genome to do the analysis on
    genome = random.sample(genomes, 1)[0]

    
    # Get the intergenic sequences
    intergenic_dict = extract_intergenic(genome)
    # Get the ORFs
    intergenic_orfs = extract_orfs(intergenic_dict)


    # Plot the distributions of CDSs and intergenic ORfs
    cds_seqs = get_cds_seqs(genome)
    cds_len = [len(seq) for seq in cds_seqs]
    orf_lens = [len(intergenic_orfs[orf]["seq"]) for orf in intergenic_orfs]
    first_quartile_len_cds = np.percentile(cds_len, 25)
    median_len_orfs = np.median(orf_lens)
    decile_len_cds = np.percentile(cds_len, 10)
    """# Calculate the common bin edges
    min_len = min(min(cds_len), min(orf_lens))
    max_len = max(max(cds_len), max(orf_lens))
    bins = np.linspace(min_len, max_len, 50)
    # Plot the histograms with the common bins
    plt.hist(orf_lens, bins=bins, alpha=1, label="intergenic ORFs")
    plt.hist(cds_len, bins=bins, alpha=0.6, label="CDSs")
    # Add a vertical line for first quartile
    plt.axvline(first_quartile_len_cds, color="red", linestyle="--", label="First quartile of CDS lengths")
    plt.legend()
    plt.xlabel("Length (bp)")
    plt.ylabel("Number")
    plt.title("Distribution of intergenic ORFs and CDS lengths")
    plt.show()


    # Entropy calculation
    cds_entropies = [get_entropy(seq) for seq in cds_seqs]
    orf_entropies = [get_entropy(intergenic_orfs[orf]["seq"]) for orf in intergenic_orfs]
    orf_entropies_small = [get_entropy(intergenic_orfs[orf]["seq"]) for orf in intergenic_orfs if len(intergenic_orfs[orf]["seq"]) < median_len_orfs]
    orf_entropies_long = [get_entropy(intergenic_orfs[orf]["seq"]) for orf in intergenic_orfs if len(intergenic_orfs[orf]["seq"]) >= median_len_orfs]
    # Calculate the common bin edges
    min_entropy = min(min(cds_entropies), min(orf_entropies))
    max_entropy = max(max(cds_entropies), max(orf_entropies))
    bins = np.linspace(min_entropy, max_entropy, 50)
    first_quartile_entropy_cds = np.percentile(cds_entropies, 25)
    # Plot the histograms with the common bins
    plt.hist(orf_entropies, bins=bins, alpha=1, label="intergenic ORFs")
    plt.hist(cds_entropies, bins=bins, alpha=0.6, label="CDSs")
    # Add a vertical line for first quartile
    plt.axvline(first_quartile_entropy_cds, color="red", linestyle="--", label="First quartile of CDS entropies")
    plt.legend()
    plt.xlabel("Shannon Entropy")
    plt.ylabel("Number")
    plt.title("Distribution of intergenic ORFs and CDS entropies")
    plt.show()

    # Split by size
    plt.hist(orf_entropies_small, bins=bins, alpha=0.5, label="Small intergenic ORFs")
    plt.hist(orf_entropies_long, bins=bins, alpha=0.5, label="Long intergenic ORFs")
    plt.legend()
    plt.xlabel("Shannon Entropy")
    plt.ylabel("Number")
    plt.title("Distribution of intergenic ORFs and CDS entropies")
    plt.show()


    # Get the interesting ORFs
    for orf in intergenic_orfs:
        entropy = get_entropy(intergenic_orfs[orf]["seq"])
        length = len(intergenic_orfs[orf]["seq"])
        intergenic_orfs[orf]["entropy"] = entropy
        if length > first_quartile_len_cds and entropy > first_quartile_entropy_cds:
            if genome not in hypotheticals:
                hypotheticals[genome] = [intergenic_orfs[orf]]
            else:
                hypotheticals[genome].append(intergenic_orfs[orf])

    # Write the number of hypothetical ORFs per genome
    if genome not in hypotheticals:
        print(f"Genome {genome} has 0 hypothetical ORFs")
    else:
        print(f"Genome {genome} has {len(hypotheticals[genome])} hypothetical ORFs")"""

    # ORF length
    cds_len = [len(seq) for seq in cds_seqs]
    orf_lens = [len(intergenic_orfs[orf]["seq"]) for orf in intergenic_orfs]
    plt.figure()
    data = [orf_lens, cds_len]
    box = plt.boxplot(data, patch_artist=True, labels=["intergenic ORFs", "CDSs"])
    colors = ['#009E73', '#E69F00']
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
    for median in box['medians']:
        median.set_color('black')
    plt.ylabel("Length (bp)")
    plt.title(f"Intergenic ORFs and CDSs lengths distribution\n(CDS length 1dt decile = {round(decile_len_cds)} nt)")
    # Save the plots
    os.makedirs(OUT_DIR, exist_ok=True)
    plt.savefig(os.path.join(OUT_DIR, f"{genome}_lengths.png"))


    ##### HETERRGENEITY METRICS ####
    ## Entropy of dna sequence
    orf_entropies_small = [get_entropy(intergenic_orfs[orf]["seq"]) for orf in intergenic_orfs if len(intergenic_orfs[orf]["seq"]) < decile_len_cds]
    orf_entropies_long = [get_entropy(intergenic_orfs[orf]["seq"]) for orf in intergenic_orfs if len(intergenic_orfs[orf]["seq"]) >= decile_len_cds]
    cds_entropies = [get_entropy(seq) for seq in cds_seqs]
    # Plot
    plt.figure()
    data = [orf_entropies_small, orf_entropies_long, cds_entropies]
    box = plt.boxplot(data, 
                      tick_labels=[f"Small intergenic ORFs\n(n = {len(orf_entropies_small)})", 
                                   f"Long intergenic ORFs\n(n = {len(orf_entropies_long)})", 
                                   f"CDSs\n(n = {len(cds_entropies)})"], 
                                   patch_artist=True)
    colors = ["#CC79A7", "#56B4E9", "#E69F00"]
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
    for median in box['medians']:
        median.set_color('black')
    plt.ylabel("Shannon Entropy")
    plt.title("Distribution of intergenic ORFs entropies (DNA sequence)")
    plt.savefig(os.path.join(OUT_DIR, f"{genome}_entropy.png"))


    ## Entropy of amino acid sequence
    orf_entropies_aa_small = [get_entropy(str(intergenic_orfs[orf]["seq"].translate())) for orf in intergenic_orfs if len(intergenic_orfs[orf]["seq"]) < decile_len_cds]
    orf_entropies_aa_long = [get_entropy(str(intergenic_orfs[orf]["seq"].translate())) for orf in intergenic_orfs if len(intergenic_orfs[orf]["seq"]) >= decile_len_cds]
    cds_entropies_aa = [get_entropy(str(seq.translate())) for seq in cds_seqs]
    # Plot
    plt.figure()
    data = [orf_entropies_aa_small, orf_entropies_aa_long, cds_entropies_aa]
    box = plt.boxplot(data, tick_labels=[f"Small intergenic ORFs\n(n = {len(orf_entropies_aa_small)})", 
                                   f"Long intergenic ORFs\n(n = {len(orf_entropies_aa_long)})", 
                                   f"CDSs\n(n = {len(cds_entropies_aa)})"],
                                   patch_artist=True)
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
    for median in box['medians']:
        median.set_color('black')
    plt.ylabel("Shannon Entropy")
    plt.title("Distribution of intergenic ORFs entropies (amino acid sequence)")
    plt.savefig(os.path.join(OUT_DIR, f"{genome}_entropy_aa.png"))

    ## Hexameric score
    # Get the hexamer frequencies
    cds_freq = get_hexamer_frequencies(cds_seqs)
    intergenic_freq = get_hexamer_frequencies([intergenic_dict[inter]["seq"] for inter in intergenic_dict])
    # Compute log-odds hexamer model
    hexamer_scores = calculate_log_odds_ratio(cds_freq, intergenic_freq)
    # Calculate the scores
    cds_scores = [score_sequence(seq, hexamer_scores) for seq in cds_seqs]
    orf_scores_small = [score_sequence(intergenic_orfs[orf]["seq"], hexamer_scores) for orf in intergenic_orfs if len(intergenic_orfs[orf]["seq"]) < decile_len_cds]
    orf_scores_long = [score_sequence(intergenic_orfs[orf]["seq"], hexamer_scores) for orf in intergenic_orfs if len(intergenic_orfs[orf]["seq"]) >= decile_len_cds]
    # Plot
    plt.figure()
    data = [orf_scores_small, orf_scores_long, cds_scores]
    box = plt.boxplot(data, tick_labels=[f"Small intergenic ORFs\n(n = {len(orf_scores_small)})",
                                   f"Long intergenic ORFs\n(n = {len(orf_scores_long)})",
                                   f"CDSs\n(n = {len(cds_scores)})"],
                                   patch_artist=True)
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
    for median in box['medians']:
        median.set_color('black')
    plt.ylabel("Hexamer score")
    plt.title("Distribution of hexamer scores")
    plt.savefig(os.path.join(OUT_DIR, f"{genome}_hexamer.png"))
