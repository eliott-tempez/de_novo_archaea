import os
import glob
from Bio import SeqIO
import re
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


GENOMES_LIST = "/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/genera_archaea/genomes_list.txt"
DATA_DIR = "/home/eliott.tempez/Documents/archaea_data/complete_122"
GENERA_DIR = "/home/eliott.tempez/Documents/archaea_data/genera/out/"


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



if __name__ == "__main__":
    hypotheticals = {}

    # Read the list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]
    for genome in genomes:
        # Get the intergenic sequences
        intergenic_dict = extract_intergenic(genome)
        # Get the ORFs
        intergenic_orfs = extract_orfs(intergenic_dict)


        # Plot the distributions of CDSs and intergenic ORfs
        cds_seqs = get_cds_seqs(genome)
        cds_len = [len(seq) for seq in cds_seqs]
        orf_lens = [len(intergenic_orfs[orf]["seq"]) for orf in intergenic_orfs]
        first_quartile_len_cds = np.percentile(cds_len, 25)
        # Calculate the common bin edges
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
        #plt.show()


        # Entropy calculation
        cds_entropies = [get_entropy(seq) for seq in cds_seqs]
        orf_entropies = [get_entropy(intergenic_orfs[orf]["seq"]) for orf in intergenic_orfs]
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
        #plt.show()


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
            print(f"Genome {genome} has {len(hypotheticals[genome])} hypothetical ORFs")

    
