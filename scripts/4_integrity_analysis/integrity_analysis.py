"""
This script, for each de novo gene, extracts all the non-coding
matches in the last and before-last outgroup. Then, it calculates
the global qcov and qcov of the longest ORF. The results are saved
in a file named "nc_coverage.tsv" in the OUT_FOLDER directory.
"""


import os
import re
import glob
import tempfile
import subprocess
import concurrent.futures
import math
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt



from my_functions.paths import *
OUT_FOLDER = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/4_integrity_analysis/"



def get_nc_placement(genome, denovo, nc_genome):
    """Get the placement of a non-coding match for a denovo gene in a genome"""
    # Read the blast_out file
    tblastn_result_file = f"{DENSE_DIR}{genome}/blast_out/TRG_multielongated_tblastn_{nc_genome}_genome.out"
    tblastn_df = pd.read_csv(tblastn_result_file, sep="\t", header=None, comment="#")
    tblastn_df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlength", "qcov", "sframe"]
    tblastn_df_denov = tblastn_df[tblastn_df["qseqid"] == f"{denovo}_elongated"].reset_index(drop=True)
    # Sort by qcov then evalue
    tblastn_df_denov = tblastn_df_denov[tblastn_df_denov["qcov"] >= 50]
    tblastn_df_denov = tblastn_df_denov.sort_values(by="evalue")
    # Extract the blast result
    strand = "+" if tblastn_df_denov.iloc[0]["sstart"] < tblastn_df_denov.iloc[0]["send"] else "-"
    contig = tblastn_df_denov.iloc[0]["sseqid"]
    start_loci = tblastn_df_denov.iloc[0]["sstart"] - 1 if strand == "+" else tblastn_df_denov.iloc[0]["send"] - 1
    end_loci = tblastn_df_denov.iloc[0]["send"] if strand == "+" else tblastn_df_denov.iloc[0]["sstart"]
    # Get the position and the qcov
    loci = (contig, int(start_loci), int(end_loci), strand)
    qcov = tblastn_df_denov.iloc[0]["qcov"]
    return loci, qcov


def get_denovo_info(genome):
    """Get the info on de novo genes for a given genome"""
    denovo_dict = {}
    denovo_file = os.path.join(DENSE_DIR, genome, "denovogenes.tsv")
    denovo_df = pd.read_csv(denovo_file, sep="\t", header=0)
    if denovo_df.empty:
        return {}

    # Get the name of de novo genes
    for denovo_gene in denovo_df["CDS"]:
        denovo_dict[denovo_gene] = {}

    # Get the de novo sequence
    fna_file = os.path.join(CDS_DIR, genome + "_CDS.faa")
    for denovo_gene in denovo_dict:
        for record in SeqIO.parse(fna_file, "fasta"):
            if record.name == denovo_gene:
                denovo_dict[denovo_gene]["sequence"] = record.seq

    # Get the NC hits
    trg_file = os.path.join(DENSE_DIR, genome, "TRG_match_matrix.tsv")
    trg_df = pd.read_csv(trg_file, sep="\t", header=0)
    for denovo in denovo_dict:
        denovo_dict[denovo]["ancestor"] = genome

        # Get the row corresponding to the denovo gene
        try:
            hit_row = trg_df[trg_df["CDS"] == f"{denovo}_elongated"]
            if hit_row.empty:
                raise ValueError(f"No corresponding row found for {denovo}_elongated in {genome}")
        except ValueError as e:
            print(e)
            continue
        # Flatten the row values to a list of strings
        row_values = hit_row.astype(str).values.flatten()
        # Get all the outgroup numbers
        gs_numbers = set()
        for value in row_values:
            matches = re.findall(r'gS(\d+)', value)
            gs_numbers.update(map(int, matches))
        gs_list = sorted(gs_numbers)
        # Add the last 2 outgroup numbers to the dict
        nc_matches_dict = {"last": {"nb": gs_list[-1], "values": {}},
                           "before_last": {"nb": gs_list[-2], "values": {}}}
        
        # get all the ancestor names
        for outgroup in nc_matches_dict:
            column_names = hit_row.columns.tolist()
            outgroup_number = nc_matches_dict[outgroup]["nb"]
            target = f"gS{outgroup_number}"
            indices = [i for i, val in enumerate(row_values) if val == target]

            for i in indices:
                nc_genome = column_names[i]
                # Get the placement of the nc match
                loci, qcov = get_nc_placement(genome, denovo, nc_genome)
                nc_matches_dict[outgroup]["values"][nc_genome] = {"loci": loci, "qcov": qcov}

        denovo_dict[denovo]["nc_matches"] = nc_matches_dict
    return denovo_dict


def get_sequence_from_loci(genome, contig, start, end, strand):
    fa_file = os.path.join(FA_DIR, genome + ".fa")

    for record in SeqIO.parse(fa_file, "fasta"):
        if str(record.name) == str(contig):
            seq = record.seq[start:end]
    return seq if strand == "+" else seq.reverse_complement()



if __name__ == "__main__":
    # Read list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]

    # Get the denovo info for each genome
    denovo_info = {}
    for genome in genomes:
        denovo_info_local = get_denovo_info(genome)
        if denovo_info_local == {}:
            continue
        denovo_info.update(denovo_info_local)
    
    qcov_dict = {}
    # For each denovo
    for denovo in denovo_info:
        qcov_dict[denovo] = {}
        # Get the denovo sequence
        denovo_seq = denovo_info[denovo]["sequence"]
        # For each species in each outgroup (last and ante-last)
        for outgroup in denovo_info[denovo]["nc_matches"]:
            qcov_dict[denovo][outgroup] = {}
            for neighbor_sp in denovo_info[denovo]["nc_matches"][outgroup]["values"]:
                # Get the nc match sequence
                nc_loci = denovo_info[denovo]["nc_matches"][outgroup]["values"][neighbor_sp]["loci"]
                nc_seq = str(get_sequence_from_loci(neighbor_sp, *nc_loci).translate())
                # Get the location of the stop codons
                stops = [pos for pos, char in enumerate(nc_seq) if char == "*"]
                # Get the longest orf
                match_stops = [0] + stops + [len(nc_seq) - 1]
                longest_orf = 0
                for i in range(len(match_stops) - 1):
                    orf_len = match_stops[i + 1] - match_stops[i] - 1
                    if orf_len > longest_orf:
                        longest_orf = orf_len

                # Get the coverages
                qcov = denovo_info[denovo]["nc_matches"][outgroup]["values"][neighbor_sp]["qcov"]
                qcov_orf = (longest_orf / len(denovo_seq)) * 100
                qcov_dict[denovo][outgroup][neighbor_sp] = {
                    "qcov": qcov,
                    "qcov_orf": qcov_orf
                }


    # Save the qcov_dict to a file
    qcov_file = os.path.join(OUT_FOLDER, "nc_coverage.tsv")
    with open(qcov_file, "w") as f:
        f.write("denovo\toutgroup\tneighbor_sp\tqcov\tqcov_orf\n")
        for denovo in qcov_dict:
            for outgroup in qcov_dict[denovo]:
                for neighbor_sp in qcov_dict[denovo][outgroup]:
                    qcov = qcov_dict[denovo][outgroup][neighbor_sp]["qcov"]
                    qcov_orf = qcov_dict[denovo][outgroup][neighbor_sp]["qcov_orf"]
                    f.write(f"{denovo}\t{outgroup}\t{neighbor_sp}\t{qcov}\t{qcov_orf}\n")
                


    