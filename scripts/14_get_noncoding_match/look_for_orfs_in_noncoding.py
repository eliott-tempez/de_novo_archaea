import os
import re
import glob
import tempfile
import subprocess
import concurrent.futures
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord




DENSE_DIR = "/home/eliott.tempez/Documents/archaea_data/dense/"
DATA_DIR = "/home/eliott.tempez/Documents/archaea_data/complete_122/"
GENOMES_LIST = "/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/6_genera_archaea/genomes_list.txt"
OUT_FOLDER = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/"



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
    fna_file = os.path.join(DATA_DIR, "CDS", genome + "_CDS.faa")
    for denovo_gene in denovo_dict:
        for record in SeqIO.parse(fna_file, "fasta"):
            if record.name == denovo_gene:
                denovo_dict[denovo_gene]["sequence"] = record.seq

    # Get the last NC hit in synteny
    ## Name of the ancestor
    trg_file = os.path.join(DENSE_DIR, genome, "TRG_match_matrix.tsv")
    trg_df = pd.read_csv(trg_file, sep="\t", header=0)
    for denovo in denovo_dict:
        try:
            matches = trg_df[trg_df["CDS"] == f"{denovo}_elongated"]
            if matches.empty:
                raise ValueError(f"No corresponding row found for {denovo}_elongated in {genome}")
        except ValueError as e:
            print(e)
            continue
        # Iterate over the row in reverse until we hit a cds
        i = len(matches.columns)
        cell = ""
        while "gS" not in cell:
            i -= 1
            cell = matches.iloc[0, i]
        denovo_dict[denovo]["ancestor_sp"] = matches.columns[i]
    
    ## Loci of the noncoding match
    unique_ancestors = set(denovo_dict[denovo]["ancestor_sp"] for denovo in denovo_dict)
    # Get the tblastn result for each ancestor
    for ancestor in unique_ancestors:
        tblastn_result_file = f"{DENSE_DIR}{genome}/blast_out/TRG_multielongated_tblastn_{ancestor}_genome.out"
        tblastn_df = pd.read_csv(tblastn_result_file, sep="\t", header=None, comment="#")
        tblastn_df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlength", "qcov", "sframe"]
        for denovo in denovo_dict:
            if denovo_dict[denovo]["ancestor_sp"] == ancestor:
                tblastn_df_denov = tblastn_df[tblastn_df["qseqid"] == f"{denovo}_elongated"].reset_index(drop=True)
                # Sort by qcov then evalue
                tblastn_df_denov = tblastn_df_denov[tblastn_df_denov["qcov"] >= 50]
                tblastn_df_denov = tblastn_df_denov.sort_values(by="evalue")
                # Extract the blast result
                strand = "+" if tblastn_df_denov.iloc[0]["sstart"] < tblastn_df_denov.iloc[0]["send"] else "-"
                contig = tblastn_df_denov.iloc[0]["sseqid"]
                start_loci = tblastn_df_denov.iloc[0]["sstart"] - 1 if strand == "+" else tblastn_df_denov.iloc[0]["send"] - 1
                end_loci = tblastn_df_denov.iloc[0]["send"] if strand == "+" else tblastn_df_denov.iloc[0]["sstart"]
                denovo_dict[denovo]["loci"] = [contig, int(start_loci), int(end_loci), strand]
                # Get the part of the de novo gene that matched
                denovo_dict[denovo]["qstart"] = tblastn_df_denov.iloc[0]["qstart"] - 1
                denovo_dict[denovo]["qend"] = tblastn_df_denov.iloc[0]["qend"]

    return denovo_dict




def get_sequence_from_loci(genome, contig, start, end, strand):
    fa_file = os.path.join(DATA_DIR, "fasta_renamed", genome + ".fa")

    for record in SeqIO.parse(fa_file, "fasta"):
        if str(record.name) == str(contig):
            seq = record.seq[start:end]
    return seq if strand == "+" else seq.reverse_complement()








if __name__ == "__main__":
    # Read list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]

    denovo_dict = {}
    n_potential_CDSs = 0
    n_denovo = 0
    good_denovo_stops = []
    for genome in genomes:
        # Get de novo info for each genome
        denovo_dict[genome] = get_denovo_info(genome)
        if denovo_dict[genome] == {}:
            continue

        # For each de novo gene, get the protein sequence of the match
        for denovo in denovo_dict[genome]:
            n_denovo += 1
            outgroup = denovo_dict[genome][denovo]["ancestor_sp"]
            contig, start, end, strand = denovo_dict[genome][denovo]["loci"]
            match_seq = get_sequence_from_loci(outgroup, contig, start, end, strand)
            match_seq_aa = str(match_seq.translate())

            # Get the location of the stop codons
            stops = [pos for pos, char in enumerate(match_seq_aa) if char == "*"]
            max_threshold_codon = len(match_seq_aa) - 4
            # Get the number of sequences that could be CDSs
            all_stops_after_threshold = all(stop > max_threshold_codon for stop in stops)
            no_stops = len(stops) == 0
            if all_stops_after_threshold or no_stops:
                n_potential_CDSs += 1
            else:
                good_denovo_stops.append(stops)


    print(f"Out of {n_denovo} de novo genes, {n_potential_CDSs} corresponding hits in the noncoding have no stop at all, or only stops in the last 3 amino acids")
    print(f"The rest have a mean of {np.mean([len(stops) for stops in good_denovo_stops])} stops in the sequence")
    print()

                

        










