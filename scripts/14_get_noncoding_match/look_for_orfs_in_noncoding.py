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
ORIGIN_FILE = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/denovo_noncoding_status.tsv"
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
    fna_file = os.path.join(CDS_DIR, genome + "_CDS.faa")
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
                denovo_dict[denovo]["qcov"] = tblastn_df_denov.iloc[0]["qcov"]

    return denovo_dict



def get_sequence_from_loci(genome, contig, start, end, strand):
    fa_file = os.path.join(FA_DIR, genome + ".fa")

    for record in SeqIO.parse(fa_file, "fasta"):
        if str(record.name) == str(contig):
            seq = record.seq[start:end]
    return seq if strand == "+" else seq.reverse_complement()



def parse_origin_file():
    origin_dict = {}
    origin_data = pd.read_csv(ORIGIN_FILE, sep="\t", header=0)
    # Replace by 0 the columns that have less than 3 codons match
    origin_data[["intergenic", "f+0", "f+1", "f+2"]] = origin_data[["intergenic", "f+0", "f+1", "f+2"]].replace(list(range(1, 10)), 0)
    for row in origin_data.iterrows():
        denovo = row[1]["denovo_gene"]
        intergenic = int(row[1]["intergenic"])
        f0 = int(row[1]["f+0"])
        f1 = int(row[1]["f+1"])
        f2 = int(row[1]["f+2"])
        origin_dict[denovo] = {"intergenic": intergenic, "f0": f0, "f1": f1, "f2": f2}
    return origin_dict



def get_orf_cov(denovo_dict, gene_list):
    orf_dict = {}
    for denovo in gene_list:
        # Get the info
        outgroup = denovo_dict[denovo]["ancestor_sp"]
        contig, start, end, strand = denovo_dict[denovo]["loci"]
        match_seq = get_sequence_from_loci(outgroup, contig, start, end, strand)
        match_seq_aa = str(match_seq.translate(table=11))
        query_aa_seq = str(denovo_dict[denovo]["sequence"])

        # Get the location of the stop codons
        stops = [pos for pos, char in enumerate(match_seq_aa) if char == "*"]
        orf_dict[denovo] = {"seq": match_seq, "seq_aa": match_seq_aa, "stops": stops}

        # Get the longest orf
        match_stops = [0] + stops + [len(match_seq_aa) - 1]
        longest_orf = 0
        for i in range(len(match_stops) - 1):
            orf_len = match_stops[i + 1] - match_stops[i] - 1
            if orf_len > longest_orf:
                longest_orf = orf_len
        orf_dict[denovo]["longest_orf"] = longest_orf

        # Get the coverage
        coverage = (longest_orf / len(query_aa_seq)) * 100
        orf_dict[denovo]["orf_qcov"] = round(coverage, 2)

    return orf_dict






if __name__ == "__main__":
    # Read list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]
    # Get origin for each de novo
    origin_dict = parse_origin_file()

    n_denovo, n_denovo_intergenic = 0, 0
    denovo_dict = {}
    intergenic_denovo = []
    for genome in genomes:

        # Get de novo info for each genome
        denovo_info = get_denovo_info(genome)
        if denovo_info == {}:
            continue

        for denovo in denovo_info:
            n_denovo += 1
            denovo_dict[denovo] = denovo_info[denovo]
            denovo_dict[denovo]["genome"] = genome

            # Origin
            n_intergenic = origin_dict[denovo]["intergenic"]
            n_f0 = origin_dict[denovo]["f0"]
            n_f1 = origin_dict[denovo]["f1"]
            n_f2 = origin_dict[denovo]["f2"]

            # Get de novo intergenic genes
            if n_intergenic > 0 and [n_f0, n_f1, n_f2] == [0, 0, 0]:
                n_denovo_intergenic += 1
                intergenic_denovo.append(denovo)
    

    # Get the genes with a qcov < 70
    low_qcov_genes_intergenic = [denovo for denovo in intergenic_denovo if denovo_dict[denovo]["qcov"] < 70]
    low_qcov_genes_rest = [denovo for denovo in denovo_dict if denovo not in intergenic_denovo and denovo_dict[denovo]["qcov"] < 70]
    print(f"\nOut of the {n_denovo} de novo genes, there are {n_denovo_intergenic} coming from an intergenic area and {n_denovo - n_denovo_intergenic} coming either from an altframe or a mix of the 2.\n")
    print(f"For the intergenic ones, {len(low_qcov_genes_intergenic)} have a qcov < 70%.")
    print(f"For the rest, {len(low_qcov_genes_rest)} have a qcov < 70%.")

    show_plots = input("\nDo you want to see the plots? (y/n) ")
    if show_plots == "y":
        qcov_intergenic = [denovo_dict[denovo]["qcov"] for denovo in intergenic_denovo]
        qcov_rest = [denovo_dict[denovo]["qcov"] for denovo in denovo_dict if denovo not in intergenic_denovo]

        # Get the common bin edges
        bin_edges = np.linspace(0, 100, 20)

        # Histograms
        plt.figure()
        plt.hist(qcov_intergenic, edgecolor = "black", bins=bin_edges)
        plt.xlim([0, 100])
        plt.ylim([0, 60])
        plt.vlines(x = 70, ymin = 0, ymax = 60, color = "red", linestyle = "dashed")
        plt.xlabel("Coverage of query sequence (%)")
        plt.ylabel("Number of hits")
        plt.title(f"Coverage of de novo genes by the tblastn hits\nfor the intergenic genes (n = {len(intergenic_denovo)})")

        plt.figure()
        plt.hist(qcov_rest, edgecolor = "black", bins=bin_edges)
        plt.xlim([0, 100])
        plt.ylim([0, 10])
        plt.vlines(x = 70, ymin = 0, ymax = 10, color = "red", linestyle = "dashed")
        plt.xlabel("Coverage of query sequence (%)")
        plt.ylabel("Number of hits")
        plt.title(f"Coverage of de novo genes by the tblastn hits\nfor the rest (n = {len(denovo_dict) - len(intergenic_denovo)})")
        plt.show()


    # Get the genes with an ORF cov < 70
    # Get the sequences and stop positions for the genes of interest
    orf_cov_dict_intergenic = get_orf_cov(denovo_dict, [denovo for denovo in intergenic_denovo if denovo not in low_qcov_genes_intergenic])
    orf_cov_dict_rest = get_orf_cov(denovo_dict, [denovo for denovo in denovo_dict if denovo not in intergenic_denovo and denovo not in low_qcov_genes_rest])
    # Get the ORF covs < 70
    low_orf_cov_genes_intergenic = [denovo for denovo in orf_cov_dict_intergenic if orf_cov_dict_intergenic[denovo]["orf_qcov"] < 70]
    low_orf_cov_genes_rest = [denovo for denovo in orf_cov_dict_rest if orf_cov_dict_rest[denovo]["orf_qcov"] < 70]
    # Print result
    print(f"\n\n\n\nFor the matches with a query coverage > 70%, we extracted the longuest ORF and calculated its coverage of the query (the de novo gene sequence)")
    print(f"\nFor the intergenic genes, {len(low_orf_cov_genes_intergenic)}/{n_denovo_intergenic - len(low_qcov_genes_intergenic)} have a longuest ORF coverage < 70%.")
    print(f"For the rest, {len(low_orf_cov_genes_rest)}/{(n_denovo - n_denovo_intergenic) - len(low_qcov_genes_rest)} have a longuest ORF coverage < 70%.")

    show_plots = input("\nDo you want to see the plots? (y/n) ")
    if show_plots == "y":
        orf_cov_intergenic = [orf_cov_dict_intergenic[denovo]["orf_qcov"] for denovo in orf_cov_dict_intergenic]
        orf_cov_rest = [orf_cov_dict_rest[denovo]["orf_qcov"] for denovo in orf_cov_dict_rest]

        # Get the common bin edges
        bin_edges = np.linspace(0, 100, 20)

        # Histograms
        plt.figure()
        plt.hist(orf_cov_intergenic, edgecolor = "black", bins=bin_edges)
        plt.xlim([0, 100])
        plt.ylim([0, 60])
        plt.vlines(x = 70, ymin = 0, ymax = 60, color = "red", linestyle = "dashed")
        plt.xlabel("Coverage of query sequence (%)")
        plt.ylabel("Number of ORFs")
        plt.title(f"Coverage of de novo genes by the longest ORF in the tblastn hits\nfor the intergenic genes (n = {len(orf_cov_dict_intergenic)})")

        plt.figure()
        plt.hist(orf_cov_rest, edgecolor = "black", bins=bin_edges)
        plt.xlim([0, 100])
        plt.ylim([0, 10])
        plt.vlines(x = 70, ymin = 0, ymax = 10, color = "red", linestyle = "dashed")
        plt.xlabel("Coverage of query sequence (%)")
        plt.ylabel("Number of ORFs")
        plt.title(f"Coverage of de novo genes by the longest ORF in the tblastn hits\nfor the rest (n = {len(orf_cov_dict_rest)})")
        plt.show()
    

    # Conclusion
    print(f"\n\n\n\nIn conclusion, we have {len(low_qcov_genes_intergenic) + len(low_orf_cov_genes_intergenic)} intergenic de novo genes with a qcov < 70% or an ORF coverage < 70%.")
    print(f"For the rest, we have {len(low_qcov_genes_rest) + len(low_orf_cov_genes_rest)} de novo genes with a qcov < 70% or an ORF coverage < 70%.\n")
    print(f"\n-> We have {len(low_qcov_genes_intergenic) + len(low_orf_cov_genes_intergenic) + len(low_qcov_genes_rest) + len(low_orf_cov_genes_rest)} good de novo gene candidates\n")

    # Export good candidates
    good_candidates = low_qcov_genes_intergenic + low_orf_cov_genes_intergenic + low_qcov_genes_rest + low_orf_cov_genes_rest
    good_candidates = list(set(good_candidates))
    with open(os.path.join(OUT_FOLDER, "good_candidates.txt"), "w") as f:
        for candidate in good_candidates:
            f.write(f"{candidate}\n")









