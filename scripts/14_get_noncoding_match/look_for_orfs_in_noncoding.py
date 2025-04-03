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



def get_good_candidates(denovo_dict, origin_dict, gene_list, threshold=5):
    no_stops, only_final_stop, stops_at_end = 0, 0, 0
    good_candidates = {}

    for denovo in gene_list:
        # Get the info
        outgroup = denovo_dict[denovo]["ancestor_sp"]
        contig, start, end, strand = denovo_dict[denovo]["loci"]
        match_seq = get_sequence_from_loci(outgroup, contig, start, end, strand)
        match_seq_aa = str(match_seq.translate(table=11))

        # Get the location of the stop codons
        stops = [pos for pos, char in enumerate(match_seq_aa) if char == "*"]

        # Get the bad and good candidates
        if len(stops) == 0:
            no_stops += 1
        elif stops == [len(match_seq_aa) - 1]:
            only_final_stop += 1
        elif min(stops) >= len(match_seq_aa) - threshold:
            stops_at_end += 1
        else:
            good_candidates[denovo] = {"seq": match_seq, "seq_aa": match_seq_aa, "stops": stops}

    return good_candidates, no_stops, only_final_stop, stops_at_end






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


    # Get the good candidates for the intergenic genes
    good_candidates_intergenic, no_stops_intergenic, only_final_stop_intergenic, stops_at_end_intergenic = get_good_candidates(denovo_dict, origin_dict, intergenic_denovo)
    # Get the good candidates for the rest
    good_candidates_rest, no_stops_rest, only_final_stop_rest, stops_at_end_rest = get_good_candidates(denovo_dict, origin_dict, [nov for nov in denovo_dict if nov not in intergenic_denovo])


            
    ############ Print results ############
    print(f"\nOut of {n_denovo} de novo genes, there are {n_denovo_intergenic} coming from an intergenic area:")
    print(f"   -> {no_stops_intergenic} have no stop codon")
    print(f"   -> {only_final_stop_intergenic} have only one stop codon at the end")
    print(f"   -> {stops_at_end_intergenic} have stop codon(s) only in the last 5 codons")
    print(f"\nWhich means {len(good_candidates_intergenic)} are good candidates")

    # Global query coverage of the tblast hit
    qcovs = [denovo_dict[denovo]["qcov"] for denovo in good_candidates_intergenic]

    # Coverage of the longuest ORF
    qcovs_orf = []
    good_good_candidates_intergenic = {}
    for denovo in good_candidates_intergenic:
        query_aa_seq_len = len(denovo_dict[denovo]["sequence"])
        aa_seq_len = len(good_candidates_intergenic[denovo]["seq_aa"])
        match_stops = [0] + good_candidates_intergenic[denovo]["stops"] + [aa_seq_len - 1]
        # Get the length of the longest ORF
        longest_orf = 0
        for i in range(len(match_stops) - 1):
            orf_len = match_stops[i + 1] - match_stops[i] - 1
            if orf_len > longest_orf:
                longest_orf = orf_len
        # Get the coverage
        coverage = (longest_orf / query_aa_seq_len) * 100
        good_candidates_intergenic[denovo]["longest_orf"] = longest_orf
        good_candidates_intergenic[denovo]["orf_qcov"] = round(coverage, 2)
        qcovs_orf.append(coverage)
        # Keep only < 70
        if coverage < 70:
            good_good_candidates_intergenic[denovo] = good_candidates_intergenic[denovo]
    
    # Show the results
    print(f"\n\nOut of {len(good_candidates_intergenic)} good candidates, there are {len(good_good_candidates_intergenic)} for which the longest ORF covers < 70% of the de novo gene")
    show_plots = input("\nDo you want to see the plots? (y/n) ")
    if show_plots == "y":
        # Get the common bin edges
        bin_edges = np.linspace(0, 100, 20)

        # Histograms
        plt.figure()
        plt.hist(qcovs, edgecolor = "black", bins=bin_edges)
        plt.xlim([0, 100])
        plt.ylim([0, 26])
        plt.xlabel("Coverage of query sequence (%)")
        plt.ylabel("Number of hits")
        plt.title(f"Coverage of de novo genes by the tblastn hits\nfor the good intergenic candidates (n = {len(good_candidates_intergenic)})")

        plt.figure()
        plt.hist(qcovs_orf, edgecolor = "black", bins=bin_edges)
        plt.xlim([0, 100])
        plt.ylim([0, 26])
        plt.vlines(x = 70, ymin = 0, ymax = 30, color = "red", linestyle = "dashed")
        plt.xlabel("Coverage of query sequence (%)")
        plt.ylabel("Number of ORFs")
        plt.title(f"Coverage of de novo genes by the longest ORF in the tblastn hits\nfor the good intergenic candidates (n = {len(good_candidates_intergenic)})")
        plt.show()


    # Display ORFs with stops
    print(f"\n\nHere are where the stops are for the {len(good_good_candidates_intergenic)} chosen candidates:\t\t\tLongest ORF length / de novo gene length (coverage)\n")
    for denovo in good_good_candidates_intergenic:
        genome = denovo_dict[denovo]["genome"]
        seq = good_candidates_intergenic[denovo]["seq_aa"]
        stops = good_candidates_intergenic[denovo]["stops"]
        relative_stops = [math.floor((s/len(seq))*100) for s in stops]
        orf_lst = ["-"] * 100
        for i in relative_stops:
            orf_lst[i] = "*"
        orf_str = "".join(orf_lst)
        print(f"{orf_str}   {good_candidates_intergenic[denovo]["longest_orf"]} / {len(denovo_dict[denovo]["sequence"])} ({good_candidates_intergenic[denovo]["orf_qcov"]}%)\n")

    
    # For the rest
    print(f"For the {n_denovo - n_denovo_intergenic} de novo genes that are not 100% intergenic, there are:")
    print(f"   -> {no_stops_rest} that have no stop codon")
    print(f"   -> {only_final_stop_rest} that have only one stop codon at the end")
    print(f"   -> {stops_at_end_rest} that have stop codon(s) only in the last 5 codons")
    print(f"\nWhich means {len(good_candidates_rest)} are good candidates")

    # Global query coverage of the tblast hit
    qcovs_rest = [denovo_dict[denovo]["qcov"] for denovo in good_candidates_rest]

    # Coverage of the longuest ORF
    qcovs_orf_rest = []
    good_good_candidates_rest = {}
    for denovo in good_candidates_rest:
        query_aa_seq_len = len(denovo_dict[denovo]["sequence"])
        aa_seq_len = len(good_candidates_rest[denovo]["seq_aa"])
        match_stops = [0] + good_candidates_rest[denovo]["stops"] + [aa_seq_len - 1]
        # Get the length of the longest ORF
        longest_orf = 0
        for i in range(len(match_stops) - 1):
            orf_len = match_stops[i + 1] - match_stops[i] - 1
            if orf_len > longest_orf:
                longest_orf = orf_len
        # Get the coverage
        coverage = (longest_orf / query_aa_seq_len) * 100
        good_candidates_rest[denovo]["longest_orf"] = longest_orf
        good_candidates_rest[denovo]["orf_qcov"] = round(coverage, 2)
        qcovs_orf_rest.append(coverage)
        # Keep only < 70
        if coverage < 70:
            good_good_candidates_rest[denovo] = good_candidates_rest[denovo]

    # Show the results
    print(f"\n\nOut of {len(good_candidates_rest)} good candidates, there are {len(good_good_candidates_rest)} for which the longest ORF covers < 70% of the de novo gene")
    show_plots_rest = input("\nDo you want to see the plots? (y/n) ")

    if show_plots_rest == "y":
        # Get the common bin edges
        bin_edges_rest = np.linspace(0, 100, 20)

        # Histograms
        plt.figure()
        plt.hist(qcovs_rest, edgecolor = "black", bins=bin_edges_rest)
        plt.xlim([0, 100])
        plt.ylim([0, 10])
        plt.xlabel("Coverage of query sequence (%)")
        plt.ylabel("Number of hits")
        plt.title(f"Coverage of de novo genes by the tblastn hits\nfor the good candidates (n = {len(good_candidates_rest)})")

        plt.figure()
        plt.hist(qcovs_orf_rest, edgecolor = "black", bins=bin_edges_rest)
        plt.xlim([0, 100])
        plt.ylim([0, 10])
        plt.vlines(x = 70, ymin = 0, ymax = 10, color = "red", linestyle = "dashed")
        plt.xlabel("Coverage of query sequence (%)")
        plt.ylabel("Number of ORFs")
        plt.title(f"Coverage of de novo genes by the longest ORF in the tblastn hits\nfor the good candidates (n = {len(good_candidates_rest)})")
        plt.show()
    
    # Display ORFs with stops
    print(f"\n\nHere are where the stops are for the {len(good_good_candidates_rest)} chosen candidates:\t\t\tLongest ORF length / de novo gene length (coverage) - origin\n")
    for denovo in good_good_candidates_rest:
        genome = denovo_dict[denovo]["genome"]
        seq = good_candidates_rest[denovo]["seq_aa"]
        stops = good_candidates_rest[denovo]["stops"]
        relative_stops = [math.floor((s/len(seq))*100) for s in stops]
        orf_lst = ["-"] * 100
        for i in relative_stops:
            orf_lst[i] = "*"
        orf_str = "".join(orf_lst)
        origin_sup0 = [[k, origin_dict[denovo][k]] for k in origin_dict[denovo] if origin_dict[denovo][k] != 0]
        origin_str = ""
        for x in origin_sup0:
            origin_str += f"{x[0]} ({str(x[1])}) + "
        print(f"{orf_str}   {good_candidates_rest[denovo]["longest_orf"]} / {len(denovo_dict[denovo]["sequence"])} ({good_candidates_rest[denovo]["orf_qcov"]}%) - {origin_str[:-3]}\n")
    


        


            


                

        










