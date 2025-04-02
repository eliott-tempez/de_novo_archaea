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



from my_functions.paths import *
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

    denovo_dict = {}
    n_denovo = 0
    good_denovo_stops = []
    no_stops = 0
    only_final_stop = 0
    stops_at_end = 0
    good_candidates = {}

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

            # Number of sequences without stop codons
            if len(stops) == 0:
                no_stops += 1
            elif stops == [len(match_seq_aa) - 1]:
                only_final_stop += 1
            elif min(stops) >= len(match_seq_aa) - 5:
                stops_at_end += 1
            else:
                good_candidates[denovo] = {"seq": match_seq_aa, "stops": stops, "genome": genome}

    print(f"\nOut of {n_denovo} de novo genes:")
    print(f"\t{no_stops} have no stop codon")
    print(f"\t{only_final_stop} have only one stop codon at the end")
    print(f"\t{stops_at_end} have stop codon(s) only in the last 5 codons")
    print(f"\n-> {len(good_candidates)} are good candidates")

    # Display ORFs with stops
    print("\nHere are where the stops are:\n")
    for denovo in good_candidates:
        genome = good_candidates[denovo]["genome"]
        seq = good_candidates[denovo]["seq"]
        stops = good_candidates[denovo]["stops"]
        """if denovo != "KOPGNMII_00670_gene_mRNA":
            continue
        print(stops)"""
        relative_stops = [math.floor((s/len(seq))*100) for s in stops]
        orf_lst = ["-"] * 100
        for i in relative_stops:
            orf_lst[i] = "*"
        orf_str = "".join(orf_lst)
        print(f"{orf_str}\n")
    

        """# Re-blast to make sure we are in the right frame
        faa_file = os.path.join(CDS_DIR, genome + "_CDS.faa")
        for record in SeqIO.parse(faa_file, "fasta"):
            if record.name == denovo:
                db_seq = record.seq
        db = tempfile.NamedTemporaryFile(delete=False)
        db.write(f">{outgroup}\n{db_seq}\n".encode())
        db.close()
        query = tempfile.NamedTemporaryFile(delete=False)
        query.write(f">{denovo}\n{seq}\n".encode())
        query.close()
        out = tempfile.NamedTemporaryFile(delete=False)
        subprocess.run(["blastp", "-query", query.name, "-subject", db.name, "-out", out.name, "-outfmt", "6"], 
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        result = pd.read_csv(out.name, sep="\t", header=None)
        if result.empty:
            print(f"No blast result for {denovo} in {outgroup}")
            break
        result.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        # Print the lowest evalue
        print(result.sort_values(by="evalue")[["evalue", "pident"]])
        os.remove(db.name)
        os.remove(query.name)
        os.remove(out.name)"""
        


            


                

        










