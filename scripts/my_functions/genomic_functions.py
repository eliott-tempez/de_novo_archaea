"""
This script contains functions that are used for genomic analysis.
"""

import os
import subprocess
import tempfile
# Paths import
from my_functions.paths import *
# Data handling
import pandas as pd
import numpy as np
from Bio import SeqIO



def extract_intergenic_segments(species):
    """"Extract all intergenic segments from a given genome"""
    # Make sure the gff file exists and read it
    gff_file = os.path.join(GFF_DIR, species + ".gff3")
    if not os.path.exists(gff_file):
        raise FileNotFoundError(f"{gff_file} does not exist")
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
            if row["type"] == "gene":
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
    fa_file = os.path.join(FA_DIR, species + ".fa")
    if not os.path.exists(fa_file):
        raise FileNotFoundError(f"{fa_file} does not exist")
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


def get_end_frame(start_frame):
    if start_frame == 0:
        return 0
    elif start_frame == 1:
        return -2
    elif start_frame == 2:
        return -1


def extract_iorfs(genome, threshold=60):
    """Extract all iORFS longer than a certain length from a given genome"""
    iorfs = []
    stops = ["TAA", "TAG", "TGA"]
    intergenic_dict = extract_intergenic_segments(genome)
    # Get the sequence in all 6 frames
    for segment in intergenic_dict:
        seq = intergenic_dict[segment]["seq"]
        # Skip small segments
        if len(seq) < threshold:
            continue
        for frame in range(3):
            seq_fr = (str(seq[frame:get_end_frame(frame)])).upper()
            seq_reverse_fr = (str(seq.reverse_complement()[frame:get_end_frame(frame)])).upper()
            # Find all segments between 2 stop codons
            for s in (seq_fr, seq_reverse_fr):
                # Find first stop
                start = -1
                for i in range(0, len(s) - 2, 3):
                    codon = s[i:i + 3]
                    if codon in stops:
                        start = i
                        break
                # If no stop codon was found, skip the segment
                if start == -1:
                    continue
                # Get each stop-separated segment
                for i in range(start + 3, len(s) - 2, 3):
                    codon = s[i:i + 3]
                    if codon in stops:
                        # If the segment is longer than the threshold, add it to the list
                        if ((i + 3) - start) >= threshold:
                            iorfs.append(s[start:i + 3])
                            start = i + 3
    return iorfs
                            
            
                
                



def extract_denovo_info(genome):
    """Extract info for all de novo genes for one genome"""
    denovo_dict = {}

    # Get the name of the de novo genes
    denovo_file = os.path.join(DENSE_DIR, genome, "denovogenes.tsv")
    if not os.path.exists(denovo_file):
        raise FileNotFoundError(f"{denovo_file} does not exist")
    denovo_df = pd.read_csv(denovo_file, sep="\t", header=0)
    if denovo_df.empty:
        return {}
    for row in denovo_df.iterrows():
        denovo_dict[row[1]["CDS"]] = {}
    
    # Get the de novo sequence
    fna_file = os.path.join(CDS_DIR, genome + "_CDS.faa")
    if not os.path.exists(fna_file):
        raise FileNotFoundError(f"{fna_file} does not exist")
    for denovo_gene in denovo_dict:
        for record in SeqIO.parse(fna_file, "fasta"):
            if record.name == denovo_gene:
                denovo_dict[denovo_gene]["sequence"] = record.seq
                pass
    
    """# Get the strand of the de novo gene
    gff_file = os.path.join(GFF_DIR, genome + ".gff3")
    if not os.path.exists(gff_file):
        raise FileNotFoundError(f"{gff_file} does not exist")
    gff_pd = pd.read_csv(gff_file, sep="\t", header=None, comment="#")
    gff_pd.columns = ["contig", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    gene_names = gff_pd["attributes"].str.split(";", expand=True)[0].str.split("=", expand=True)
    gff_pd["gene_name"] = gene_names[1]
    for denovo_gene in denovo_dict:
        if "gene_mRNA" in denovo_gene:
            denovo_name = denovo_gene.split("_gene_mRNA")[0]
        # Get the strand of the de novo gene
        try:
            matches = gff_pd[gff_pd["gene_name"] == denovo_name]
            if matches.empty:
                raise ValueError(f"No corresponding row found for {denovo_gene} in {genome}")
        except ValueError as e:
            print(e)
            continue
        # Get the strand and contig
        denovo_dict[denovo_gene]["strand"] = matches.iloc[0]["strand"]"""

    # Get the last NC hit in synteny
    ## Name of the ancestor
    trg_file = os.path.join(DENSE_DIR, genome, "TRG_match_matrix.tsv")
    if not os.path.exists(trg_file):
        raise FileNotFoundError(f"{trg_file} does not exist")
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
        tblastn_result_file = os.path.join(DENSE_DIR, genome, "blast_out", f"TRG_multielongated_tblastn_{ancestor}_genome.out")
        if not os.path.exists(tblastn_result_file):
            raise FileNotFoundError(f"{tblastn_result_file} does not exist")
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



def extract_cds_info(genome):
    """Get the list of all CDSs for a given genome"""
    cdss = {}
    gff_file = os.path.join(GFF_DIR, genome + ".gff3")
    if not os.path.exists(gff_file):
        raise FileNotFoundError(f"{gff_file} does not exist")
    gff_df = pd.read_csv(gff_file, sep="\t", header=None, comment="#")
    gff_df.columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    cds_df = gff_df[gff_df["type"] == "CDS"].reset_index(drop=True)
    for row in cds_df.itertuples():
        cds_name = row.attributes.split(";")[0].split("=")[1]
        strand = row.strand
        contig = row.seqid
        start = row.start - 1
        end = row.end
        seq = ""
        fa_file = os.path.join(FA_DIR, genome + ".fa")
        if not os.path.exists(fa_file):
            raise FileNotFoundError(f"{fa_file} does not exist")
        for record in SeqIO.parse(fa_file, "fasta"):
            if str(record.name) == str(contig):
                if strand == "+":
                    seq = record.seq[start:end]
                else:
                    seq = record.seq[start:end].reverse_complement()
        cdss[cds_name] = {"contig": contig, "strand": strand, "start": start, "end": end, "seq": seq}

    return cdss



def get_sequence_from_loci(genome, contig, start, end, strand):
    """Extract the nucleotide sequence from a given loci"""
    fa_file = os.path.join(FA_DIR, genome + ".fa")
    if not os.path.exists(fa_file):
        raise FileNotFoundError(f"{fa_file} does not exist")
    for record in SeqIO.parse(fa_file, "fasta"):
        if str(record.name) == str(contig):
            if strand == "+":
                return record.seq[start:end]
            else:
                return record.seq[start:end].reverse_complement()
    return None