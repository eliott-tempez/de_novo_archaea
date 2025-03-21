import os
import re
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import tempfile
import subprocess
import numpy as np




DENSE_DIR = "/home/eliott.tempez/Documents/archaea_data/dense/"
DATA_DIR = "/home/eliott.tempez/Documents/archaea_data/complete_122/"
GENOMES_LIST = "/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/genera_archaea/genomes_list.txt"
OUT_FOLDER = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/explore_dense_results/"



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



def get_CDS_info(genome):
    """Get the list of all CDSs for a given genome"""
    sp_genes = []
    gff_file = os.path.join(DATA_DIR, "gff3_no_fasta", genome + ".gff3")
    gff_df = pd.read_csv(gff_file, sep="\t", header=None, comment="#")
    gff_df.columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    cds_df = gff_df[gff_df["type"] == "CDS"].reset_index(drop=True)
    for row in cds_df.itertuples():
        strand = row.strand
        contig = row.seqid
        start = row.start - 1
        end = row.end
        seq = ""
        fa_file = os.path.join(DATA_DIR, "fasta_renamed", genome + ".fa")
        for record in SeqIO.parse(fa_file, "fasta"):
            if record.name == contig:
                if strand == "+":
                    seq = record.seq[start:end]
                else:
                    seq = record.seq[start:end].reverse_complement
        sp_genes.append([contig, strand, int(start), int(end), seq])
    return sp_genes



def get_sequence_from_loci(genome, contig, start, end):
    fa_file = os.path.join(DATA_DIR, "fasta_renamed", genome + ".fa")
    for record in SeqIO.parse(fa_file, "fasta"):
        if str(record.name) == str(contig):
            return record.seq[start:end]
    return None



def get_frame_from_blast(query, subject):
    """Get the frame of the query sequence in the subject sequence"""
    # Create temp files
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as query_file:
        query_file.write(f">query\n{query}")
        query_file_name = query_file.name
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as subject_file:
        subject_file.write(f">subject\n{subject}")
        subject_file_name = subject_file.name
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as output_file:
        output_file_name = output_file.name
    
    # Run tblastn
    cmd = f'tblastn -query {query_file_name} -subject {subject_file_name} -out {output_file_name} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs sframe"'
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(e)
    
    # Read the output
    blast_df = pd.read_csv(output_file_name, sep="\t", header=None)
    blast_df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "qcov", "sframe"]
    blast_df = blast_df.sort_values(by="evalue")
    frame = blast_df.iloc[0]["sframe"]
    # Remove temp files
    os.remove(query_file_name)
    os.remove(subject_file_name)
    os.remove(output_file_name)
    return int(frame)


def frame_start(x, origin, strand):
    if strand == "+":
        if np.abs(x - origin) % 3 == 0:
            return x
        elif np.abs(x - origin) % 3 == 1 and x > origin:
            return x + 2
        elif np.abs(x - origin) % 3 == 1 and x < origin:
            return x + 1
        elif np.abs(x - origin) % 3 == 2 and x > origin:
            return x + 1
        elif np.abs(x - origin) % 3 == 2 and x < origin:
            return x + 2
    elif strand == "-":
        if np.abs(x - origin) % 3 == 0:
            return x
        elif np.abs(x - origin) % 3 == 1 and x > origin:
            return x - 1
        elif np.abs(x - origin) % 3 == 1 and x < origin:
            return x - 2
        elif np.abs(x - origin) % 3 == 2 and x > origin:
            return x - 2
        elif np.abs(x - origin) % 3 == 2 and x < origin:
            return x - 1



def frame_end(end_match, frame_start, strand):
    if strand == "+":
        if (end_match - frame_start) % 3 == 0:
            return end_match
        elif (end_match - frame_start) % 3 == 1:
            return end_match - 1
        elif (end_match - frame_start) % 3 == 2:
            return end_match - 2
    elif strand == "-":
        if (frame_start - end_match) % 3 == 0:
            return end_match
        elif (frame_start - end_match) % 3 == 1:
            return end_match + 2
        elif (frame_start - end_match) % 3 == 2:
            return end_match + 1



def get_ancestor_seq_in_frame(seq, start, end, start_match, end_match, strand):
    if strand == "+":
        begin_at = frame_start(start_match, start, strand)
        end_at =  frame_end(end_match, begin_at, strand) - start_match
        begin_at -= start_match
    elif strand == "-":
        end_at = frame_start(end_match, end, strand)
        begin_at = frame_end(start_match, end_at, strand) - start_match
        end_at -= start_match
    return seq[begin_at:end_at]
        



def get_nc_origin(genome, denovo_dict):
    origin_frames = {}
    for denovo in denovo_dict:
        origin_frames[denovo] = {}
        seq_ancestor_match = get_sequence_from_loci(denovo_dict[denovo]["ancestor_sp"], denovo_dict[denovo]["loci"][0], denovo_dict[denovo]["loci"][1], denovo_dict[denovo]["loci"][2])
        if seq_ancestor_match is None:
            print(f"Could not find the sequence for {denovo} ({genome}) in {denovo_dict[denovo]['ancestor_sp']}")
            continue

        # Get all the CDss in the ancestor genome
        ancestor_cdss = get_CDS_info(denovo_dict[denovo]["ancestor_sp"])
        contig_match, start_match, end_match, strand_match = denovo_dict[denovo]["loci"]
        # Keep the CDSs that contain the de novo gene
        de_novo_is_on_plus = strand_match == "+"
        all_loci = list(range(start_match, end_match+1))
        loci_in_gene = []
        for cds in ancestor_cdss:
            contig, strand, start, end, seq = cds

            # For the + strand
            if de_novo_is_on_plus:
                if strand == "+" and contig == contig_match:
                    if any(start <= i <= end for i in all_loci):
                        # Get the nucleotides in the CDS
                        loci_in_gene += [i for i in all_loci if start <= i <= end]

                        # Get the frame
                        # Query protein sequence
                        qstart = denovo_dict[denovo]["qstart"]
                        qend = denovo_dict[denovo]["qend"]
                        query = denovo_dict[denovo]["sequence"][qstart:qend]
                        # Subject nucleotide sequence : in frame with the subject gene
                        subject = get_ancestor_seq_in_frame(seq_ancestor_match, start, end, start_match, end_match, strand_match)
                        frame = get_frame_from_blast(query, subject)
                        if f"f{frame - 1}" not in origin_frames[denovo]:
                            origin_frames[denovo][f"f+{frame - 1}"] = len([i for i in all_loci if start <= i <= end])
                        else:
                            origin_frames[denovo][f"f+{frame - 1}"] += len([i for i in all_loci if start <= i <= end])

            # For the - strand
            else:
                if strand == "-" and contig == contig_match:
                    if any(start <= i <= end for i in all_loci):
                        loci_in_gene += [i for i in all_loci if start <= i <= end]
                        # Get the frame
                        # Query protein sequence
                        qstart = denovo_dict[denovo]["qstart"]
                        qend = denovo_dict[denovo]["qend"]
                        query = denovo_dict[denovo]["sequence"][qstart:qend]
                        # Subject nucleotide sequence
                        seq_ancestor_match = seq_ancestor_match.reverse_complement()
                        subject = get_ancestor_seq_in_frame(seq_ancestor_match, start, end, start_match, end_match, strand_match)
                        frame = get_frame_from_blast(query, subject)
                        if f"f{frame - 1}" not in origin_frames[denovo]:
                            origin_frames[denovo][f"f+{frame - 1}"] = len([i for i in all_loci if start <= i <= end])
                        else:
                            origin_frames[denovo][f"f+{frame - 1}"] += len([i for i in all_loci if start <= i <= end])
        # Get the loci that are in the intergenic
        loci_in_intergenic = [i for i in all_loci if i not in loci_in_gene]
        if len(loci_in_intergenic) > 0:
            origin_frames[denovo]["intergenic"] = len(loci_in_intergenic) if de_novo_is_on_plus else len(loci_in_intergenic)

    return origin_frames



if __name__ == "__main__":
    # Read list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]

    # Get the de novo genes info
    denovo_dict = {}
    origin_frames = {}
    results_list = []
    n_denovo = 0

    for genome in genomes:
        denovo_dict[genome] = get_denovo_info(genome)
        # Keep only the genomes with de novo genes
        if denovo_dict[genome] == {}:
            continue
            

        # Get the corresponding area in the ancestor genome
        origin_frames[genome] = get_nc_origin(genome, denovo_dict[genome])
    
    
        # Store in pandas dataframe
        for denovo in origin_frames[genome]:
            outgroup = denovo_dict[genome][denovo]["ancestor_sp"]
            noncoding_match_contig, noncoding_match_start, noncoding_match_end, noncoding_match_strand = denovo_dict[genome][denovo]["loci"]
            intergenic = origin_frames[genome][denovo]["intergenic"] if "intergenic" in origin_frames[genome][denovo] else 0
            f0 = origin_frames[genome][denovo]["f+0"] if "f+0" in origin_frames[genome][denovo] else 0
            f1 = origin_frames[genome][denovo]["f+1"] if "f+1" in origin_frames[genome][denovo] else 0
            f2 = origin_frames[genome][denovo]["f+2"] if "f+2" in origin_frames[genome][denovo] else 0
            results_list.append({"genome": genome, "denovo_gene": denovo, "outgroup": outgroup, "noncoding_match_contig": noncoding_match_contig, "noncoding_match_start": noncoding_match_start, "noncoding_match_end": noncoding_match_end, "noncoding_match_strand": noncoding_match_strand, "intergenic": intergenic, "f+0": f0, "f+1": f1, "f+2": f2})
            n_denovo += 1
        
    print(f"{n_denovo} de novo genes found in total")
    results = pd.DataFrame(results_list, columns=["genome", "denovo_gene", "outgroup", "noncoding_match_contig", "noncoding_match_start", "noncoding_match_end", "noncoding_match_strand", "intergenic", "f+0", "f+1", "f+2"])
    # Save to file
    results.to_csv(os.path.join(OUT_FOLDER, "denovo_noncoding_status.tsv"), sep="\t", index=False)

        


    







