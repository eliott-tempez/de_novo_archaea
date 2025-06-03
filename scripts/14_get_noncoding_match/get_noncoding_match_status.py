"""
This script retrieves the origin of a given de novo gene
by looking at its non-coding match in a given genome.

Inputs:
* -d / --gene_id: ID of the de novo gene.
* -f / --focal_genome: focal genome name.
* -n / --nc_genome: genome to check for non-coding matches.
"""


##### Imports #####
import sys
import os
import argparse
import tempfile
import subprocess
import concurrent.futures
import numpy as np
import pandas as pd
from Bio import SeqIO

# Filenames
from my_functions.paths import DENSE_DIR, GENOMES_LIST, CDS_DIR, GFF_DIR, FA_DIR


##### Functions #####
def get_args():
    parser = argparse.ArgumentParser(
        description="Get the origin of a de novo gene based on its non-coding match in a genome."
    )
    parser.add_argument("-d", "--gene_id", type=str, required=True, help="The ID of the de novo gene to check.")
    parser.add_argument("-f", "--focal_genome", type=str, required=True, help="The name of the focal genome with the de novo gene.")
    parser.add_argument("-n", "--nc_genome", type=str, required=True, help="The genome to check against.")
    return parser.parse_args()


def get_sequence_from_loci(genome, contig, start, end):
    fa_file = os.path.join(FA_DIR, genome + ".fa")
    for record in SeqIO.parse(fa_file, "fasta"):
        if str(record.name) == str(contig):
            return record.seq[start:end]
    return None


def get_cds_info(genome):
    """Get the list of all CDSs for a given genome"""
    sp_genes = []
    gff_file = os.path.join(GFF_DIR, genome + ".gff3")
    gff_df = pd.read_csv(gff_file, sep="\t", header=None, comment="#")
    gff_df.columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    cds_df = gff_df[gff_df["type"] == "CDS"].reset_index(drop=True)
    for row in cds_df.itertuples():
        strand = row.strand
        contig = row.seqid
        start = row.start - 1
        end = row.end
        seq = ""
        cds_name = row.attributes.split(";")[0].split("=")[1]
        fa_file = os.path.join(FA_DIR, genome + ".fa")
        for record in SeqIO.parse(fa_file, "fasta"):
            if record.name == contig:
                if strand == "+":
                    seq = record.seq[start:end]
                else:
                    seq = record.seq[start:end].reverse_complement
        sp_genes.append([contig, strand, int(start), int(end), seq, cds_name])
    return sp_genes


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


def get_nc_seq_in_frame(seq, start, end, start_match, end_match, strand):
    if strand == "+":
        begin_at = frame_start(start_match, start, strand)
        end_at =  frame_end(end_match, begin_at, strand) - start_match
        begin_at -= start_match
    elif strand == "-":
        end_at = frame_start(end_match, end, strand)
        begin_at = frame_end(start_match, end_at, strand) - start_match
        end_at -= start_match
    return seq[begin_at:end_at]


def run_blast(query_file_name, subject_file_name, output_file_name):
    cmd = f'tblastn -query {query_file_name} -subject {subject_file_name} -word_size 2 -out {output_file_name} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovs sframe"'
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(e)


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
    
    # Run tblastn in parallel
    with concurrent.futures.ThreadPoolExecutor() as executor:
        future = executor.submit(run_blast, query_file_name, subject_file_name, output_file_name)
        future.result()  # Wait for the BLAST search to complete
    
    # Read the output
    blast_df = pd.read_csv(output_file_name, sep="\t", header=None)
    blast_df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "qcov", "sframe"]
    blast_df = blast_df.sort_values(by=["qcov", "evalue"], ascending=[False, True])
    frame = blast_df.iloc[0]["sframe"]
    # Remove temp files
    os.remove(query_file_name)
    os.remove(subject_file_name)
    os.remove(output_file_name)
    return int(frame)





def extract_denovo_info(gene_id, focal_genome, nc_genome):
    denovo_dict = {
        "gene_id": gene_id,
        "focal_genome": focal_genome,
        "nc_genome": nc_genome}

    # Get the de novo sequence
    fna_file = os.path.join(CDS_DIR, focal_genome + "_CDS.faa")
    for record in SeqIO.parse(fna_file, "fasta"):
        if record.name == gene_id:
            denovo_dict["sequence"] = record.seq
            break

    # Get the NC hit position
    tblastn_result_file = f"{DENSE_DIR}{focal_genome}/blast_out/TRG_multielongated_tblastn_{nc_genome}_genome.out"
    tblastn_df = pd.read_csv(tblastn_result_file, sep="\t", header=None, comment="#")
    tblastn_df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlength", "qcov", "sframe"]
    # Get the rows corresponding to the denovo gene
    tblastn_df_denov = tblastn_df[tblastn_df["qseqid"] == f"{gene_id}_elongated"].reset_index(drop=True)
    # Sort by qcov then evalue
    tblastn_df_denov = tblastn_df_denov[tblastn_df_denov["qcov"] >= 50]
    tblastn_df_denov = tblastn_df_denov.sort_values(by="evalue")

    # Extract the blast result
    strand = "+" if tblastn_df_denov.iloc[0]["sstart"] < tblastn_df_denov.iloc[0]["send"] else "-"
    contig = tblastn_df_denov.iloc[0]["sseqid"]
    start_loci = tblastn_df_denov.iloc[0]["sstart"] - 1 if strand == "+" else tblastn_df_denov.iloc[0]["send"] - 1
    end_loci = tblastn_df_denov.iloc[0]["send"] if strand == "+" else tblastn_df_denov.iloc[0]["sstart"]
    # Add the information to the dictionary
    denovo_dict["loci"] = [contig, int(start_loci), int(end_loci), strand]
    # Get the part of the de novo gene that matched
    denovo_dict["qstart"] = tblastn_df_denov.iloc[0]["qstart"] - 1
    denovo_dict["qend"] = tblastn_df_denov.iloc[0]["qend"]

    return denovo_dict


def get_nc_origin(denovo_info):
    origin_frames = {}
    nc_genome = denovo_info["nc_genome"]

    # Get the nc sequence
    nc_seq = get_sequence_from_loci(denovo_info["nc_genome"], *denovo_info["loci"][:3])
    if nc_seq is None:
        print(f"No sequence found for {denovo_info['nc_genome']} at {denovo_info['loci']}")
        return origin_frames
    
    # Get all the CDSs in the nc genome
    nc_cdss = get_cds_info(nc_genome)
    contig_match, start_match, end_match, strand_match = denovo_info["loci"]
    # Keep the CDSs that contain the de novo gene
    de_novo_is_on_plus = strand_match == "+"
    all_loci = list(range(start_match, end_match))
    loci_in_gene = []
    n_different_genes = 0
    different_genes = []
    for cds in nc_cdss:
        contig, strand, start, end, seq, cds_name = cds

        # For the + strand
        if de_novo_is_on_plus:
            if strand == "+" and contig == contig_match:
                if any(start <= i < end for i in all_loci):
                    n_different_genes += 1
                    different_genes.append(cds_name)
                    # Get the nucleotides in the CDS
                    loci_in_gene += [i for i in all_loci if start <= i < end]

                    # Get the frame
                    # Query protein sequence
                    qstart = denovo_info["qstart"]
                    qend = denovo_info["qend"]
                    query = denovo_info["sequence"][qstart:qend]
                    # Subject nucleotide sequence : in frame with the subject gene
                    subject = get_nc_seq_in_frame(nc_seq, start, end, start_match, end_match, strand_match)
                    frame = get_frame_from_blast(query, subject)
                    if f"f{frame - 1}" not in origin_frames:
                        origin_frames[f"f+{frame - 1}"] = len([i for i in all_loci if start <= i < end])
                    else:
                        origin_frames[f"f+{frame - 1}"] += len([i for i in all_loci if start <= i < end])

        # For the - strand
        else:
            if strand == "-" and contig == contig_match:
                if any(start <= i < end for i in all_loci):
                    n_different_genes += 1
                    different_genes.append(cds_name)
                    loci_in_gene += [i for i in all_loci if start <= i < end]

                    # Get the frame
                    # Query protein sequence
                    qstart = denovo_info["qstart"]
                    qend = denovo_info["qend"]
                    query = denovo_info["sequence"][qstart:qend]
                    # Subject nucleotide sequence
                    nc_seq = nc_seq.reverse_complement()
                    subject = get_nc_seq_in_frame(nc_seq, start, end, start_match, end_match, strand_match)
                    frame = get_frame_from_blast(query, subject)
                    if f"f{frame - 1}" not in origin_frames:
                        origin_frames[f"f+{frame - 1}"] = len([i for i in all_loci if start <= i < end])
                    else:
                        origin_frames[f"f+{frame - 1}"] += len([i for i in all_loci if start <= i < end])

    # Keep the number of genes we found a match in
    origin_frames["n_different_genes"] = n_different_genes
    origin_frames["different_genes"] = different_genes

    # Get the loci that are in the intergenic
    loci_in_intergenic = [i for i in all_loci if i not in loci_in_gene]
    if len(loci_in_intergenic) > 0:
        origin_frames["intergenic"] = len(loci_in_intergenic) if de_novo_is_on_plus else len(loci_in_intergenic)

    return origin_frames










if __name__ == "__main__":
    args = get_args()
    gene_id = args.gene_id
    focal_genome = args.focal_genome
    nc_genome = args.nc_genome

    # Extract the info for the de novo gene
    denovo_info = extract_denovo_info(gene_id, focal_genome, nc_genome)

    # Extract the origin frames
    origin_frames = get_nc_origin(denovo_info)

    # Print to stdout
    overlapping_genes = origin_frames["different_genes"]
    n_overlapping_genes = origin_frames["n_different_genes"]
    f0 = origin_frames.get("f+0", 0)
    f1 = origin_frames.get("f+1", 0)
    f2 = origin_frames.get("f+2", 0)
    intergenic = origin_frames.get("intergenic", 0)
    nc_match_contig, nc_match_start, nc_match_end, nc_match_strand = denovo_info["loci"]

    result_dtf = pd.DataFrame({
        "gene_id": [gene_id],
        "focal_genome": [focal_genome],
        "nc_genome": [nc_genome],
        "nc_match_contig": nc_match_contig,
        "nc_match_strand": nc_match_strand,
        "nc_match_start": nc_match_start,
        "nc_match_end": nc_match_end,
        "n_overlapping_genes": [n_overlapping_genes],
        "overlapping_genes": [";".join(overlapping_genes)],
        "f0": [f0],
        "f1": [f1],
        "f2": [f2],
        "intergenic": [intergenic]
    })
    result_dtf.to_csv(sys.stdout, sep="\t", index=False)

