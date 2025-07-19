"""
This script calculates the conservation of 100 random ORFs for each genome,
compared with all 115 other genomes.
"""

import re
import subprocess
import tempfile
import os
import sys
import random
import concurrent.futures
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from my_functions.genomic_functions import get_sequence_from_loci
from my_functions.paths import GENOMES_LIST, CDS_DIR, FA_DIR, GFF_DIR
OUT_FOLDER = "../results/11_conservation_analysis/"


def extract_100_orfs(genome):
    fa_file = os.path.join(FA_DIR, genome + ".fa")
    gff_file = os.path.join(GFF_DIR, genome + ".gff3")
    orftrack_file = f"mapping_orf_{genome}.gff"
    orfget_file = f"mapping_orf_{genome}_nc_intergenic.pfasta"
    log_file = "orftrack.log"

    # Use orftrack to annotate the genome
    result_orftrack = subprocess.run(
        ["orftrack", "--fna", fa_file, "--gff", gff_file],
        capture_output=True,
        text=True
    )
    if result_orftrack.returncode != 0:
        print(f"Error in orftrack: {result_orftrack.stderr}")
    
    # Use orfget to extract ORFs
    result_orfget = subprocess.run(
        ["orfget", "--fna", fa_file, "--gff", orftrack_file, "--features_include", "nc_intergenic"],
        capture_output=True,
        text=True
    )
    if result_orfget.returncode != 0:
        print(f"Error in orfget: {result_orfget.stderr}")
    
    # Read the ORFs
    all_orfs = []
    with open(orfget_file, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq = line.strip()
                all_orfs.append(seq)

    # Sample 100 random ORFs and create dict
    orfs = {}
    sampled_orfs = random.sample(all_orfs, 100)
    for i, seq in enumerate(sampled_orfs):
        orf_id = f"{genome}_orf_{i+1}"
        orfs[orf_id] = {
            "sequence": seq,
            "genome": genome
        }

    # Delete temporary files
    os.remove(orftrack_file)
    os.remove(orfget_file)
    os.remove(log_file)
    return orfs
                    


def get_qcov(blast_results, subject_genome):
    n_integral_nc = 0
    n_fragmented_nc = 0

    # For each row
    for index, row in blast_results.iterrows():
        # Get the coordinates
        strand = "+" if int(row["sframe"]) >= 0 else "-"
        sstart = int(row["sstart"]) - 1 if strand == "+" else int(row["send"]) - 1
        send = int(row["send"]) if strand == "+" else int(row["sstart"])
        contig = row["sseqid"]
        # Get the sequence 
        subject_seq = get_sequence_from_loci(
            genome=subject_genome,
            contig=contig,
            start=sstart,
            end=send,
            strand=strand
        )
        # Extract the length of the longest orf
        orfs = subject_seq.translate().split("*")
        longest_orf_length = max(len(orf) for orf in orfs)
        # Calculate the qcov
        qcov = (longest_orf_length / row["qlen"]) * 100

        if qcov >= 70:
            n_integral_nc += 1
        else:
            n_fragmented_nc += 1
    return n_integral_nc, n_fragmented_nc
        



def blast(query_fasta_file, subject_fasta_file, blast_type, subject_genome):
    # Create temporary output file
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.tsv') as temp_output:
        output_file = temp_output.name
    # Run the blast command
    out_command = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovhsp sframe"
    command = [
        blast_type,
        "-query", query_fasta_file,
        "-subject", subject_fasta_file,
        "-outfmt", out_command,
        "-out", output_file,
        "-evalue", "1e-3"
    ]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Read the output
    # Check the file isn't blank
    if os.path.getsize(output_file) == 0:
        return 0, []
    blast_results = pd.read_csv(output_file, sep="\t", header=None)
    os.remove(output_file)
    blast_results.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                             "qstart", "qend", "sstart", "send", "evalue", "bitscore", 
                             "qlen", "qcovhsp", "sframe"]
    # Keep only if qcovhsp is above 50%
    blast_results = blast_results[blast_results["qcovhsp"] > 50]
    # Keep only one row per query sequence
    blast_results = blast_results.sort_values(by=["evalue", "qcovhsp"], ascending=[True, False])
    blast_results = blast_results.drop_duplicates(subset=["qseqid"], keep='first')

    if blast_results.empty:
        return (0, []) if blast_type == "blastp" else (0, 0)

    # If non-coding, check the integrity
    if blast_type == "tblastn":
        n_integral_nc, n_fragmented_nc = get_qcov(blast_results, subject_genome)
        return n_integral_nc, n_fragmented_nc
    
    # If coding, return the names of the matched orfs as well
    matched_orfs = list(blast_results["qseqid"])
    return len(blast_results), matched_orfs







def process_neighbour_genome(orfs, genome, neighbour_genome):
    """Process a single neighbour genome and return its conservation status."""
    neighbour_cds_fasta = f"{CDS_DIR}{neighbour_genome}_CDS.faa"
    neighbour_genome_fasta = f"{FA_DIR}{neighbour_genome}.fa"

    # Create a temporary fasta file with the ORFs
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as query_file:
        for orf_id in orfs:
            seq = orfs[orf_id]["sequence"]
            query_file.write(f">{orf_id}\n{seq}\n")
    query_fasta_file = query_file.name

    # Start with blastp to get cds matches
    n_cds, orfs_to_remove = blast(
        query_fasta_file=query_fasta_file,
        subject_fasta_file=neighbour_cds_fasta,
        blast_type="blastp",
        subject_genome=neighbour_genome
    )

    # Remove matched orfs
    os.remove(query_fasta_file)
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as query_file:
        for orf_id in orfs:
            if orf_id not in orfs_to_remove:
                seq = orfs[orf_id]["sequence"]
                query_file.write(f">{orf_id}\n{seq}\n")
    query_fasta_file = query_file.name

    # Do tblastn to get nc matches
    n_integral_nc, n_fragmented_nc = blast(
        query_fasta_file=query_fasta_file,
        subject_fasta_file=neighbour_genome_fasta,
        blast_type="tblastn",
        subject_genome=neighbour_genome
    )

    n_no_match = 100 - (n_cds + n_integral_nc + n_fragmented_nc)

    # Delete temporary files
    os.remove(query_fasta_file)

    n = (n_cds, n_integral_nc, n_fragmented_nc, n_no_match)
    return neighbour_genome, n




def blast_orfs_and_get_conserv(orfs, genome, max_workers=None):
    """
    Parallelized version that processes neighbour genomes simultaneously
    """
    conservations = {}
    
    # Get list of neighbour genomes (exclude current genome)
    neighbour_genomes = [ng for ng in genomes if ng != genome]
    
    # Use ThreadPoolExecutor to process neighbour genomes in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_genome = {
            executor.submit(process_neighbour_genome, orfs, genome, neighbour_genome): neighbour_genome
            for neighbour_genome in neighbour_genomes
        }
        
        # Collect results as they complete
        for future in concurrent.futures.as_completed(future_to_genome):
            try:
                neighbour_genome, n_status = future.result()
                conservations[neighbour_genome] = n_status
                print(f"  Completed: {neighbour_genome} -> {n_status}")
            except Exception as exc:
                neighbour_genome = future_to_genome[future]
                print(f"  Error processing {neighbour_genome}: {exc}")
                conservations[neighbour_genome] = "error"

    return conservations
    




if __name__ == "__main__":
    # Read list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]
    
    random_orfs = {}
    conservations = {}

    # For each genome
    for genome in genomes:
        print(f"Processing genome: {genome}...")
        # Extract 100 random ORFs
        random_orfs = extract_100_orfs(genome)
        
        # Blast the sequences and extract the results
        conservations[genome] = blast_orfs_and_get_conserv(
            orfs=random_orfs, 
            genome=genome,
            max_workers=None  # Uses all available CPU cores
        )
        print("\n")

      
    # Print the results to a file
    with open(f"{OUT_FOLDER}orf_conservation.tsv", "w") as f:
        f.write("genome\tneighbour_genome\tn_cds\tn_integral_nc\tn_fragmented_nc\tn_no_match\n")
        for genome in conservations:
            for neighbour_genome in conservations[genome]:
                n_cds, n_integral_nc, n_fragmented_nc, n_no_match = conservations[genome][neighbour_genome]
                f.write(f"{genome}\t{neighbour_genome}\t{n_cds}\t{n_integral_nc}\t{n_fragmented_nc}\t{n_no_match}\n")