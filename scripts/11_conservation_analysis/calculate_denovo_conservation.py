"""
This script calculates the conservation of de novo genes across the different genomes.
It also prints the fasta files of all homolog CDSs for each de novo gene.
"""

import re
import subprocess
import tempfile
import os
import concurrent.futures
import pandas as pd
from Bio import SeqIO
from my_functions.genomic_functions import extract_denovo_info, get_sequence_from_loci


from my_functions.paths import GENOMES_LIST, CDS_DIR, FA_DIR
GOOD_CANDIDATES_FILE = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/18_clustering/good_candidates_clustering.tsv"
OUT_FOLDER = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/20_de_novo_conservation/"




def get_qcov(blast_results, subject_genome):
    # Get the coordinates
    strand = "+" if int(blast_results.iloc[0]["sframe"]) >= 0 else "-"
    sstart = int(blast_results.iloc[0]["sstart"]) - 1 if strand == "+" else int(blast_results.iloc[0]["send"]) - 1
    send = int(blast_results.iloc[0]["send"]) if strand == "+" else int(blast_results.iloc[0]["sstart"])
    contig = blast_results.iloc[0]["sseqid"]
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
    qcov = (longest_orf_length / blast_results.iloc[0]["qlen"]) * 100
    return qcov



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
        return False, None, None
    blast_results = pd.read_csv(output_file, sep="\t", header=None)
    blast_results.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                             "qstart", "qend", "sstart", "send", "evalue", "bitscore", 
                             "qlen", "qcovhsp", "sframe"]
    # Keep only if qcovhsp is above 50%
    blast_results = blast_results[blast_results["qcovhsp"] > 50]
    # If one row or more, we have a homology signal
    is_homology = not blast_results.empty
    # Get the best hit sseqid
    sseqid = None
    if is_homology:
        # Keep only the best result
        blast_results = blast_results.sort_values(by=["evalue", "qcovhsp"], ascending=[True, False])
        blast_results = blast_results.iloc[[0]]
        sseqid = blast_results.iloc[0]["sseqid"]
    
    # Remove the temporary output file
    os.remove(output_file)

    # If nc match, check the integrity
    qcov = None
    if blast_type == "tblastn" and is_homology:
        qcov = get_qcov(blast_results, subject_genome)

    return is_homology, qcov, sseqid



import concurrent.futures

def process_neighbour_genome(neighbour_genome, genome, seq):
    if neighbour_genome == genome:
        return neighbour_genome, None, None

    conservation_status = None

    ## Start with a blastp
    # Create a temporary file for the query sequence
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as query_file:
        query_file.write(f">query_sequence\n{seq}\n")
    # Perform a blastp against the neighbour genome's CDS
    neighbour_cds_fasta = f"{CDS_DIR}{neighbour_genome}_CDS.faa"
    is_cds, _, sseqid = blast(
        query_fasta_file=query_file.name,
        subject_fasta_file=neighbour_cds_fasta,
        blast_type="blastp",
        subject_genome=neighbour_genome
    )

    ## If no match, try with a tblastn
    if not is_cds:
        neighbour_genome_fasta = f"{FA_DIR}{neighbour_genome}.fa"
        is_nc, qcov, sseqid = blast(
            query_fasta_file=query_file.name,
            subject_fasta_file=neighbour_genome_fasta,
            blast_type="tblastn",
            subject_genome=neighbour_genome
        )
    else:
        is_nc, qcov, sseqid = False, None, sseqid

    os.remove(query_file.name)

    # Determine conservation status
    if is_cds:
        conservation_status = "cds"
    elif not is_nc:
        conservation_status = "no match"
    else:
        is_integral = qcov >= 70
        conservation_status = "integral nc" if is_integral else "fragmented nc"
    return neighbour_genome, conservation_status, sseqid



def look_for_conservation(genome, seq):
    conservation = {}
    homologs = {}
    # Process genomes parallelly
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = {
            executor.submit(process_neighbour_genome, neighbour_genome, genome, seq): neighbour_genome
            for neighbour_genome in genomes
        }
        for future in concurrent.futures.as_completed(futures):
            neighbour_genome = futures[future]
            try:
                result_genome, status, sseqid = future.result()
                if status:
                    conservation[result_genome] = status
                    if status == "cds":
                        homologs[result_genome] = sseqid

            except Exception as e:
                print(f"Error processing {neighbour_genome}: {e}")
    return conservation, homologs




if __name__ == "__main__":
    # Read list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]

    # Read list of good candidates
    clusters = pd.read_csv(GOOD_CANDIDATES_FILE, sep="\t", header=None)
    clusters.columns = ["cluster_name", "gene_name"]
    good_candidates = list(clusters["gene_name"])

    # Extract denovo sequences
    denovo_genes = {}
    for genome in genomes:
        local_dict = extract_denovo_info(genome)
        for denovo in local_dict:
            local_dict[denovo]["genome"] = genome
            denovo_genes[denovo] = local_dict[denovo]

    conservations = {}
    homologs = {}
    # For each de novo, look for conservation in all neighbours
    for denovo in good_candidates:
        genome = denovo_genes[denovo]["genome"]
        seq = denovo_genes[denovo]["sequence"]

        conservations[denovo], homologs[denovo] = look_for_conservation(
            genome=genome,
            seq=seq
        )
        
    # Print the results to a file
    with open(f"{OUT_FOLDER}de_novo_conservation.tsv", "w") as f:
        f.write("denovo\tgenome\tneighbour_genome\tstatus\n")
        for denovo in conservations:
            genome = denovo_genes[denovo]["genome"]
            for neighbour_genome in conservations[denovo]:
                status = conservations[denovo][neighbour_genome]
                f.write(f"{denovo}\t{genome}\t{neighbour_genome}\t{status}\n")

    # Print the fasta files
    os.makedirs(f"{OUT_FOLDER}coding_homologs", exist_ok=True)
    for denovo in homologs:
        genome = denovo_genes[denovo]["genome"]
        file_name = f"{OUT_FOLDER}coding_homologs/{denovo}.fasta"
        with open(file_name, "w") as f:
            # Print the denovo sequence
            f.write(f">{denovo}|{genome}|denovo\n{denovo_genes[denovo]['sequence']}\n")
            # Print the homologs
            for neighbour_genome in homologs[denovo]:
                # Get the homolog sequence
                cds_file_neighbour = f"{CDS_DIR}{neighbour_genome}_CDS.faa"
                homolog_name = homologs[denovo][neighbour_genome]
                for record in SeqIO.parse(cds_file_neighbour, "fasta"):
                    if record.name == homolog_name:
                        f.write(f">{homolog_name}|{neighbour_genome}|homolog\n{str(record.seq)}\n")
                        break