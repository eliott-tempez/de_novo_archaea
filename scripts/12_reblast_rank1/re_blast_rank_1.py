import os
import glob
import tempfile
import subprocess
import pandas as pd
from Bio import SeqIO

GENERA_DIR = "/datas/ELIOTT/archaea_data/genera/out/"
CDS_DIR = "/datas/ELIOTT/archaea_data/genera/CDS/"
FOCAL_SPECIES = "GCA_001433455@Thermococcus_barophilus_CH5"
FOCAL_TAXID = 55802
NR = "/datas/NR/nr_2.0.13.dmnd"


def run_diamond(query_sequences):
    # Create a temporary file for the query sequences
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as query_file:
        for trg_id, trg_seq in query_sequences.items():
            query_file.write(f">{trg_id}\n{trg_seq}\n")
        query_file_path = query_file.name
    # Create a file for the BLAST output
    output_file_path = "diamond_output_re_rank_1.tsv"

    # Run the BLAST (diamond)
    print("(Running diamond...)")
    output_format = "6 qseqid sseqid staxids evalue qcovhsp"
    diamond_command = f"diamond blastp -d {NR} -p 8 -q {query_file_path} -o {output_file_path} -f {output_format} --very-sensitive --evalue 1e-3 --max-target-seqs 0 --quiet" 
    diamond_std = subprocess.run(diamond_command, shell=True)
    if diamond_std.returncode != 0:
        raise RuntimeError(f"diamond command failed with return code {diamond_std.returncode}")
    os.remove(query_file_path)
    print(f"Output written to {output_file_path}")




if __name__ == "__main__":
    # Load the genera results
    gene_age_pattern = os.path.join(GENERA_DIR, FOCAL_SPECIES, "*gene_ages.tsv")
    gene_age_file = glob.glob(gene_age_pattern)[0]
    gene_age_df = pd.read_csv(gene_age_file, sep="\t", header=0)
    # Keep only the rank 1 genes
    gene_age_df = gene_age_df[gene_age_df["rank"] == 1]
    # Sample at most 1000 genes
    gene_age_df = gene_age_df.sample(min(1000, gene_age_df.shape[0]))
    print(f"Number of rank 1 genes sampled: {gene_age_df.shape[0]}")
    # Get the gene names
    gene_names = gene_age_df["#gene"].tolist()

    # Get the gene sequences from the fasta file
    fasta_file = os.path.join(CDS_DIR, FOCAL_SPECIES + "_CDS.faa")
    gene_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.name in gene_names:
            gene_sequences[record.name] = record.seq
    print(f"Number of gene sequences retrieved: {len(gene_sequences)}")
    # Run diamond and capture output
    run_diamond(gene_sequences)