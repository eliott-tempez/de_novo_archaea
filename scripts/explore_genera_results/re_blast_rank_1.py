import os
import glob
import tempfile
import subprocess
import pandas as pd
from Bio import SeqIO

GENERA_DIR = "/home/eliott.tempez/Documents/archaea_data/genera/out/"
DATA_DIR = "/home/eliott.tempez/Documents/archaea_data/complete_122/"
FOCAL_SPECIES = "GCA_001433455@Thermococcus_barophilus_CH5"


def run_blast(query_sequences, db_faa_file):
    # Create a temporary file for the query sequences
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as query_file:
        for trg_id, trg_seq in query_sequences.items():
            query_file.write(f">{trg_id}\n{trg_seq}\n")
        query_file_path = query_file.name
    # Create a temporary file for the BLAST output
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as output_file:
        output_file_path = output_file.name

    # Run the BLAST
    result = subprocess.run(["blastp", "-query", query_file_path, "-subject", db_faa_file, "-out", output_file_path, "-outfmt", "6 qseqid sseqid qlen evalue qcovs", "-num_threads", "12"], capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(f"BLAST command failed with return code {result.returncode}: {result.stderr.decode()}")
    # Read the output file
    result = pd.read_csv(output_file_path, sep="\t", header=None)
    result.columns = ["qseqid", "sseqid", "qlen", "evalue", "qcov"]
    # Keep only rows for which qcov > 50 % and evalue < 1e-3
    result = result[(result["qcov"] > 50) & (result["evalue"] < 1e-3)]
    # Keep only the best hit for each query
    result = result.sort_values("evalue").drop_duplicates("qseqid")
    # Remove the temporary files
    os.remove(query_file_path)
    os.remove(output_file_path)
    # Return number of matches
    return len(result)


if __name__ == "__main__":
    # Load the genera results
    gene_age_pattern = os.path.join(GENERA_DIR, FOCAL_SPECIES, "*gene_ages.tsv")
    gene_age_file = glob.glob(gene_age_pattern)[0]
    gene_age_df = pd.read_csv(gene_age_file, sep="\t", header=0)
    # Keep only the rank 1 genes
    gene_age_df = gene_age_df[gene_age_df["rank"] == 1]
    # Sample at most 1000 genes
    gene_age_df = gene_age_df.sample(min(1000, gene_age_df.shape[0]))
    # Get the gene names
    gene_names = gene_age_df["#gene"].tolist()

    # Get the gene sequences from the fasta file
    fasta_file = os.path.join(DATA_DIR, "CDS", FOCAL_SPECIES + "_CDS.faa")
    gene_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.name in gene_names:
            gene_sequences[record.name] = record.seq