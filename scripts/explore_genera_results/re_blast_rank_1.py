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


def parse_diamond_output(output_file):
    result = pd.read_csv(output_file, sep="\t", header=None)
    result.columns = ["qseqid", "sseqid", "taxids", "evalue", "qcov"]
    # Keep only rows for which qcov > 50 % and evalue < 1e-3
    result = result[(result["qcov"] > 50) & (result["evalue"] < 1e-3)]
    unique_genes = result["qseqid"].nunique()
    print(f"Out of these, {unique_genes} genes have found 1+ match")
    # for each query gene, get the highest node in the taxonomy
    lca_dict = {}
    for gene, group in result.groupby("qseqid"):
        taxids = group["taxids"].tolist()
        # Get the highest node in the taxonomy
        taxid_str = " ".join(set([str(FOCAL_TAXID)] + taxids))
        get_lca_cmd = f"echo {taxid_str} | taxonkit lca"
        lca_std = subprocess.run(get_lca_cmd, shell=True, capture_output=True)
        if lca_std.returncode != 0:
            raise RuntimeError(f"taxonkit lca command failed with return code {lca_std.returncode}")
        lca = int(lca_std.stdout.decode("utf-8").strip().split("\t")[1])
        lca_dict[gene] = lca
    # Print result to stdout
    lcas = pd.Series(lca_dict)
    print("Number of genes per LCA:")
    print(lcas.value_counts())


def run_diamond(query_sequences):
    # Create a temporary file for the query sequences
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as query_file:
        for trg_id, trg_seq in query_sequences.items():
            query_file.write(f">{trg_id}\n{trg_seq}\n")
        query_file_path = query_file.name
    # Create a temporary file for the BLAST output
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as output_file:
        output_file_path = output_file.name

    # Run the BLAST (diamond)
    print("(Running diamond...)")
    output_format = "6 qseqid sseqid staxids evalue qcovs"
    diamond_command = f"diamond blastp -d {NR} -p 8 -q {query_file_path} -o {output_file_path} -f {output_format} --very-sensitive"
    diamond_std = subprocess.run(diamond_command, shell=True)
    if diamond_std.returncode != 0:
        raise RuntimeError(f"diamond command failed with return code {diamond_std.returncode}")
    # Parse output
    parse_diamond_output(output_file_path)
    # Remove temporary files
    os.remove(query_file_path)
    os.remove(output_file_path)





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