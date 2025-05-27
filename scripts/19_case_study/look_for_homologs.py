import os
import tempfile
import subprocess
import pandas as pd
from Bio import SeqIO


from my_functions.paths import DENSE_DIR
genome = "GCA_001484685@Thermococcus_sp_2319x1"
good_candidates_file = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/good_candidates.txt"
dense_dir = os.path.join(DENSE_DIR, genome)



def build_global_proteome(genome):
    # Create temp file
    temp_fasta = tempfile.NamedTemporaryFile(delete=False, suffix=".faa")
    temp_fasta_file = temp_fasta.name
    cds_dir = os.path.join(dense_dir, "CDS")
    # Get list of all .faa files in dir
    faa_files = [f for f in os.listdir(cds_dir) if f.endswith(".faa") and genome not in f]
    for faa_file in faa_files:
        genome_name = faa_file.split("_CDS.faa")[0]
        fasta_path = os.path.join(cds_dir, faa_file)
        # Write the sequences to the temp file
        with open(temp_fasta_file, "a") as temp_fasta:
            for record in SeqIO.parse(fasta_path, "fasta"):
                record.id = f"{genome_name}_{record.name}"
                SeqIO.write(record, temp_fasta, "fasta")
    return temp_fasta_file



def blast_denovo(denovo_gene, genome, temp_fasta):
    # Create temp file that is the copy of the temp_fasta file
    temp_fasta_denovo = tempfile.NamedTemporaryFile(delete=False, suffix=".faa")
    temp_fasta_denovo_file = temp_fasta_denovo.name
    with open(temp_fasta, "r") as tfasta, open(temp_fasta_denovo_file, "w") as tdenovo:
        for line in tfasta:
            tdenovo.write(line)

    # Open the focal fasta file
    fasta_file = os.path.join(dense_dir, "CDS", f"{genome}_CDS.faa")
    with open(temp_fasta_denovo_file, "a") as tfasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # If the record is the denovo gene, save it and skip writing it to the temp file
            if record.name == denovo_gene:
                denovo_seq = record.seq
            # Otherwise, write the record to the temp file
            else:
                record.id = f"{genome}_{record.name}"
                SeqIO.write(record, tfasta, "fasta")
    
    # Create tempfile for query and output
    query_file = tempfile.NamedTemporaryFile(delete=False, suffix=".faa")
    query_file_name = query_file.name
    with open(query_file_name, "w") as qfile:
        qfile.write(f">{denovo_gene}\n{denovo_seq}\n")
    blast_output = tempfile.NamedTemporaryFile(delete=False, suffix=".tsv")
    blast_output_file = blast_output.name
    
    # Run blastp
    blast_cmd = [
        "blastp",
        "-query", query_file_name,
        "-subject", temp_fasta_denovo_file,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp",
        "-out", blast_output_file,
        "-evalue", "1e-3"
    ]
    subprocess.run(blast_cmd, check=True)
    # Read the blast output
    # Check the output isn't empty
    if os.stat(blast_output_file).st_size == 0:
        os.remove(temp_fasta_denovo_file)
        os.remove(query_file_name)
        os.remove(blast_output_file)
        return []
    blast_df = pd.read_csv(blast_output_file, sep="\t", header=None)
    blast_df.columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp"
    ]

    # Remove temp files
    os.remove(temp_fasta_denovo_file)
    os.remove(query_file_name)
    os.remove(blast_output_file)
    # Return the blast dataframe
    return blast_df



if __name__ == "__main__":
    # Read the denovo genes
    denovo_file = os.path.join(dense_dir, "denovogenes.tsv")
    denovo_df = pd.read_csv(denovo_file, sep="\t", header=0)
    denovo_genes = denovo_df.loc[:, "CDS"].tolist()
    # Read the good candidates
    with open(good_candidates_file, "r") as f:
        good_candidates = [line.strip() for line in f.readlines()]
    good_candidates = list(set(denovo_genes) & set(good_candidates))

    # Get temp fasta file of all proteomes except focal
    temp_fasta = build_global_proteome(genome)

    # For each denovo gene, blast against the global proteome
    print(f"This is the result of a blastp against all CDSs for the denovo genes of {genome}.")
    print("#" * 60)
    print()
    for dnv in denovo_genes:
        good_bad = "good" if dnv in good_candidates else "bad"
        match_df = blast_denovo(dnv, genome, temp_fasta)
        print(f"Denovo gene: {dnv} ({good_bad}), matches found: {len(match_df)}\n")
        if len(match_df) > 0:
            print(match_df.to_string(index=False))
        print(f"\n\n{"_" * 20}\n\n")
    
    # Delete the temp fasta file
    os.remove(temp_fasta)
