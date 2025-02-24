import os
import glob
import re
import random
import subprocess
import tempfile
import pandas as pd
import Bio.SeqIO as SeqIO


OUTPUT_DIR = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/calculate_conservation/"
GENOMES_LIST = "/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/genera_archaea/genomes_list.txt"
GENERA_DIR = "/home/eliott.tempez/Documents/archaea_data/genera/out/"
DATA_DIR = "/home/eliott.tempez/Documents/archaea_data/complete_122/"
FOCAL_SPECIES = "GCA_001433455@Thermococcus_barophilus_CH5"
TRG_RANK = 7.0



############################## FUNCTIONS ##############################
def initialise_conservation_df(genomes, genera_dir):
    """Initialise the dataframe with empty rows and columns"""
    #col = ["n_trg", "n_cds", "n_noncoding", "n_f1", "n_f2"]
    col = ["n_trg"]
    ncol = len(col)
    conservation_df = pd.DataFrame(columns=col, index=genomes)
    genera_species = os.listdir(genera_dir)
    for g in genomes:
        # Initialise with 0 for the species for which we have the data
        if g in genera_species:
            conservation_df.loc[g] = [0]*ncol
        # Initialise with None for the species for which we don't have the data
        else:
            conservation_df.loc[g] = [None]*ncol
    return conservation_df



def extract_focal_TRGs(focal_species, trg_threshold, genera_dir):
    """Extract 1000 TRGs of the focal species"""
    # Make sure the file exists and read it
    trg_pattern = os.path.join(genera_dir, focal_species, "*gene_ages.tsv")
    trg_files = glob.glob(trg_pattern)
    if trg_files:
        trg_file = trg_files[0]
        trg_df = pd.read_csv(trg_file, sep="\t", header=0)
    else:
        raise FileNotFoundError(f"No file matching pattern {trg_pattern}")
    # Extract the TRGs
    focal_TRGs = trg_df[trg_df["rank"] >= trg_threshold]["#gene"].tolist()
    # Keep 1000 at most
    if len(focal_TRGs) > 1000:
        focal_TRGs = random.sample(focal_TRGs, 1000)
    return focal_TRGs



def run_blast(query_sequences, db_faa_file, blast_type):
    """Run BLAST for the query sequences against the database and return the number of matches"""
    # Create a temporary file for the query sequences
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as query_file:
        for trg_id, trg_seq in query_sequences.items():
            query_file.write(f">{trg_id}\n{trg_seq}\n")
        query_file_path = query_file.name
    # Create a temporary file for the BLAST output
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as output_file:
        output_file_path = output_file.name

    # Run the BLAST
    result = subprocess.run([blast_type, "-query", query_file_path, "-subject", db_faa_file, "-out", output_file_path, "-outfmt", "6 qseqid sseqid qlen evalue qcovs"], capture_output=True)
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

    return len(result)



def calculate_TRG_conservation(focal_sp, conservation_df, focal_TRGs, data_dir):
    """Calculate the conservation of the TRGs for all the species"""
    # Get the protein sequences for the TRGs
    TRG_sequences = {}
    CDS_fasta_file = os.path.join(data_dir, "CDS/" + focal_sp + "_CDS.faa")
    CDS_fasta = SeqIO.parse(CDS_fasta_file, "fasta")
    for record in CDS_fasta:
        if record.id in focal_TRGs:
            TRG_sequences[record.id] = record.seq
    if len(TRG_sequences) != len(focal_TRGs):
        print("Some TRGs are missing in the faa file")

    # For each species
    for species in conservation_df.index:
        # Skip species for which we don't have the data
        if conservation_df.loc[species, "n_trg"] is None:
            continue
        # For the focal species, add the number of TRGs
        if species == focal_sp:
            conservation_df.loc[species, "n_trg"] = len(focal_TRGs)
            continue

        # For all the other species, run BLAST for all the TRGs
        db_fa_file = os.path.join(data_dir, "fasta_renamed/" + species + ".fa")
        nb_matches = run_blast(TRG_sequences, db_fa_file, "tblastn")
        print(f"Species {species}: {nb_matches} matches")
        conservation_df.loc[species, "n_trg"] = nb_matches
    return conservation_df













        















############################## MAIN ##############################
if __name__ == "__main__":
    # Read the list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]
    # Initiate the dataframe
    conservation_df = initialise_conservation_df(genomes, GENERA_DIR)


    ########################################
    ################# TRGs #################
    ########################################
    print("Calculating TRG conservation...")
    # Extract the name of 1000 TRGs of the focal species
    focal_TRGs = extract_focal_TRGs(FOCAL_SPECIES, TRG_RANK, GENERA_DIR)
    print(f"Extracted {len(focal_TRGs)} TRGs for the focal species {FOCAL_SPECIES}\n")
    # Calculate the TRG conservation and add to dataframe
    conservation_df = calculate_TRG_conservation(FOCAL_SPECIES, conservation_df, focal_TRGs, DATA_DIR)
    print("\n")






    # Write result to file
    conservation_df.to_csv(os.path.join(OUTPUT_DIR, f"conservation_df_{FOCAL_SPECIES}.tsv"), sep="\t")

