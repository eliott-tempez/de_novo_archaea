import os
import glob
import re
import random
import subprocess
import tempfile
import pandas as pd
import Bio.SeqIO as SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool
from functools import partial


OUTPUT_DIR = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/calculate_conservation/"
GENOMES_LIST = "/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/genera_archaea/genomes_list.txt"
GENERA_DIR = "/home/eliott.tempez/Documents/archaea_data/genera/out/"
DATA_DIR = "/home/eliott.tempez/Documents/archaea_data/complete_122/"
FOCAL_SPECIES = "GCA_001433455@Thermococcus_barophilus_CH5"
TRG_RANK = 7.0
NCPUS = 12



############################## FUNCTIONS ##############################
def initialise_conservation_df(genomes):
    """Initialise the dataframe with empty rows and columns"""
    col = ["n_trg", "n_cds", "n_noncoding", "n_f1", "n_f2"]
    conservation_df = pd.DataFrame(columns=col, index=genomes)
    return conservation_df



def get_seqs_from_gene_names(focal_sp, gene_names, extension_type):
    """Get the sequences of the genes from the gene names"""
    fasta_file = os.path.join(DATA_DIR, "CDS/" + focal_sp + "_CDS." + extension_type)
    fasta = SeqIO.parse(fasta_file, "fasta")
    seqs = {}
    for record in fasta:
        if record.name in gene_names:
            seqs[record.name] = record.seq
    return seqs



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
    result = subprocess.run([blast_type, "-query", query_file_path, "-subject", db_faa_file, "-out", output_file_path, "-outfmt", "6 qseqid sseqid qlen evalue qcovs", "-num_threads", NCPUS], capture_output=True)
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



def extract_focal_TRGs(focal_species, trg_threshold):
    """Extract at most 1000 TRGs of the focal species"""
    # Make sure the file exists and read it
    trg_pattern = os.path.join(GENERA_DIR, focal_species, "*gene_ages.tsv")
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



def calculate_conservation_parallel(species, focal_sp, conservation_df, sequences, blast_type):
    if species == focal_sp:
        conservation_df.loc[species, blast_type] = len(sequences)
        return conservation_df.loc[species]
    
    db_fa_file = os.path.join(DATA_DIR, "fasta_renamed/" + species + ".fa")
    nb_matches = run_blast(sequences, db_fa_file, "tblastn")
    print(f"Species {species}: {nb_matches} matches for {blast_type}")
    conservation_df.loc[species, blast_type] = nb_matches
    return conservation_df.loc[species]



def process_conservation(focal_sp, conservation_df, sequences, blast_type):
    with Pool(processes=NCPUS) as pool:
        func = partial(calculate_conservation_parallel, focal_sp=focal_sp, 
                       conservation_df=conservation_df, sequences=sequences, blast_type=blast_type)
        results = pool.map(func, conservation_df.index)
    for res in results:
        conservation_df.loc[res.name] = res
    return conservation_df
























def calculate_TRG_conservation(focal_sp, conservation_df, focal_TRGs):
    """Calculate the conservation of the TRGs for all the species"""
    # Get the protein sequences for the TRGs
    TRG_sequences = {}
    CDS_fasta_file = os.path.join(DATA_DIR, "CDS/" + focal_sp + "_CDS.faa")
    CDS_fasta = SeqIO.parse(CDS_fasta_file, "fasta")
    for record in CDS_fasta:
        if record.name in focal_TRGs:
            TRG_sequences[record.name] = record.seq
    if len(TRG_sequences) != len(focal_TRGs):
        print("Some TRGs are missing in the faa file")

    # For each species
    for species in conservation_df.index:
        # For the focal species, add the number of TRGs
        if species == focal_sp:
            conservation_df.loc[species, "n_trg"] = len(focal_TRGs)
            continue

        # For all the other species, run BLAST for all the TRGs
        db_fa_file = os.path.join(DATA_DIR, "fasta_renamed/" + species + ".fa")
        nb_matches = run_blast(TRG_sequences, db_fa_file, "tblastn")
        print(f"Species {species}: {nb_matches} matches")
        conservation_df.loc[species, "n_trg"] = nb_matches
    return conservation_df



def extract_focal_CDSs(focal_species):
    """Extract 1000 CDSs of the focal species"""
    # Extract all CDSs of the focal species
    CDS_fasta_file = os.path.join(DATA_DIR, "CDS/" + focal_species + "_CDS.faa")
    CDS_fasta = SeqIO.parse(CDS_fasta_file, "fasta")
    focal_CDSs = {}
    for record in CDS_fasta:
        focal_CDSs[record.name] = record.seq
    # Keep 1000 at most
    focal_CDSs_names = list(focal_CDSs.keys())
    if len(focal_CDSs) > 1000:
        focal_CDSs_names = random.sample(focal_CDSs_names, 1000)
    focal_CDSs = {k: focal_CDSs[k] for k in focal_CDSs_names}
    return focal_CDSs



def calculate_CDS_conservation(focal_sp, conservation_df, focal_CDSs):
    """Calculate the conservation of the CDSs for all the species"""
    # For each species
    for species in conservation_df.index:
        # For the focal species, add the number of CDSs
        if species == focal_sp:
            conservation_df.loc[species, "n_cds"] = len(focal_CDSs)
            continue
        # For all the other species, run BLAST for all the CDSs
        db_fa_file = os.path.join(DATA_DIR, "fasta_renamed/" + species + ".fa")
        nb_matches = run_blast(focal_CDSs, db_fa_file, "tblastn")
        print(f"Species {species}: {nb_matches} matches")
        conservation_df.loc[species, "n_cds"] = nb_matches
    return conservation_df



def extract_focal_noncoding(focal_sp):
        """Extract 1000 non-coding sequences of the focal species"""
        # Make sure the file exists and read it
        gbk_pattern = os.path.join(DATA_DIR, "annotation/" + focal_sp + ".gb*")
        gbk_files = glob.glob(gbk_pattern)
        if gbk_files:
            gbk_file = gbk_files[0]
            gbk_content = list(SeqIO.parse(gbk_file, "genbank"))
        else:
            raise FileNotFoundError(f"No file matching pattern {gbk_pattern}")

        # Extract data for each contig
        noncoding_segments = []
        for record in gbk_content:
            # Get all coding loci from the genbank file (cds ranges)
            coding_loci = []
            for feature in record.features:
                if feature.type == "CDS":
                    coding_loci.append((feature.location.start, feature.location.end))
            # Get the starting positions for all codons
            start_point_codons = []
            for cds_range in coding_loci:
                start_point_codons += list(range(cds_range[0], cds_range[1], 3))
            # order and delete duplicates
            start_point_codons.sort()
            full_set = set(range(start_point_codons[0], start_point_codons[-1] + 1))
            start_point_codons = set(start_point_codons)
            # Get all the integers that are not in the list of first codons
            noncoding_positions = list(full_set - start_point_codons)
            noncoding_positions.sort()
            # Keep only 300-nucl-long segments that are not in frame with a cds
            for i in noncoding_positions[:-300]:
                considered_segment_codons = [c for c in range(i, i + 300, 3)]
                if not any([c in start_point_codons for c in considered_segment_codons]):
                    noncoding_segments.append((record.name, i, i + 300))
        # Extract 1000 noncoding segments at most
        if len(noncoding_segments) > 1000:
            noncoding_segments = random.sample(noncoding_segments, 1000)

        # Extract the corresponding sequences from fasta file
        # convert list to dict
        noncoding_dict_loc = {}
        for i in range(len(noncoding_segments)):
            contig, start, end = noncoding_segments[i]
            if contig not in noncoding_dict_loc:
                noncoding_dict_loc[contig] = {}
            noncoding_dict_loc[contig][f"segment_{i}"] = (start, end)
        # Extract the sequences
        noncoding_dict_seq = {}
        fa_file = os.path.join(DATA_DIR, "fasta_renamed/" + focal_sp + ".fa")
        fa_content = list(SeqIO.parse(fa_file, "fasta"))
        for contig in fa_content:
            if contig.name in noncoding_dict_loc:
                for segment in noncoding_dict_loc[contig.name]:
                    start, end = noncoding_dict_loc[contig.name][segment]
                    # translate the sequence
                    nucl_seq = contig.seq[start:end]
                    prot_seq = nucl_seq.translate()
                    noncoding_dict_seq[segment] = prot_seq
        return noncoding_dict_seq



def calculate_noncoding_conservation(focal_sp, conservation_df, focal_noncoding):
    # For each species
    for species in conservation_df.index:
        # For the focal species, add the number of CDSs
        if species == focal_sp:
            conservation_df.loc[species, "n_noncoding"] = len(focal_noncoding)
            continue
        # For all the other species, run BLAST for all the non-coding sequences
        db_fa_file = os.path.join(DATA_DIR, "fasta_renamed/" + species + ".fa")
        nb_matches = run_blast(focal_noncoding, db_fa_file, "tblastn")
        print(f"Species {species}: {nb_matches} matches")
        conservation_df.loc[species, "n_noncoding"] = nb_matches
    return conservation_df
            


def extract_focal_altframes(focal_sp):
    # Extract CDS sequences in nucleotides form
    CDS_fasta_file = os.path.join(DATA_DIR, "CDS/" + focal_sp + "_CDS.fna")
    CDS_fasta = SeqIO.parse(CDS_fasta_file, "fasta")
    focal_CDSs = {}
    for record in CDS_fasta:
        focal_CDSs[record.name] = record.seq
    # Keep 1000 at most
    focal_CDSs_names = list(focal_CDSs.keys())
    if len(focal_CDSs) > 1000:
        focal_CDSs_names = random.sample(focal_CDSs_names, 1000)
    focal_CDSs = {k: focal_CDSs[k] for k in focal_CDSs_names}
    # Extract the proteic sequences for the +1 and +2 frames
    focal_f1 = {k: focal_CDSs[k][1:].translate() for k in focal_CDSs}
    focal_f2 = {k: focal_CDSs[k][2:].translate() for k in focal_CDSs}
    return focal_f1, focal_f2



def calculate_altframes_conservation(focal_sp, conservation_df, focal_f1, focal_f2):
    # For each species
    for species in conservation_df.index:
        # For the focal species, add the number of CDSs
        if species == focal_sp:
            conservation_df.loc[species, "n_f1"] = len(focal_f1)
            conservation_df.loc[species, "n_f2"] = len(focal_f2)
            continue
        # For all the other species, run BLAST for all the +1 and +2 frames
        db_fa_file = os.path.join(DATA_DIR, "fasta_renamed/" + species + ".fa")
        nb_matches_f1 = run_blast(focal_f1, db_fa_file, "tblastn")
        nb_matches_f2 = run_blast(focal_f2, db_fa_file, "tblastn")
        print(f"Species {species}: {nb_matches_f1} matches for +1 frame and {nb_matches_f2} matches for +2 frame")
        conservation_df.loc[species, "n_f1"] = nb_matches_f1
        conservation_df.loc[species, "n_f2"] = nb_matches_f2
    return conservation_df




############################## MAIN ##############################
if __name__ == "__main__":
    # Read the list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]
    # Initiate the dataframe
    conservation_df = initialise_conservation_df(genomes)
    print("\n")


    ########################################
    ################# TRGs #################
    ########################################
    print("Calculating TRG conservation...")
    # Extract the name of 1000 (at most) TRGs of the focal species
    focal_TRG_genes = extract_focal_TRGs(FOCAL_SPECIES, TRG_RANK)
    # Get the protein sequences for each TRG
    focal_TRGs = get_seqs_from_gene_names(FOCAL_SPECIES, focal_TRG_genes, "faa")
    print(f"Extracted {len(focal_TRGs)} TRGs for the focal species {FOCAL_SPECIES}\n")
    # Calculate the TRG conservation and add to dataframe
    conservation_df = process_conservation(FOCAL_SPECIES, conservation_df, focal_TRGs, "blastp")
    print("\n\n")



    ########################################
    ################# CDSs #################
    ########################################
    print("Calculating CDS conservation...")
    # Extract 1000 CDss of the focal species
    focal_CDSs = extract_focal_CDSs(FOCAL_SPECIES)



    print(f"Extracted {len(focal_CDSs)} CDSs for the focal species {FOCAL_SPECIES}\n")
    # Calculate the CDS conservation and add to dataframe
    conservation_df = calculate_CDS_conservation(FOCAL_SPECIES, conservation_df, focal_CDSs)
    print("\n\n")



    ########################################
    ############## Non-coding ##############
    ########################################
    print("Calculating non-coding conservation...")
    # Extract 1000 non-coding sequences of the focal species
    focal_noncoding = extract_focal_noncoding(FOCAL_SPECIES)
    print(f"Extracted {len(focal_noncoding)} non-coding sequences for the focal species {FOCAL_SPECIES}\n")
    # Calculate the non-coding conservation and add to dataframe
    conservation_df = calculate_noncoding_conservation(FOCAL_SPECIES, conservation_df, focal_noncoding)
    print("\n\n")



    ########################################
    ########### +1 and +2 frames ###########
    ########################################
    print("Calculating conservation for +1 and +2 frames...")
    # Extract 1000 sequences for +1 and +2 frames
    focal_f1, focal_f2 = extract_focal_altframes(FOCAL_SPECIES)
    print(f"Extracted {len(focal_f1)} sequences for the +1 frame and +2 frames for the focal species {FOCAL_SPECIES}\n")
    # Calculate the conservation for +1 and +2 frames and add to dataframe
    conservation_df = calculate_altframes_conservation(FOCAL_SPECIES, conservation_df, focal_f1, focal_f2)



    # Write result to file
    conservation_df.to_csv(os.path.join(OUTPUT_DIR, f"conservation_df_{FOCAL_SPECIES}.tsv"), sep="\t")

