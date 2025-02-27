import os
import glob
import re
import random
import subprocess
import tempfile
import pandas as pd
import Bio.SeqIO as SeqIO
from Bio.Seq import Seq


OUTPUT_DIR = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/calculate_conservation/"
GENOMES_LIST = "/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/genera_archaea/genomes_list.txt"
GENERA_DIR = "/home/eliott.tempez/Documents/archaea_data/genera/out/"
DATA_DIR = "/home/eliott.tempez/Documents/archaea_data/complete_122/"
FOCAL_SPECIES = "GCA_001433455@Thermococcus_barophilus_CH5"
TRG_RANK = 7.0
NCPUS = 8



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
    if len(focal_TRGs) > 100:
        focal_TRGs = random.sample(focal_TRGs, 100)
    return focal_TRGs


def extract_focal_CDSs(focal_species):
    """Extract at most 1000 CDSs of the focal species"""
    # Extract all CDSs >=302 nt of the focal species
    CDS_fasta_file = os.path.join(DATA_DIR, "CDS/" + focal_species + "_CDS.fna")
    CDS_fasta = SeqIO.parse(CDS_fasta_file, "fasta")
    focal_CDSs = {}
    for record in CDS_fasta:
        if len(record.seq) >= 302:
            focal_CDSs[record.name] = record.seq
    # Keep 1000 at most
    focal_CDSs_names = list(focal_CDSs.keys())
    if len(focal_CDSs) > 100:
        focal_CDSs_names = random.sample(focal_CDSs_names, 100)
    focal_CDSs = {k: focal_CDSs[k] for k in focal_CDSs_names}
    return focal_CDSs


def cut_chunks(focal_CDSs):
    """Cut the CDSs in 302-long chunks"""
    for k, v in focal_CDSs.items():
        len_over_302 = len(v) - 302
        # Chose chunk randomly in-frame
        possible_starts = list(range(0, len_over_302, 3))
        start = random.choice(possible_starts)
        focal_CDSs[k] = v[start:start + 302]
    return focal_CDSs


def get_db(species, colname, db):
    """Get the database fasta file for the species depending on the type of sequences"""
    if colname == "n_trg" or colname == "n_cds":
        db_fasta_file = os.path.join(DATA_DIR, "CDS/" + species + "_CDS.faa")
    elif colname == "n_f1":
        # Get the dna sequences for the matching genes
        gene_names = list(db.values())
        gene_dict = get_seqs_from_gene_names(species, gene_names, "fna")
        # Keep +1 frame only and translate to protein
        gene_dict = {k: v[1:-2].translate(table=11) for k, v in gene_dict.items()}
        # Write the sequences to a temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as db_file:
            for gene, seq in gene_dict.items():
                db_file.write(f">{gene}\n{seq}\n")
            db_fasta_file = db_file.name
    elif colname == "n_f2":
        # Get the dna sequences for the matching genes
        gene_names = list(db.values())
        gene_dict = get_seqs_from_gene_names(species, gene_names, "fna")
        # Keep +2 frame only and translate to protein
        gene_dict = {k: v[2:-1].translate(table=11) for k, v in gene_dict.items()}
        # Write the sequences to a temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as db_file:
            for gene, seq in gene_dict.items():
                db_file.write(f">{gene}\n{seq}\n")
            db_fasta_file = db_file.name
    
    return db_fasta_file


def run_blast(query_sequences, db_faa_file, blast_type, keep_homologs, db):
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
    try:
        output = subprocess.run([blast_type, "-query", query_file_path, "-subject", db_faa_file, "-out", output_file_path, "-outfmt", "6 qseqid sseqid qlen evalue qcovs", "-num_threads", str(NCPUS)], capture_output=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"BLAST command failed with return code {e.returncode}: {e.stderr.decode()}")
        raise
    # Read the output file
    result = pd.read_csv(output_file_path, sep="\t", header=None)
    result.columns = ["qseqid", "sseqid", "qlen", "evalue", "qcov"]
    # Keep only rows for which qcov > 50 % and evalue < 1e-3
    result = result[(result["qcov"] > 50) & (result["evalue"] < 1e-3)]
    # If we have a match couple already known
    if db:
        # Keep only the matches in the right CDS
        result = result[result.apply(lambda row: db.get(row ["qseqid"]) == row ["sseqid"], axis=1)]
    # Keep only the best hit for each query
    result = result.sort_values("evalue").drop_duplicates("qseqid")
    # Remove the temporary files
    os.remove(query_file_path)
    os.remove(output_file_path)
    # Keep homolog sequences
    if keep_homologs:
        homologs = {row["qseqid"]: row["sseqid"] for i, row in result.iterrows()}
        return len(result), homologs
    # Return number of matches
    return len(result), None


def process_conservation(focal_sp, conservation_df, colname, query_sequences, db=None):
    """Process the conservation for the focal species and add to the dataframe"""
    # Keep homolog sequences for the CDSs
    homologs, keep_homologs = None, False
    if colname == "n_cds":
        keep_homologs = True
        homologs = {}

    # For each species
    for species in conservation_df.index:
        # For the focal species, add the number of sequences sampled
        if species == focal_sp:
            conservation_df.loc[species, colname] = len(query_sequences)
            continue

        # For all the other species, run BLAST for all the sequences
        # Get pre-established db if exists
        db_sp = db[species] if db else None
        db_fasta_file = get_db(species, colname, db_sp)
        nb_matches, hom = run_blast(query_sequences, db_fasta_file, "blastp", keep_homologs, db_sp)
        print(f"Species {species}: {nb_matches} matches")
        conservation_df.loc[species, colname] = nb_matches
        if keep_homologs:
            homologs[species] = hom
        
        # Delete db file if it was created
        if colname in ["n_f1", "n_f2", "n_intergenic"]:
            os.remove(db_fasta_file)

    return conservation_df, homologs 



    










































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
                    prot_seq = nucl_seq.translate(table=11)
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
            








############################## MAIN ##############################
if __name__ == "__main__":
    # Read the list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]
    # Initiate the dataframe
    conservation_df = initialise_conservation_df(genomes)
    print("\n")


    """########################################
    ################# TRGs #################
    ########################################
    print("Calculating TRG conservation...")
    # Extract the name of 1000 (at most) TRGs of the focal species
    focal_TRG_genes = extract_focal_TRGs(FOCAL_SPECIES, TRG_RANK)
    # Get the protein sequences for each TRG
    focal_TRGs = get_seqs_from_gene_names(FOCAL_SPECIES, focal_TRG_genes, "faa")
    print(f"Extracted {len(focal_TRGs)} TRGs for the focal species {FOCAL_SPECIES}\n")
    # Calculate the TRG conservation and add to dataframe
    conservation_df, _ = process_conservation(FOCAL_SPECIES, conservation_df, "n_trg", focal_TRGs)
    print("\n\n")"""



    ########################################
    ############ CDSs / +1 / +2 ############
    ########################################
    print("Calculating CDS conservation...")
    # Extract 1000 CDss of the focal species
    focal_CDSs_nt = extract_focal_CDSs(FOCAL_SPECIES)
    print(f"Extracted {len(focal_CDSs_nt)} CDSs for the focal species {FOCAL_SPECIES}\n")
    # Cut 302-long chunks randomly
    focal_CDSs_nt = cut_chunks(focal_CDSs_nt)
    # Translate to get 100 aa-long protein sequences
    focal_CDSs = {k: v[:300].translate(table=11) for k, v in focal_CDSs_nt.items()}
    # Calculate the CDS conservation and add to dataframe + keep homolog sequences
    conservation_df, homologs = process_conservation(FOCAL_SPECIES, conservation_df, "n_cds", focal_CDSs)
    print("\n\n")


    # Frames 1 and 2
    # Translate alt frames to get 100 aa-long protein sequences
    focal_CDS_frame_1 = {k: v[1:301].translate(table=11) for k, v in focal_CDSs_nt.items()}
    focal_CDS_frame_2 = {k: v[2:302].translate(table=11) for k, v in focal_CDSs_nt.items()}
    print("Calculating conservation for +1 frame...")
    # Calculate the conservation for +1 frame and add to dataframe
    conservation_df, _ = process_conservation(FOCAL_SPECIES, conservation_df, "n_f1", focal_CDS_frame_1, homologs)
    print("\n\n")
    print("Calculating conservation for +2 frame...")
    # Calculate the conservation for +2 frame and add to dataframe
    conservation_df, _ = process_conservation(FOCAL_SPECIES, conservation_df, "n_f2", focal_CDS_frame_2, homologs)
    print("\n\n")


















    """########################################
    ############## Non-coding ##############
    ########################################
    print("Calculating non-coding conservation...")
    # Extract 1000 non-coding sequences of the focal species
    focal_noncoding = extract_focal_noncoding(FOCAL_SPECIES)
    print(f"Extracted {len(focal_noncoding)} non-coding sequences for the focal species {FOCAL_SPECIES}\n")
    # Calculate the non-coding conservation and add to dataframe
    conservation_df = calculate_noncoding_conservation(FOCAL_SPECIES, conservation_df, focal_noncoding)
    print("\n\n")"""



    # Write result to file
    conservation_df.to_csv(os.path.join(OUTPUT_DIR, f"conservation_df_{FOCAL_SPECIES}.tsv"), sep="\t")

