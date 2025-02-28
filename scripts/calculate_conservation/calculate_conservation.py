import os
import glob
import re
import random
import subprocess
import tempfile
import pandas as pd
import Bio.SeqIO as SeqIO
from Bio.Seq import Seq
import concurrent.futures


OUTPUT_DIR = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/calculate_conservation/"
GENOMES_LIST = "/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/genera_archaea/genomes_list.txt"
GENERA_DIR = "/home/eliott.tempez/Documents/archaea_data/genera/out/"
DATA_DIR = "/home/eliott.tempez/Documents/archaea_data/complete_122/"
FOCAL_SPECIES = "GCA_020386975@Thermococcus_bergensis_T7324"
TRG_RANK = 7.0
NCPUS = 12



############################## FUNCTIONS ##############################
def initialise_conservation_df(genomes):
    """Initialise the dataframe with empty rows and columns"""
    col = ["n_trg", "n_cds", "n_intergenic", "n_f1", "n_f2"]
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
    if len(focal_TRGs) > 1000:
        focal_TRGs = random.sample(focal_TRGs, 1000)
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
    if len(focal_CDSs) > 1000:
        focal_CDSs_names = random.sample(focal_CDSs_names, 1000)
    focal_CDSs = {k: focal_CDSs[k] for k in focal_CDSs_names}
    return focal_CDSs


def extract_intergenic(species):
        """Extract all intergenic sequences"""
        # Make sure the file exists and read it
        gbk_pattern = os.path.join(DATA_DIR, "annotation", species + ".gb*")
        gbk_files = [f for f in glob.glob(gbk_pattern) if not f.endswith(".fai")]
        if gbk_files:
            gbk_file = gbk_files[0]
            gbk_content = list(SeqIO.parse(gbk_file, "genbank"))
        else:
            raise FileNotFoundError(f"No file matching pattern {gbk_pattern}")

        # Extract data for each contig
        intergenic_segments = {}
        for record in gbk_content:
            contig_length = len(record.seq)
            intergenic_segments[record.name] = []
            # Get all nucl positions in CDSs from the genbank file
            coding_loci = []
            for feature in record.features:
                if feature.type == "gene":
                    # Make sure we don't include the bugged CDS that is as long as the genome
                    if feature.location.end - feature.location.start < contig_length:
                        coding_loci += list(range(feature.location.start, feature.location.end))
            # Get all the integers that are not in the list of coding loci
            full_set = set(range(min(coding_loci), max(coding_loci) + 1))
            noncoding_positions = list(full_set - set(coding_loci))
            # Convert to consecutive ranges
            noncoding_positions.sort()
            start = noncoding_positions[0]
            for i in range(1, len(noncoding_positions)):
                if noncoding_positions[i] != noncoding_positions[i - 1] + 1:
                    intergenic_segments[record.name].append((start, noncoding_positions[i - 1] + 1))
                    start = noncoding_positions[i]
            intergenic_segments[record.name].append((start, noncoding_positions[-1] + 1))

        # Extract the corresponding sequences from fasta file
        fa_file = os.path.join(DATA_DIR, "fasta_renamed/" + species + ".fa")
        fa_content = list(SeqIO.parse(fa_file, "fasta"))
        intergenic_dict = {}
        i = 0
        for contig in fa_content:
            if contig.name in intergenic_segments:
                for j in range(len(intergenic_segments[contig.name])):
                    i += 1
                    start, end = intergenic_segments[contig.name][j]
                    intergenic_dict[f"{contig.name}_{i}"] = contig.seq[start:end]
        return intergenic_dict


def cut_chunks(seqs, length):
    """Cut the CDSs in lengths-long chunks randomly"""
    for k, v in seqs.items():
        len_over = len(v) - length
        if len_over == 0:
            len_over = 1
        # Chose chunk randomly in-frame
        possible_starts = list(range(0, len_over, 3))
        start = random.choice(possible_starts)
        seqs[k] = v[start:(start + length)]
    return seqs


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
    elif colname == "n_intergenic":
        # Extract intergenic sequences
        intergenic_dict = extract_intergenic(species)
        # Write the sequences to a temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as db_file:
            for gene, seq in intergenic_dict.items():
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

    # Check the file isn't blank
    with open(output_file_path, "r") as f:
        content = f.read().strip()
    if content == "":
        os.remove(query_file_path)
        os.remove(output_file_path)
        return 0, None

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
    # Keep homolog sequences
    if keep_homologs:
        homologs = {row["qseqid"]: row["sseqid"] for i, row in result.iterrows()}
        return len(result), homologs
    # Remove the temporary files
    os.remove(query_file_path)
    os.remove(output_file_path)
    # Return number of matches
    return len(result), None


def process_conservation_for_species(species, focal_sp, conservation_df, colname, query_sequences, db=None):
    """Process the conservation for the focal species and add to the dataframe"""
    # Keep homolog sequences for the CDSs
    homologs, keep_homologs = None, False
    if colname == "n_cds":
        keep_homologs = True

    if species == focal_sp:
        conservation_df.loc[species, colname] = len(query_sequences)
        return conservation_df, homologs
    
    # For all the other species, run BLAST for all the sequences
    # Get pre-established db if exists
    db_sp = db[species] if db else None
    db_fasta_file = get_db(species, colname, db_sp)
    blast_type = "blastn" if colname == "n_intergenic" else "blastp"
    # run blast
    nb_matches, homologs = run_blast(query_sequences, db_fasta_file, blast_type, keep_homologs, db_sp)
    print(f"Species {species}: {nb_matches} matches")
    conservation_df.loc[species, colname] = nb_matches
    
    # Delete db file if it was created
    if colname in ["n_f1", "n_f2", "n_intergenic"]:
        os.remove(db_fasta_file)

    return conservation_df, homologs 


def process_conservation_parallel(focal_sp, conservation_df, colname, query_sequences, db=None):
    """Process the conservation for all species in parallel"""
    homologs, keep_homologs = None, False
    if colname == "n_cds":
        keep_homologs = True
        homologs = {}

    with concurrent.futures.ThreadPoolExecutor(max_workers=NCPUS) as executor:
        futures = {executor.submit(process_conservation_for_species, species, focal_sp, conservation_df, colname, query_sequences, db): species for species in conservation_df.index}
        for future in concurrent.futures.as_completed(futures):
            species = futures[future]
            try:
                conservation_df, hom = future.result()
                if keep_homologs:
                    homologs[species] = hom
            except Exception as exc:
                print(f"Species {species} generated an exception: {exc}")

    return conservation_df, homologs







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
    """print("Calculating TRG conservation...")
    # Extract the name of 1000 (at most) TRGs of the focal species
    focal_TRG_genes = extract_focal_TRGs(FOCAL_SPECIES, TRG_RANK)
    # Get the protein sequences for each TRG
    focal_TRGs = get_seqs_from_gene_names(FOCAL_SPECIES, focal_TRG_genes, "faa")
    print(f"Extracted {len(focal_TRGs)} TRGs for the focal species {FOCAL_SPECIES}\n")
    # Calculate the TRG conservation and add to dataframe
    conservation_df, _ = process_conservation_parallel(FOCAL_SPECIES, conservation_df, "n_trg", focal_TRGs)
    print("\n\n")"""



    ########################################
    ############ CDSs / +1 / +2 ############
    ########################################
    print("Calculating CDS conservation...")
    # Extract 1000 CDss of the focal species
    focal_CDSs_nt = extract_focal_CDSs(FOCAL_SPECIES)
    print(f"Extracted {len(focal_CDSs_nt)} CDSs for the focal species {FOCAL_SPECIES}\n")
    # Cut 302-long chunks randomly
    focal_CDSs_nt = cut_chunks(focal_CDSs_nt, 302)
    # Translate to get 100 aa-long protein sequences
    focal_CDSs = {k: v[:300].translate(table=11) for k, v in focal_CDSs_nt.items()}
    # Calculate the CDS conservation and add to dataframe + keep homolog sequences
    conservation_df, homologs = process_conservation_parallel(FOCAL_SPECIES, conservation_df, "n_cds", focal_CDSs)
    print("\n\n")


    # Frames 1 and 2
    # Translate alt frames to get 100 aa-long protein sequences
    focal_CDS_frame_1 = {k: v[1:301].translate(table=11) for k, v in focal_CDSs_nt.items()}
    focal_CDS_frame_2 = {k: v[2:302].translate(table=11) for k, v in focal_CDSs_nt.items()}
    print("Calculating conservation for +1 frame...")
    # Calculate the conservation for +1 frame and add to dataframe
    conservation_df, _ = process_conservation_parallel(FOCAL_SPECIES, conservation_df, "n_f1", focal_CDS_frame_1, homologs)
    print("\n\n")
    print("Calculating conservation for +2 frame...")
    # Calculate the conservation for +2 frame and add to dataframe
    conservation_df, _ = process_conservation_parallel(FOCAL_SPECIES, conservation_df, "n_f2", focal_CDS_frame_2, homologs)
    print("\n\n")



    ########################################
    ############## Intergenic ##############
    ########################################
    """print("Calculating intergenic conservation...")
    # Extract 1000 intergenic 100 sequences of the focal species
    focal_intergenic_nt = extract_intergenic(FOCAL_SPECIES)
    # Keep only the intergenic sequences of at least 100 nt and cut chunks
    focal_intergenic = {k: v for k, v in focal_intergenic_nt.items() if len(v) >= 100}
    focal_intergenic = cut_chunks(focal_intergenic, 100)
    # Keep 1000 at most
    focal_intergenic_names = list(focal_intergenic.keys())
    if len(focal_intergenic) > 1000:
        focal_intergenic_names = random.sample(focal_intergenic_names, 1000)
    focal_intergenic = {k: focal_intergenic[k] for k in focal_intergenic_names}
    print(f"Extracted {len(focal_intergenic)} intergenic sequences for the focal species {FOCAL_SPECIES}\n")
    # Calculate the intergenic conservation and add to dataframe
    conservation_df, _ = process_conservation_parallel(FOCAL_SPECIES, conservation_df, "n_intergenic", focal_intergenic)
    print("\n\n")"""



    # Write result to file
    conservation_df.to_csv(os.path.join(OUTPUT_DIR, f"conservation_df_{FOCAL_SPECIES}.tsv"), sep="\t")

