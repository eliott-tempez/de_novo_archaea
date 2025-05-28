
import sys
import os
import glob
import re
import subprocess
import tempfile
import random
import multiprocessing
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.SeqUtils.ProtParam import ProteinAnalysis


SORTED_AA = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
OUT_DIR = "out/"
from my_functions.paths import DENSE_DIR, GENERA_DIR, GENOMES_LIST, CDS_DIR, FA_DIR
TRG_RANK = 7.0
GOOD_CANDIDATES_ONLY = True


def extract_cds_sequences(genome):
    cds_dict = {}
    faa_file = os.path.join(CDS_DIR, f"{genome}_CDS.fna")
    if not os.path.exists(faa_file):
        raise FileNotFoundError(f"No file {faa_file}")
    for record in SeqIO.parse(faa_file, "fasta"):
        if record.name not in cds_dict:
            cds_dict[record.name] = {}
        cds_dict[record.name]["sequence"] = record.seq
    return cds_dict


def extract_trg_names(focal_species, trg_threshold):
    """Extract TRGs of the focal species"""
    # Make sure the file exists and read it
    trg_pattern = os.path.join(GENERA_DIR, focal_species, "*gene_ages.tsv")
    trg_files = glob.glob(trg_pattern)
    if trg_files:
        trg_file = trg_files[0]
        trg_df = pd.read_csv(trg_file, sep="\t", header=0)
    else:
        raise FileNotFoundError(f"No file matching pattern {trg_pattern}")
    # Extract the TRGs
    focal_TRGs = trg_df[trg_df["rank"] >= trg_threshold]["#gene_mRNA"].tolist()
    return focal_TRGs


def extract_denovo_names(focal_species, use_good_candidates=False):
    """Extract the names of the de novo genes"""
    # Make sure the file exists and read it
    if not use_good_candidates:
        denovo_file = os.path.join(DENSE_DIR, focal_species, "denovogenes.tsv")
        if os.path.exists(denovo_file):
            denovo_df = pd.read_csv(denovo_file, sep="\t", header=0)
        else:
            raise FileNotFoundError(f"No file {denovo_file}")
        # Extract the names
        denovo_names = denovo_df["CDS"].tolist()
    else:
        denovo_names = []
        denovo_file = "good_candidates.txt"
        with open(denovo_file, "r") as f:
            for line in f:
                denovo_names.append(line.strip())
    return denovo_names


def get_species_gc_content(genome):
    fa_file = os.path.join(FA_DIR, genome + ".fa")
    if not os.path.exists(fa_file):
        raise FileNotFoundError(f"No file {fa_file}")
    seq = ""
    for record in SeqIO.parse(fa_file, "fasta"):
        seq += str(record.seq)
    return GC(seq)


def get_species_iorf_gc(genome):
    iorfs_dict = {}
    concat_seq = ""

    # Get iorfs from file
    with open("iorfs.txt", "r") as f:
        for line in f:
            if line.startswith(f">{genome}"):
                iorfs_str = f.readline().strip()
                break
    iorfs = iorfs_str.split()

    for i, segment in enumerate(iorfs):
        # Get the GC %
        seq = str(segment)
        concat_seq += seq
        # Add the iorfs to the dict
        iorfs_dict[f"{genome}_iorf_{i}"] = {"sequence": Seq(seq)}
    
    return iorfs_dict, GC(concat_seq)


def get_hcas(cds_names, all_cdss):
    #### HCA ####
    # Create temp fasta file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix="faa") as faa:
        for cds in cds_names:
            faa.write(f">{cds}\n{all_cdss[cds]['sequence'].translate(table=11)}\n")
        faa_file_path = faa.name
    # Create empty temp result file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix="tab") as result:
        result_file_path = result.name

    # Run pyHCA
    result = subprocess.run([
        "hcatk", "segment", "-i",
        faa_file_path, "-o", result_file_path,
    ])
    # Keep only useful lines and columns
    result_lines = []
    with open(result_file_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                args = line[1:].strip().split()
                cds_name = args[0]
                hca = float(args[-1])
                result_lines.append([cds_name, hca])

    # Read the result file
    hca_df = pd.DataFrame(result_lines, columns=["Seq_ID", "HCA"])

    # Remove the temp files
    os.remove(faa_file_path)
    os.remove(result_file_path)

    return hca_df


def calculate_proportion_of_seq_disordered(iupred_scores):
    # Imported from ORFmine's ORFold
    count_seg_tmp = 0
    count_agg_seg = 0
    for i,pos in enumerate(iupred_scores):
        if pos > 0.5:
            count_seg_tmp += 1 
        elif pos <= 0.5 and count_seg_tmp >=5:
            count_agg_seg = count_agg_seg + count_seg_tmp
            count_seg_tmp = 0
        elif pos <= 0.5 and count_seg_tmp <5:
            count_seg_tmp = 0 
            continue
        if i == len(iupred_scores) -1:
            if count_seg_tmp >=5:
                count_agg_seg = count_agg_seg + count_seg_tmp
            else:
                continue
    return(round(count_agg_seg/len(iupred_scores),3)) 


def calculate_proportion_of_seq_aggregable(b_aggregation):
    # Imported from ORFmine's ORFold
    count_seg_tmp = 0
    count_agg_seg = 0
    for i,pos in enumerate(b_aggregation):
        if pos > 5.0:
            count_seg_tmp += 1 
        elif pos <= 5.0 and count_seg_tmp >=5:
            count_agg_seg = count_agg_seg + count_seg_tmp
            count_seg_tmp = 0
        elif pos <= 5.0 and count_seg_tmp <5:
            count_seg_tmp = 0 
            continue
        if i == len(b_aggregation) -1:
            if count_seg_tmp >=5:
                count_agg_seg = count_agg_seg + count_seg_tmp
            else:
                continue
    return(round(count_agg_seg/len(b_aggregation),3))


def get_iupred(cds, all_cdss):
    aa_seq = str(all_cdss[cds]["sequence"].translate(table=11))
    # Create temp fasta file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix="faa") as faa:
        faa.write(f">{cds}\n{aa_seq}\n")
        faa_file_path = faa.name
    # Create empty temp result file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix="tab") as result:
        result_file_path = result.name

    # Run iupred
    with open(result_file_path, "w") as result_file:
        result = subprocess.run([
            "python", "iupred2a.py", "-a", faa_file_path, "short"
        ], stdout=result_file)
    # Get the score
    iupred_scores = []
    with open(result_file_path, "r") as f:
        for line in f:
            if not line.startswith("#"):
                args = line.strip().split("\t")
                if len(args) > 1:
                    iupred_scores.append(float(args[2]))
    disord = calculate_proportion_of_seq_disordered(iupred_scores)

    # Remove the temp files
    os.remove(faa_file_path)
    os.remove(result_file_path)
    return disord
    


def get_hca(all_hcas, cds_name):
    orfold_line = all_hcas[all_hcas["Seq_ID"] == cds_name]
    hca = orfold_line["HCA"].values[0]
    return hca


def get_tango(cds, all_cdss):
    # Paths
    genome = all_cdss[cds]["genome"]
    output_prefix = f"tango_results_{cds}_{genome}"
    current_dir = os.getcwd()
    tango_path = os.path.join(current_dir, "tango_x86_64_release")
    # Command arguments
    aa_seq = re.sub(r"[\*]", "", str(all_cdss[cds]["sequence"].translate(table=11)))
    args = f'ct="N" nt="N" ph="7.4" te="298" io="0.1" seq="{aa_seq}"'
    # Run Tango
    result = subprocess.run(f"{tango_path} {output_prefix} {args}", capture_output=True, text=True, shell=True, errors="replace")
    # Extract the score
    aggreg_scores = []

    if not os.path.exists(f"{output_prefix}.txt"):
        return np.nan
    with open(f"{output_prefix}.txt", "r") as f:
        for line in f:
            # Remove all blankspaces
            line = re.sub(r"\s+", " ", line)
            line_lst = line.strip().split()
            # read values
            b_aggreg = line_lst[5]
            # Keep only if floats
            try:
                b_aggreg = float(b_aggreg)
                aggreg_scores.append(b_aggreg)
            except ValueError:
                continue
    # Remove the temp file
    os.remove(f"{output_prefix}.txt")
    # Calculate the proportion of aggregable sequences
    aggreg = calculate_proportion_of_seq_aggregable(aggreg_scores)
    return aggreg


def get_common_aa_use(aa_use):
    # Types of aa
    polar_aa = ["S", "T", "N", "Q"]
    hydrophobic_aa = ["V", "I", "L", "M", "F", "W", "Y"]
    positive_aa = ["K", "R", "H"]
    negative_aa = ["D", "E"]
    proline_glycine_aa = ["P", "G"]
    alanine_aa = ["A"]
    cysteine_aa = ["C"]
    # Calculate the proportion of each type of aa
    polar_use = sum([aa_use[aa] for aa in polar_aa if aa in aa_use])
    hydrophobic_use = sum([aa_use[aa] for aa in hydrophobic_aa if aa in aa_use])
    positive_use = sum([aa_use[aa] for aa in positive_aa if aa in aa_use])
    negative_use = sum([aa_use[aa] for aa in negative_aa if aa in aa_use])
    proline_glycine_use = sum([aa_use[aa] for aa in proline_glycine_aa if aa in aa_use])
    alanine_use = sum([aa_use[aa] for aa in alanine_aa if aa in aa_use])
    cysteine_use = sum([aa_use[aa] for aa in cysteine_aa if aa in aa_use])
    return [polar_use, hydrophobic_use, positive_use, negative_use, proline_glycine_use, alanine_use, cysteine_use]


def process_cds(cds):
    genome = all_cdss[cds]["genome"]
    # Extract nucleotide sequence and calculate GC-related measures
    nuc_seq = all_cdss[cds]["sequence"]
    gc_seq = GC(nuc_seq)
    gc_species = all_cdss[cds]["gen_gc"]
    inter_gc_species = all_cdss[cds]["intergenic_gc"]
    gc_rate = gc_seq / gc_species
    inter_gc_rate = gc_seq / inter_gc_species
    # Translate nucleotide to protein sequence and remove stop/unknowns
    prot_seq = re.sub(r"[\*X]", "", str(nuc_seq.translate(table=11)))
    # Analyse protein properties
    analysis = ProteinAnalysis(prot_seq)
    aromaticity = analysis.aromaticity()
    instability = analysis.instability_index()
    flexibility = analysis.flexibility()
    mean_flexibility = (sum(flexibility) / len(flexibility)) if len(flexibility) != 0 else np.nan
    hydropathy = analysis.gravy()
    # Get sequence length
    length = len(nuc_seq)
    # HCA score from precomputed DataFrame
    hca = get_hca(all_hcas, cds)
    # Calculate disorder and aggregation
    disord = get_iupred(cds, all_cdss)
    aggreg = get_tango(cds, all_cdss)
    # Collect the basic result values
    result = [genome, cds, gc_rate, aromaticity, instability, mean_flexibility, hydropathy, length, hca, disord, aggreg, inter_gc_rate, gc_species, inter_gc_species]
    # Amino acid usage
    aa_use = analysis.amino_acids_percent
    sorted_aa_use = {key: value for key, value in sorted(aa_use.items())}
    for aa in sorted_aa_use:
        result.append(aa_use[aa])
    # Amino acid group proportions (e.g. polar, hydrophobic)
    result += get_common_aa_use(sorted_aa_use)
    # Annotate type of CDS (cds, trg, denovo, iorf)
    if cds in cds_names:
        result.append("cds")
    elif cds in trg_names:
        result.append("trg")
    elif cds in denovo_names:
        result.append("denovo")
    elif cds in iorf_names:
        result.append("iorf")
    return result





if __name__ == "__main__":
    # Get the list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]
    denovo_names, trg_names, cds_names, iorf_names = [], [], [], []
    all_cdss = {}
    all_values = []


    # Extract the cds info
    print("Extracting all CDSs...\n")
    for genome in genomes:
        # Get the species gc content
        genome_gc = get_species_gc_content(genome)
        iorfs_dict, intergenic_gc = get_species_iorf_gc(genome)
        # Extract de novo names
        if not GOOD_CANDIDATES_ONLY:
            denovo_names += extract_denovo_names(genome)
        # Extract TRG names
        trg_names += extract_trg_names(genome, TRG_RANK)
        # Extract info for all CDSs
        all_cds_gen = extract_cds_sequences(genome)
        # Extract info for iorfs
        iorf_names += iorfs_dict.keys()
        all_cds_gen.update(iorfs_dict)
        # Add the genome name
        for cds in all_cds_gen:
            all_cds_gen[cds]["genome"] = genome
            all_cds_gen[cds]["gen_gc"] = genome_gc
            all_cds_gen[cds]["intergenic_gc"] = intergenic_gc
        all_cdss.update(all_cds_gen)
        cds_names += all_cds_gen.keys()

    if GOOD_CANDIDATES_ONLY:
        denovo_names = extract_denovo_names(genome, True)

    # Drop duplicates
    cds_names = list(set(cds_names) - set(trg_names) - set(iorf_names))
    trg_names = list(set(trg_names) - set(denovo_names))

    # Sample
    """cds_names = random.sample(cds_names, 5000)
    trg_names = random.sample(trg_names, 5000)
    #denovo_names = random.sample(denovo_names, 10)
    iorf_names = random.sample(iorf_names, 5000)"""

    # Calculate descriptors for all cdss
    all_cds_names = denovo_names + trg_names + cds_names + iorf_names
    n = len(all_cdss)
    i = 0

    # Extract all hcas
    all_hcas = get_hcas(all_cds_names, all_cdss)

    results = []
    print("Starting parallel processing...\n")
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        for idx, result in enumerate(pool.imap(process_cds, all_cds_names), start=1):
            results.append(result)
            if idx % (n // 10) == 0:  # Print progress every 10% of total sequences
                print(f"{idx}/{n} sequences analysed")

    print("\nDone!")

    # Save the results
    df = pd.DataFrame(results, columns=["genome", "cds", "gc_rate", "aromaticity", "instability", "mean_flexibility", "hydropathy", "length", "hca", "disord", "aggreg", "inter_gc_rate", "gc_species", "inter_gc_species"] + [f"{a}_use" for a in SORTED_AA] + ["polar_use", "hydrophobic_use", "positive_use", "negative_use", "proline-glycine_use", "alanine_use", "cysteine_use"] + ["type"])
    df.to_csv(os.path.join(OUT_DIR, "sequence_features_good_candidates_all.csv"), sep="\t", index=False)