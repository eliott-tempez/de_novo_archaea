
import sys
import os
import glob
import re
import subprocess
import tempfile
import random
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils.ProtParam import ProteinAnalysis


sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
OUT_DIR = "out/"
from my_functions.paths import DENSE_DIR, GENERA_DIR, GENOMES_LIST, CDS_DIR, FA_DIR
from my_functions.genomic_functions import extract_iorfs
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
    iorfs = extract_iorfs(genome)
    concat_seq = ""
    for segment in iorfs:
        seq = str(segment)
        concat_seq += seq
    return GC(concat_seq)


def get_hcas(cds_names, all_cdss):
    current_dir = os.getcwd()
    # Get the fold potential
    # Create temp fasta file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix="faa") as faa:
        for cds in cds_names:
            faa.write(f">{cds}\n{all_cdss[cds]['sequence']}\n")
        faa_file_path = faa.name
    # Create empty temp result file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix="tab") as result:
        result_file_path = result.name

    # Run pyHCA
    result = subprocess.run([
        "hcatk", "segment", "-i",
        faa_file_path, "-o", result_file_path,
    ])
    with open(result_file_path, "r") as f:
        content = f.read()
        print(f"pyHCA output:\n{content}\n")
    # Keep only useful lines and columns
    result_lines = []
    with open(result_file_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                args = line[1:].strip().split()
                cds_name = args[0]
                hca = args[-1]
                result_lines.append([cds_name, hca])

    # Read the result file
    hca_df = pd.DataFrame(result_lines, columns=["Seq_ID", "HCA"])

    return hca_df


def get_orfold_descript(all_hcas, cds_name):
    orfold_line = all_hcas[all_hcas["Seq_ID"] == cds]
    hca = orfold_line["HCA"].values[0]
    #disord = orfold_line["Disord"].values[0]
    #aggreg = orfold_line["Aggreg"].values[0]
    disord, aggreg = 0, 0
    return hca, disord, aggreg




if __name__ == "__main__":
    # Get the list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]
    denovo_names, trg_names, cds_names = [], [], []
    all_cdss = {}
    all_values = []


    # Extract the cds info
    print("Extracting all CDSs...\n")
    for genome in genomes:
        # Get the species gc content
        genome_gc = get_species_gc_content(genome)
        intergenic_gc = get_species_iorf_gc(genome)

        # Extract de novo names
        if not GOOD_CANDIDATES_ONLY:
            denovo_names += extract_denovo_names(genome)
        # Extract TRG names
        trg_names += extract_trg_names(genome, TRG_RANK)
        # Extract info for all CDSs
        all_cds_gen = extract_cds_sequences(genome)
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
    cds_names = list(set(cds_names) - set(trg_names))
    trg_names = list(set(trg_names) - set(denovo_names))

    # Sample
    cds_names = random.sample(cds_names, 2)
    trg_names = random.sample(trg_names, 2)
    denovo_names = random.sample(denovo_names, 2)

    # Calculate descriptors for all cdss
    all_cds_names = denovo_names + trg_names + cds_names
    n = len(all_cdss)
    i = 0

    # Extract all hcas
    all_hcas = get_hcas(all_cds_names, all_cdss)
    print(all_hcas)

    results = []
    for cds in all_cds_names:
        i += 1
        genome = all_cdss[cds]["genome"]
        # Extract GC rate
        nuc_seq = all_cdss[cds]["sequence"]
        gc_seq = GC(nuc_seq)
        gc_species = all_cdss[cds]["gen_gc"]
        inter_gc_species = all_cdss[cds]["intergenic_gc"]
        gc_rate = gc_seq / gc_species
        inter_gc_rate = gc_seq / inter_gc_species
        # Extract protein info
        prot_seq = re.sub(r"[\*X]", "", str(nuc_seq.translate()))
        analysis = ProteinAnalysis(prot_seq)
        aromaticity = analysis.aromaticity()
        instability = analysis.instability_index()
        flexibility = analysis.flexibility()
        mean_flexibility = sum(flexibility) / len(flexibility)
        hydropathy = analysis.gravy()
        # Extract sequence length
        length = len(nuc_seq)
        # Extract hca, disorder and aggregation
        #hca, disord, aggreg = get_orfold_descript(all_hcas, cds)
        hca = get_orfold_descript(all_hcas, cds)
        print(cds, hca)
        result = [genome, cds, gc_rate, aromaticity, instability, mean_flexibility, hydropathy, length, hca, disord, aggreg, inter_gc_rate, gc_species, inter_gc_species]

        # Extract aa use
        aa_use = analysis.amino_acids_percent
        # Sort dict alphabetically
        sorted_aa_use = {key: value for key, value in sorted(aa_use.items())}
        for aa in sorted_aa_use:
            result.append(aa_use[aa])

        # Add the type of cds
        if cds in cds_names:
            result.append("cds")
        elif cds in trg_names:
            result.append("trg")
        elif cds in denovo_names:
            result.append("denovo")
        
        if i % round(n/20) == 0:
            print(f"{i}/{n} cds analysed...")
        
        results.append(result)

    print("\nDone!")

    # Save the results
    df = pd.DataFrame(results, columns=["genome", "cds", "gc_rate", "aromaticity", "instability", "mean_flexibility", "hydropathy", "length", "hca", "disord", "aggreg", "inter_gc_rate", "gc_species", "inter_gc_species"] + [f"{a}_use" for a in list(sorted_aa_use.keys())] + ["type"])
    df.to_csv(os.path.join(OUT_DIR, "sequence_features_good_candidates_all.csv"), sep="\t", index=False)