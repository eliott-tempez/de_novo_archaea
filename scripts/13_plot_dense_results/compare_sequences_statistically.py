
import sys
import os
import glob
import re
import random
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction as GC
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np


sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from my_functions.paths import DENSE_DIR, GENERA_DIR, GENOMES_LIST, CDS_DIR, FA_DIR
from my_functions.genomic_functions import extract_iorfs
OUT_DIR = "out/"
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


def calculate_descriptors(descriptors, all_cdss, cds_names):
    descriptors_of_interest = {}
    current_dir = os.getcwd()
    # Get the fold potential
    # Create temp fasta file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix="faa") as faa:
        for cds in cds_names:
            faa.write(f">{cds}\n{all_cdss[cds]['sequence']}\n")
        faa_file_path = faa.name
    # Get the temp file directory
    faa_dir = os.path.dirname(faa_file_path)
    faa_basename = os.path.basename(faa_file_path)
    # Full container path to the temp file (inside the container)
    container_faa_path = f"/tmpdata/{faa_basename}"
    # Create output dir
    orfold_output_dir = os.path.join(current_dir, "orfold")
    os.makedirs(orfold_output_dir, exist_ok=True)
    # Get programs dir
    programs_dir = os.path.join(current_dir, "softwares")


    # Run ORFold
    container_path = "orfmine_latest.sif"
    # Bind both the fasta temp dir and the output dir
    result = subprocess.run([
    "singularity", "exec",
    "--bind", f"{faa_dir}:/database/tmpdata",
    "--bind", f"{orfold_output_dir}:/workdir/orfold",
    "--bind", f"{programs_dir}:/ORFmine/orfold_v1/orfold/softwares/",
    container_path,
    "orfold", "-faa", container_faa_path, "-options", "HIT"], text=True, capture_output=True)
    # Fix the broken output
    result_file = os.path.join(orfold_output_dir, faa_basename + ".tab")
    subprocess.run(["sed", "-i", r"s/[[:space:]]\+/;/g", result_file])

    # Result path (in the output dir now)
    orfold_result = pd.read_csv(result_file, sep=";", header=0)
    # Delete temp file
    os.remove(faa_file_path)
    os.remove(result_file)


    for cds in cds_names:
        # Add the info to the global dict if we don't already have it
        if cds not in descriptors:
            descriptors[cds] = {}
            # Extract GC rate
            nuc_seq = all_cdss[cds]["sequence"]
            gc_seq = GC(nuc_seq)
            gc_species = all_cdss[cds]["gen_gc"]
            inter_gc_species = all_cdss[cds]["intergenic_gc"]
            descriptors[cds]["gc_rate"] = gc_seq / gc_species
            descriptors[cds]["inter_gc_rate"] = gc_seq / inter_gc_species
            # Extract protein info
            prot_seq = re.sub(r"[\*X]", "", str(nuc_seq.translate()))
            analysis = ProteinAnalysis(prot_seq)
            descriptors[cds]["aromaticity"] = analysis.aromaticity()
            descriptors[cds]["instability"] = analysis.instability_index()
            flexibility = analysis.flexibility()
            descriptors[cds]["mean_flexibility"] = sum(flexibility) / len(flexibility)
            descriptors[cds]["hydropathy"] = analysis.gravy()
            # Extract sequence length
            descriptors[cds]['length'] = len(nuc_seq)
            # Extract hca, iupred, tango
            orfold_line = orfold_result[orfold_result["Seq_ID"] == cds]
            descriptors[cds]['hca'] = orfold_line["HCA"].values[0]
            descriptors[cds]["disord"] = orfold_line["Disord"].values[0]
            descriptors[cds]["aggreg"] = orfold_line["Aggreg"].values[0]
            # Extract aa use
            aa_use = analysis.amino_acids_percent
            for aa in aa_use:
                descriptors[cds][f"{aa}_use"] = aa_use[aa]

        # Add info to descriptor dict for samples
        for descriptor_name in descriptors[cds]:
            if descriptor_name not in descriptors_of_interest:
                descriptors_of_interest[descriptor_name] = []
            descriptors_of_interest[descriptor_name].append(descriptors[cds][descriptor_name])
    return descriptors, descriptors_of_interest


def get_median_diff(first_descript, second_descript):
    median_diff = {}
    for descriptor in first_descript:
        median_1 = np.median(first_descript[descriptor])
        median_2 = np.median(second_descript[descriptor])
        median_diff[descriptor] = median_1 - median_2
    return median_diff


def pool_cdss(first_descript, second_descript):
    first_pool = {}
    second_pool = {}
    for descriptor in first_descript:
        # Pool values
        all_values = first_descript[descriptor] + second_descript[descriptor]
        # Split randomly in 2
        random.shuffle(all_values)
        split = len(all_values) // 2
        first_pool[descriptor] = all_values[:split]
        second_pool[descriptor] = all_values[split:]
    return first_pool, second_pool


def compare_medians(signif_dict, obs_diff_dict, pool_diff_dict):
    for descriptor in obs_diff_dict:
        if descriptor not in signif_dict:
            signif_dict[descriptor] = []
        obs_diff = obs_diff_dict[descriptor]
        pool_diff = pool_diff_dict[descriptor]
        if np.abs(obs_diff) > np.abs(pool_diff):
            signif_dict[descriptor].append(1)
        else:
            signif_dict[descriptor].append(0)
    return signif_dict


def get_pvals(signif_dict, val1, val2):
    results = []
    for descriptor in signif_dict:
        pval = (np.array(signif_dict[descriptor]) == 0).sum() / len(signif_dict[descriptor])
        results.append([val1, val2, descriptor, pval])
    return results


def export_p_val(denovo_trg_signif, denovo_cds_signif, trg_cds_signif, file_name):
    results = (get_pvals(denovo_trg_signif, "denovo", "trg"))
    results += get_pvals(denovo_cds_signif, "denovo", "cds")
    results += get_pvals(trg_cds_signif, "trg", "cds")
    df = pd.DataFrame(results, columns=["descriptor", "type_1", "type_2", "p_value"])
    df.to_csv(file_name, sep="\t", index=False)


def process_sample(k, n_denovo, trg_names, cds_names, denovo_descript, descriptors, all_cdss):
    """Function to process a single sample iteration."""
    # Sample TRGs and CDSs
    trg_names_sampled = random.sample(list(trg_names), n_denovo)
    cds_names_sampled = random.sample(list(cds_names), n_denovo)

    # Calculate descriptors for all sampled cdss and trgs
    descriptors, trg_descript = calculate_descriptors(descriptors, all_cdss, trg_names_sampled)
    descriptors, cds_descript = calculate_descriptors(descriptors, all_cdss, cds_names_sampled)

    # Get the median differences
    denovo_trg_diff = get_median_diff(denovo_descript, trg_descript)
    denovo_cds_diff = get_median_diff(denovo_descript, cds_descript)
    trg_cds_diff = get_median_diff(trg_descript, cds_descript)

    # Pool sampled cdss and get median diff
    denovo_trg_pool_1, denovo_trg_pool_2 = pool_cdss(denovo_descript, trg_descript)
    denovo_cds_pool_1, denovo_cds_pool_2 = pool_cdss(denovo_descript, cds_descript)
    trg_cds_pool_1, trg_cds_pool_2 = pool_cdss(trg_descript, cds_descript)
    denovo_trg_diff_pool = get_median_diff(denovo_trg_pool_1, denovo_trg_pool_2)
    denovo_cds_diff_pool = get_median_diff(denovo_cds_pool_1, denovo_cds_pool_2)
    trg_cds_diff_pool = get_median_diff(trg_cds_pool_1, trg_cds_pool_2)

    # Compare medians
    denovo_trg_signif = compare_medians({}, denovo_trg_diff, denovo_trg_diff_pool)
    denovo_cds_signif = compare_medians({}, denovo_cds_diff, denovo_cds_diff_pool)
    trg_cds_signif = compare_medians({}, trg_cds_diff, trg_cds_diff_pool)

    return denovo_trg_signif, denovo_cds_signif, trg_cds_signif

















if __name__ == "__main__":
    # Get the list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]
    denovo_names, trg_names, cds_names = [], [], []
    all_cdss = {}
    descriptors = {}
    denovo_trg_signif, denovo_cds_signif, trg_cds_signif = {}, {}, {}


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
    cds_names = set(cds_names) - set(trg_names)
    trg_names = set(trg_names) - set(denovo_names)

    # Calculate descriptors for the denovo genes
    descriptors, denovo_descript = calculate_descriptors(descriptors, all_cdss, denovo_names)

    # Repeat the process n times
    n_denovo = len(denovo_names)
    n = 1
    num_workers = 8 if n < 100000 else 32

    print(f"Starting parallel processing with {num_workers} workers...")
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Use partial to pass fixed arguments to the function
        process_func = partial(process_sample, n_denovo=n_denovo, trg_names=trg_names, cds_names=cds_names,
                               denovo_descript=denovo_descript, descriptors=descriptors, all_cdss=all_cdss)

        # Submit tasks to the executor
        results = list(executor.map(process_func, range(n)))

    # Aggregate results
    for denovo_trg_signif_part, denovo_cds_signif_part, trg_cds_signif_part in results:
        for descriptor in denovo_trg_signif_part:
            if descriptor not in denovo_trg_signif:
                denovo_trg_signif[descriptor] = []
            denovo_trg_signif[descriptor].extend(denovo_trg_signif_part[descriptor])

        for descriptor in denovo_cds_signif_part:
            if descriptor not in denovo_cds_signif:
                denovo_cds_signif[descriptor] = []
            denovo_cds_signif[descriptor].extend(denovo_cds_signif_part[descriptor])

        for descriptor in trg_cds_signif_part:
            if descriptor not in trg_cds_signif:
                trg_cds_signif[descriptor] = []
            trg_cds_signif[descriptor].extend(trg_cds_signif_part[descriptor])


    # Calculate and export the p values
    print("\nDone!")
    export_p_val(denovo_trg_signif, denovo_cds_signif, trg_cds_signif, OUT_DIR + "pvalues.tsv")