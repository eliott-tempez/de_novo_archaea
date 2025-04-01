
import os
import glob
import re
import pandas as pd
from my_functions.genomic_functions import extract_cds_info
from Bio.SeqUtils import GC
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np


OUT_DIR = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/13_plot_dense_results/"
from my_functions.paths import DENSE_DIR, GENERA_DIR, GENOMES_LIST
TRG_RANK = 7.0


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


def extract_denovo_names(focal_species):
    """Extract the names of the de novo genes"""
    # Make sure the file exists and read it
    denovo_file = os.path.join(DENSE_DIR, focal_species, "denovogenes.tsv")
    if os.path.exists(denovo_file):
        denovo_df = pd.read_csv(denovo_file, sep="\t", header=0)
    else:
        raise FileNotFoundError(f"No file {denovo_file}")
    # Extract the names
    denovo_names = denovo_df["CDS"].tolist()
    return denovo_names






if __name__ == "__main__":
    # Get the list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]
    all_values = []

    i = 0
    for genome in genomes:
        i += 1
        # Extract info for all CDSs
        cds_dict_tmp = extract_cds_info(genome)
        cds_dict = {}
        #Rename
        for cds in cds_dict_tmp:
            cds_dict[f"{cds}_gene_mRNA"] = cds_dict_tmp[cds]
        cds_names = list(cds_dict.keys())
        # Extract TRG names
        trg_names = extract_trg_names(genome, TRG_RANK)
        # Extract de novo names
        denovo_names = extract_denovo_names(genome)

        # Extract GC rate
        gc_dict = {}
        aa_percent_dict = {}
        aromaticity_dict = {}
        instability_dict = {}
        flexibility_dict = {}
        for cds in cds_dict:
            # Extract GC rate
            gc_dict[cds] = GC(cds_dict[cds]["seq"])
            # Extract protein info
            prot_seq = re.sub(r"[\*X]", "", str(cds_dict[cds]["seq"].translate()))
            analysis = ProteinAnalysis(prot_seq)
            aa_percent_dict[cds] = analysis.get_amino_acids_percent()
            aromaticity_dict[cds] = analysis.aromaticity()
            instability_dict[cds] = analysis.instability_index()
            flexibility_dict[cds] = analysis.flexibility()

        # Get the mean for all the values
        ## GC
        gc_mean_cds = sum(list(gc_dict.values())) / len(gc_dict)
        gc_mean_trg = sum([v for k, v in gc_dict.items() if k in trg_names]) / len(trg_names)
        gc_mean_denovo = sum([v for k, v in gc_dict.items() if k in denovo_names]) / len(denovo_names) if denovo_names else np.nan
        all_values.append([genome, "gc", gc_mean_cds, gc_mean_trg, gc_mean_denovo])

        ## Aromaticity
        aromaticity_mean_cds = sum(list(aromaticity_dict.values())) / len(aromaticity_dict)
        aromaticity_mean_trg = sum([v for k, v in aromaticity_dict.items() if k in trg_names]) / len(trg_names)
        aromaticity_mean_denovo = sum([v for k, v in aromaticity_dict.items() if k in denovo_names]) / len(denovo_names) if denovo_names else np.nan
        all_values.append([genome, "aromaticity", aromaticity_mean_cds, aromaticity_mean_trg, aromaticity_mean_denovo])

        ## Instability
        instability_mean_cds = sum(list(instability_dict.values())) / len(instability_dict)
        instability_mean_trg = sum([v for k, v in instability_dict.items() if k in trg_names]) / len(trg_names)
        instability_mean_denovo = sum([v for k, v in instability_dict.items() if k in denovo_names]) / len(denovo_names) if denovo_names else np.nan
        all_values.append([genome, "instability", instability_mean_cds, instability_mean_trg, instability_mean_denovo])

        if i % 10 == 0:
            print(f"{i}/{len(genomes)}...")
    print("Done!")

    # Save the results
    df = pd.DataFrame(all_values, columns=["genome", "feature", "cds", "trg", "denovo"])
    df.to_csv(os.path.join(OUT_DIR, "sequence_features.csv"), sep="\t", index=False)