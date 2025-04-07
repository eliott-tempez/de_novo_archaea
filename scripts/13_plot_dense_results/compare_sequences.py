
import os
import glob
import re
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np



OUT_DIR = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/13_plot_dense_results/sequences/"
from my_functions.paths import DENSE_DIR, GENERA_DIR, GENOMES_LIST, CDS_DIR
TRG_RANK = 7.0


def extract_cds_sequences(genome):
    cds_dict = {}
    faa_file = os.path.join(CDS_DIR, f"{genome}_CDS.fna")
    if not os.path.exists(faa_file):
        raise FileNotFoundError(f"No file {faa_file}")
    for record in SeqIO.parse(faa_file, "fasta"):
        cds_dict[record.name] = record.seq
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
        denovo_file = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/good_candidates.txt"
        with open(denovo_file, "r") as f:
            for line in f:
                denovo_names.append(line.strip())
    return denovo_names



if __name__ == "__main__":
    # Get the list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]
    all_values = []
    n_cds = 0
    n_trg = 0
    n_denovo = 0

    i = 0
    for genome in genomes:
        i += 1
        # Extract info for all CDSs
        cds_dict = extract_cds_sequences(genome)
        cds_names = list(cds_dict.keys())
        n_cds += len(cds_names)
        # Extract TRG names
        trg_names = extract_trg_names(genome, TRG_RANK)
        n_trg += len(trg_names)
        # Extract de novo names
        denovo_names = extract_denovo_names(genome, True)
        n_denovo += len(denovo_names)


        for cds in cds_dict:
            # Extract GC rate
            gc_content = GC(cds_dict[cds])
            # Extract protein info
            nuc_seq = cds_dict[cds]
            prot_seq = re.sub(r"[\*X]", "", str(cds_dict[cds].translate()))
            analysis = ProteinAnalysis(prot_seq)
            aromaticity = analysis.aromaticity()
            instability = analysis.instability_index()
            flexibility = analysis.flexibility()
            mean_flexibility = sum(flexibility) / len(flexibility)
            hydropathy = analysis.gravy()
            # Extract sequence length
            len_nu = len(nuc_seq)

            # Add to results
            results = [genome, cds, gc_content, aromaticity, instability, mean_flexibility, hydropathy, len_nu]
            all_values.append(results + ["cds"])
            if cds in trg_names:
                all_values.append(results + ["trg"])
            if cds in denovo_names:
                all_values.append(results + ["denovo"])

        if i % 10 == 0:
            print(f"{i}/{len(genomes)}...")

    print("Done!")

    # Save the results
    df = pd.DataFrame(all_values, columns=["genome", "cds", "gc_content", "aromaticity", "instability", "mean_flexibility", "hydropathy", "len_nu", "type"])
    df.to_csv(os.path.join(OUT_DIR, "sequence_features_good_candidates.csv"), sep="\t", index=False)