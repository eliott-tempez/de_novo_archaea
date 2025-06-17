import re
import os


from my_functions.paths import DENSE_DIR, GENOMES_LIST, CDS_DIR, GFF_DIR, FA_DIR
from my_functions.genomic_functions import extract_denovo_info

good_fasta = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/18_clustering/good_candidates.fasta"
all_denovo_file = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/all_denovo.txt"
good_denovo_file = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/good_candidates.txt"
out_folder = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/my_functions"





# Read the denovo names
all_denovo, good_denovo, bad_denovo = [], [], []
with open(all_denovo_file, "r") as f:
    for line in f:
        all_denovo.append(line.strip())
with open(good_denovo_file, "r") as f:
    for line in f:
        good_denovo.append(line.strip())
bad_denovo = list(set(all_denovo) - set(good_denovo))

# Get the list of genomes
with open(GENOMES_LIST, "r") as f:
    genomes = f.readline().split()
genomes = [re.sub('"', '', g) for g in genomes]
denovo_names = []

# Extract all denovo info
all_denovo_dict = {}
for genome in genomes:
    local_dict = extract_denovo_info(genome)
    for denovo in local_dict:
        local_dict[denovo]["genome"] = genome
    all_denovo_dict.update(local_dict)

# Write all denovo info to file
# Chosen denovo
with open(os.path.join(out_folder, "chosen_denovo.tsv"), "w") as f:
    f.write("denovo\tgenome\tcontig\tstart\tend\tstrand\n")
    for denovo in good_denovo:
        genome = all_denovo_dict[denovo]["genome"]
        loci = all_denovo_dict[denovo]["denovo_loci"]
        contig = loci[0]
        start = loci[1]
        end = loci[2]
        strand = loci[3]
        f.write(f"{denovo}\t{genome}\t{contig}\t{start}\t{end}\t{strand}\n")
# Discarded denovo
with open(os.path.join(out_folder, "discarded_denovo.tsv"), "w") as f:
    f.write("denovo\tgenome\tcontig\tstart\tend\tstrand\n")
    for denovo in bad_denovo:
        genome = all_denovo_dict[denovo]["genome"]
        loci = all_denovo_dict[denovo]["denovo_loci"]
        contig = loci[0]
        start = loci[1]
        end = loci[2]
        strand = loci[3]
        f.write(f"{denovo}\t{genome}\t{contig}\t{start}\t{end}\t{strand}\n")


# Write fasta file
with open(good_fasta, "w") as f:
    for denovo in good_denovo:
        seq = all_denovo_dict[denovo]["sequence"]
        f.write(f">{denovo}\n{seq}\n")