import os
import glob
import re
from Bio import SeqIO
from Bio.SeqUtils import GC

DATA_DIR = "/home/eliott.tempez/Documents/archaea_data/complete_122/"
GENOMES_LIST = "/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/6_genera_archaea/genomes_list.txt"
OUTPUT_FILE = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/10_analyse_intergenic/intergenic_lengths.tsv"
OUTPUT_FILE_GC = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/10_analyse_intergenic/intergenic_gc.tsv"




def extract_intergenic(species):
        """Extract all intergenic sequences"""
        # Make sure the file exists and read it
        gbk_pattern = os.path.join(DATA_DIR, "reannotated_gbk_75", species + ".gbk")
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
                if feature.type == "CDS":
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
    

# Read the list of genomes
with open(GENOMES_LIST, "r") as f:
    genomes = f.readline().split()
genomes = [re.sub('"', '', g) for g in genomes]


mean_intergenic_lengths = {}
intergenic_gc = {}
for genome in genomes:
    all_intergenic_seqs = ""
    intergenic_dict = extract_intergenic(genome)
    # Get the mean length of the intergenic sequences
    intergenic_lengths = [len(seq) for seq in intergenic_dict.values()]
    mean_intergenic_length = sum(intergenic_lengths) / len(intergenic_lengths)
    mean_intergenic_lengths[genome] = mean_intergenic_length

    # Concatenate all intergenic seqs
    for seq in intergenic_dict.values():
        all_intergenic_seqs += str(seq)
    # Calculate GC content
    intergenic_gc[genome] = GC(all_intergenic_seqs)

# Write the results to a file
with open(OUTPUT_FILE, "w") as f:
    f.write("genome\tmean_intergenic_length\n")
    for genome, mean_intergenic_length in mean_intergenic_lengths.items():
        f.write(f"{genome}\t{mean_intergenic_length}\n")

with open(OUTPUT_FILE_GC, "w") as f:
    f.write("genome\tintergenic_gc\n")
    for genome, gc_content in intergenic_gc.items():
        f.write(f"{genome}\t{gc_content}\n")
    