"""
The data of the 116 Thermococcacae did not contain associated taxids, which are necessary for Dense.
This scripts uses Diamond to Blast 100 random CDS of each genome against the NR database.
It returns a m8 files with the results.

Author: Eliott Tempez.
M2 internship - 2025.
"""


import os
import re
import random
import subprocess
from Bio import SeqIO

from my_functions.paths import GENOMES_LIST, CDS_DIR

OUTPUT_FOLD = "/datas/ELIOTT/retrieve_taxid/"
NR = "/datas/NR/nr_2.0.13.dmnd"

# Get name of genomes
with open(GENOMES_LIST, "r") as f:
    genomes = f.readline().split()
genomes = [re.sub('"', '', g) for g in genomes]

# Read all CDSs for each genome
cdss_dict = {}
for genome in genomes:
    cds_file = f"{CDS_DIR}{genome}.faa"
    for cds in SeqIO.parse(cds_file, "fasta"):
        seq = str(cds.seq)
        cdss_dict[cds.name] = seq

    # Sample 100 random CDS
    sampled_cdss = random.sample(list(cdss_dict.keys()), 100)

    # Export in fasta format
    if not os.path.exists(OUTPUT_FOLD + "fasta/"):
        os.makedirs(OUTPUT_FOLD + "fasta/")
    fasta_file = f"{OUTPUT_FOLD}fasta/{genome}.fasta"
    with open(fasta_file, "w") as f_out:
        for cds in sampled_cdss:
            f_out.write(f">{genome}|{cds}\n{cdss_dict[cds]}\n")
        
    # Run diamond in bash
    if not os.path.exists(OUTPUT_FOLD + "diamond/"):
        os.makedirs(OUTPUT_FOLD + "diamond/")
    output_format = "6 qseqid sseqid staxids pident pident bitscore"
    diamond_command = f"diamond blastp -d {NR} -p 32 -q {fasta_file} -o {OUTPUT_FOLD}diamond/{genome}.m8 --taxonlist 2157 -k 1 -f {output_format} --very-sensitive --quiet" 
    subprocess.run(diamond_command, shell=True)


