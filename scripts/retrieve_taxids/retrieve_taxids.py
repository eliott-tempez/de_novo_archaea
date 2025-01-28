import os
import random
import subprocess
from Bio import SeqIO


INPUT_FOLD = "/home/eliott.tempez/Documents/archaea_data/complete_122/annotation/"
OUTPUT_FOLD = "/home/eliott.tempez/Documents/archaea_data/retrieve_taxid/results/"
DIAMOND_EXEC = "/home/eliott.tempez/Programs/diamond/diamond"
NR = "/datas/NR/nr_2.0.13.dmnd"

# Get the files
files = os.listdir(INPUT_FOLD)

f = files[0]
print(f)

# For each species
# Get all CDS
cds_list = []
for gb_record in SeqIO.parse(open(INPUT_FOLD + f, "r"), "genbank"):
    cds_list += [f for f in gb_record.features if f.type == "CDS"]
nb_cds = len(cds_list)

# Sample 100 random CDS
if nb_cds > 100:
    cds_list_100 = random.sample(cds_list, 100)
else:
    cds_list_100 = cds_list

# Export in fasta format
if not os.path.exists(OUTPUT_FOLD + "fasta/"):
    os.makedirs(OUTPUT_FOLD + "fasta/")
fasta_file = f"{OUTPUT_FOLD + "fasta/"}{".".join(f.split('.')[:-1])}.fasta"
with open(fasta_file, "w") as f_out:
    for cds in cds_list_100:
        f_out.write(f">{cds.qualifiers['locus_tag'][0]}\n{cds.extract(gb_record.seq)}\n")
    
# Run diamond in bash
if not os.path.exists(OUTPUT_FOLD + "diamond/"):
    os.makedirs(OUTPUT_FOLD + "diamond/")
diamond_command = f"{DIAMOND_EXEC} -d {NR} -p 32 -q {fasta_file} -o {OUTPUT_FOLD + "diamond/"}{'.'.join(f.split('.')[:-1])}.m8"
subprocess.run(diamond_command, shell=True)


