import os
import random
import subprocess
import glob
from Bio import SeqIO

TAXID = "/datas/ELIOTT/archaea_data/taxid.tsv"
INPUT_FOLD = "/datas/ELIOTT/archaea_data/genome/gbk/"
OUTPUT_FOLD = "/datas/ELIOTT/archaea_data/retrieve_taxid/"
NR = "/datas/NR/nr_2.0.13.dmnd"

# Get name of unknown organisms
organisms = []
with open(TAXID, "r") as f_taxid:
    for line in f_taxid:
        org = line.split("\t")[0]
        status = line.split("\t")[2]
        if int(status) == 0:
            organisms.append(org)

# For each species
for name in organisms:
    print(name, ":")
    input_file_path = glob.glob(INPUT_FOLD + name + ".*")[0]
    # Get all CDS with a protein sequence
    cds_list = []
    with open(input_file_path, "r") as input_handle:
        for gb_record in SeqIO.parse(input_handle, "genbank"):
            cds_list += [f for f in gb_record.features if f.type == "CDS" and "translation" in f.qualifiers]
    nb_cds = len(cds_list)

    # Sample 100 random CDS
    print("Retrieving 100 random CDS...")
    if nb_cds > 100:
        cds_list_100 = random.sample(cds_list, 100)
    else:
        cds_list_100 = cds_list

    # Export in fasta format
    print("Exporting in fasta format...")
    if not os.path.exists(OUTPUT_FOLD + "fasta/"):
        os.makedirs(OUTPUT_FOLD + "fasta/")
    fasta_file = f"{OUTPUT_FOLD}fasta/{name}.fasta"
    with open(fasta_file, "w") as f_out:
        for cds in cds_list_100:
            f_out.write(f">{name}|{cds.qualifiers['locus_tag'][0]}\n{cds.qualifiers['translation'][0]}\n")
        
    # Run diamond in bash
    print("Running diamond...")
    if not os.path.exists(OUTPUT_FOLD + "diamond/"):
        os.makedirs(OUTPUT_FOLD + "diamond/")
    output_format = "6 qseqid sseqid staxids pident pident bitscore"
    diamond_command = f"diamond blastp -d {NR} -p 32 -q {fasta_file} -o {OUTPUT_FOLD}diamond/{name}.m8 --taxonlist 2157 -k 1 -f {output_format} --very-sensitive --quiet" 
    print(f"{OUTPUT_FOLD}diamond/{name}.m8")
    subprocess.run(diamond_command, shell=True)

    print("done!\n")


