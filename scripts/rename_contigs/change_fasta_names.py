import os
import pandas as pd
import Bio.SeqIO


INPUT_FOLDER_GBK = "/home/eliott.tempez/Documents/archaea_data/complete_122/annotation/"
INPUT_FOLDER_FASTA = "/home/eliott.tempez/Documents/archaea_data/complete_122/fasta/"
OUTPUT_FOLDER_FASTA = "/home/eliott.tempez/Documents/archaea_data/complete_122/fasta_renamed/"
CORR_MAT = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/rename_contigs/correspondance_table.txt"


# Get the basename of a file
def get_basename(filename):
    return ".".join(filename.split(".")[:-1])


# Get the filenames
files_gbk = []
for file in os.listdir(INPUT_FOLDER_GBK):
    if file.endswith(".gbk") or file.endswith(".gbff"):
        files_gbk.append(file)



### Change the name of the genomes with 1 contig ###
for fgbk in files_gbk:
    name = get_basename(fgbk)
    ffasta = INPUT_FOLDER_FASTA + name + ".fa"
    with open(os.path.join(INPUT_FOLDER_GBK, fgbk), "r") as f_in:
        records = list(Bio.SeqIO.parse(f_in, "genbank"))
    # only if we have 1 record
    if len(records) == 1:
        contig_name = records[0].id
    

        # Create new fasta file with accurate name
        with open(ffasta, "r") as f_in:
            records = list(Bio.SeqIO.parse(f_in, "fasta"))
            new_file_content = (f">{contig_name}\n{records[0].seq}\n")
            with open(os.path.join(OUTPUT_FOLDER_FASTA, f"{name}.fa"), "w") as f_out:
                f_out.write(new_file_content)



### Change the name of the genomes with multiple contigs ###
# Populate dict with correspondance table
dict_contigs = {}
with open(CORR_MAT, "r") as f:
    for line in f:
        if len(line.split()) == 1:
            name = line.split()[0]
            dict_contigs[name] = {}
        elif len(line.split()) == 2:
            gbk_contig_name = line.split()[0]
            fasta_contig_name = line.split()[1]
            dict_contigs[name][fasta_contig_name] = gbk_contig_name
# Delete organisms for which we don't have all contigs
keys_to_delete = []
for name in dict_contigs:
    if "0" in dict_contigs[name].keys():
        keys_to_delete.append(name)
        continue
for key in keys_to_delete:
    del dict_contigs[key]
print(dict_contigs)
    

# Create new fasta
for name in dict_contigs:
    old_fasta = INPUT_FOLDER_FASTA + name + ".fa"
    new_fasta = OUTPUT_FOLDER_FASTA + name + ".fa"
    with open(old_fasta, "r") as f_in:
        records = list(Bio.SeqIO.parse(f_in, "fasta"))
    new_file_content = ""
    for record in records:
        new_file_content += f">{dict_contigs[name][record.id]}\n{record.seq}\n"
    with open(new_fasta, "w") as f_out:
        f_out.write(new_file_content)


        
        
        






