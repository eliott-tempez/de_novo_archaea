import os
import pandas as pd
import Bio.SeqIO


INPUT_FOLDER_GBK = "/home/eliott.tempez/Documents/archaea_data/complete_122/annotation/"
INPUT_FOLDER_FASTA = "/home/eliott.tempez/Documents/archaea_data/complete_122/fasta/"
OUTPUT_FOLDER_FASTA = "/home/eliott.tempez/Documents/archaea_data/complete_122/fasta_renamed/"


# Get the basename of a file
def get_basename(filename):
    return ".".join(filename.split(".")[:-1])


# Get the filenames
"""files_gbk = []
for file in os.listdir(INPUT_FOLDER_GBK):
    if file.endswith(".gbk") or file.endswith(".gbff"):
        files_gbk.append(file)"""


"""### Change the name of the genomes with 1 contig ###
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
                f_out.write(new_file_content)"""





        
        
        






