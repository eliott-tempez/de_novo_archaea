import os
import Bio.SeqIO

GBK_FOLDER = "/home/eliott.tempez/Documents/archaea_data/complete_122/annotation/"
FASTA_FOLDER = "/home/eliott.tempez/Documents/archaea_data/complete_122/fasta/"
OUT_FILE = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/rename_contigs/contigs_names.txt"

# Get the files
def process_files(input_folder, file_type):
    dict = {}
    # Get the files
    files = os.listdir(input_folder)
    for file in files:
        # Check we have a file
        if not os.path.isfile(os.path.join(input_folder, file)):
            continue
        # Check we have the right files
        if "Thermococcus_celer_SH1" in file or "Thermococcus_sp_690" in file or "correction" in file:
            continue
        file_path = os.path.join(input_folder, file)
        species_name = ".".join(file.split(".")[:-1])
        dict[species_name] = []
        # read with seqIO
        with open(file_path, "r") as f_in:
            records = Bio.SeqIO.parse(f_in, file_type)
            # Add contig name to the dict
            for record in records:
                dict[species_name].append(record.name)
    return dict



gbk_dict = process_files(GBK_FOLDER, "genbank")
fasta_dict = process_files(FASTA_FOLDER, "fasta")

# Print the names in a file
with open(OUT_FILE, "w") as f_out:
    for species in gbk_dict:
        f_out.write(f"{species}\n")
        f_out.write(f"{gbk_dict[species]}\n")
        f_out.write(f"{fasta_dict[species]}\n")
        f_out.write("\n")
