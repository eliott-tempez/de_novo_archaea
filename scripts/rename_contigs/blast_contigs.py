import os
import subprocess
import Bio.SeqIO


INPUT_FOLDER_GBK = "/home/eliott.tempez/Documents/archaea_data/complete_122/annotation/"
INPUT_FOLDER_FASTA = "/home/eliott.tempez/Documents/archaea_data/complete_122/fasta/"
OUTPUT_FILE = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/rename_contigs/blast_result.txt"


# Get the basename of a file
def get_basename(filename):
    return ".".join(filename.split(".")[:-1])


# Get the filenames
files_gbk = []
files_fasta = []
for file in os.listdir(INPUT_FOLDER_GBK):
    if file.endswith(".gbk") or file.endswith(".gbff"):
        files_gbk.append(file)
for file in os.listdir(INPUT_FOLDER_FASTA):
    if file.endswith(".fa"):
        files_fasta.append(file)


# Blast each protein in the gbk file against the fasta file
for fgbk in files_gbk:
    name = get_basename(fgbk)
    ffasta = INPUT_FOLDER_FASTA + name + ".fa"
    with open(os.path.join(INPUT_FOLDER_GBK, fgbk), "r") as f_in:
        records = list(Bio.SeqIO.parse(f_in, "genbank"))
        # only if we have more than 1 record
        if len(records) > 1:
            print(name)
            for record in records:
                print(record.id)
                for feature in record.features:
                    if feature.type == "CDS":
                        if "translation" in feature.qualifiers:
                            tmp_file_content = (f">{record.id}\n{feature.qualifiers['translation'][0]}\n")
                            with open("tmp.fa", "w") as f_out:
                                f_out.write(tmp_file_content)
                            # Create the blast command
                            blast_command = f"tblastn -num_threads 10 -query tmp.fa -subject {ffasta} -outfmt '6 qseqid sseqid evalue length' | sort -k3,3g | head -n 5"
                            # Run the blast command and capture the output
                            result = subprocess.run(blast_command, shell=True, capture_output=True, text=True)
                            result_list = result.stdout.split("\n")
                            result_updated = [name + "\t" + l for l in result_list if l]
                            if result_updated:
                                with open(OUTPUT_FILE, "a") as f_out:
                                    for line in result_updated:
                                        f_out.write((line) + "\n")
                            os.remove("tmp.fa")
            print("\n")
