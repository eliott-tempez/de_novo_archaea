"""
This script extracts all intergenic open reading frames (iorfs) from a list of
genomes using the orftrack and orfget tools. It saves the iorfs in a text file
iorfs.txt.
"""

import re
import os
import subprocess


from my_functions.paths import GENOMES_LIST, FA_DIR, GFF_DIR
OUT_PATH = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/7_characterise_denovo"



def extract_iorfs(species):
    fa_file = os.path.join(FA_DIR, species + ".fa")
    gff_file = os.path.join(GFF_DIR, species + ".gff3")
    orftrack_file = f"mapping_orf_{species}.gff"
    orfget_file = f"mapping_orf_{species}_nc_intergenic.nfasta"
    log_file = "orftrack.log"

    iorfs = []
    # Use orftrack to annotate the genome
    result_orftrack = subprocess.run(
        ["orftrack", "--fna", fa_file, "--gff", gff_file],
        capture_output=True,
        text=True
    )
    if result_orftrack.returncode != 0:
        print(f"Error in orftrack: {result_orftrack.stderr}")

    result_orfget = subprocess.run(
        ["orfget", "--fna", fa_file, "--gff", orftrack_file, "--features_include", "nc_intergenic", "--type", "nucl"],
        capture_output=True,
        text=True
    )
    if result_orfget.returncode != 0:
        print(f"Error in orfget: {result_orfget.stderr}")

    # Read the iorfs
    with open(orfget_file, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq = line.strip()
                iorfs.append(seq)
    
    # Remove the temp files
    os.remove(orftrack_file)
    os.remove(orfget_file)
    os.remove(log_file)

    # Return the iorfs
    return iorfs




if __name__ == "__main__":
    # Get the list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]

    out_file = os.path.join(OUT_PATH, "iorfs.txt")
    # Get the iorfs for each genome
    for genome in genomes:
        iorfs = extract_iorfs(genome)
        # Write the iorfs to the end of the file
        with open(out_file, "a") as f:
            f.write(f">{genome}\n")
            f.write(" ".join(iorfs) + "\n")
        print(f"Extracted {len(iorfs)} iorfs for {genome}")