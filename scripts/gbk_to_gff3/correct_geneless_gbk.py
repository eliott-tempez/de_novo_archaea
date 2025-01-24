import os

INPUT_FOLD = "/home/eliott.tempez/Documents/archaea_data/complete_122/annotation/"
ORGANISMS = ["Thermococcus_33N", "Thermococcus_celer_SH1", "Thermococcus_sp_690"]










# Iterate on the 3 files
for arch in ORGANISMS:
    for file in os.listdir(INPUT_FOLD):
        if arch not in file:
            break
        with open(INPUT_FOLD + file, "r") as f_in:
            chunk = ""
            # Open the output file (with the gene before each CDS)
            with open(INPUT_FOLD + ".".join(file.split(".")[:-1]) + "_correction.gbk", "w") as f_out:

                for line in f_in:
                    # Get the relevant info
                    # gene range
                    if line.startswith("     CDS"):
                        if chunk:
                            f_out.write(chunk)
                            chunk = ""
                        range = line.split(" ")[-1]
                    # locus tag
                    elif "/locus_tag=" in line:
                        locus_tag = line.split(" ")[-1]
                    # gene name
                    elif "/gene=" in line:
                        gene = line.split(" ")[-1]

                    chunk += line


                    



