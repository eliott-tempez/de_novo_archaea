import os

INPUT_FOLD = "/home/eliott.tempez/Documents/archaea_data/complete_122/annotation/"
ORGANISMS = ["Thermococcus_33N", "Thermococcus_celer_SH1", "Thermococcus_sp_690"]

# Get the files
files = os.listdir(INPUT_FOLD)
files_filtered = []
is_found = False
for arch in ORGANISMS:
    is_found = False
    while not is_found:
        for file in files:
            if arch in file and "correction" not in file:
                files_filtered.append(file)
                is_found = True

has_seen_genes = False
# Iterate on the 3 files
for f in files_filtered:
    with open(INPUT_FOLD + f, "r") as f_in:
        chunk = ""
        current_chunk = ""
        range = ""
        locus_tag = ""
        gene = ""
        # Open the output file (with the gene before each CDS)
        with open(INPUT_FOLD + ".".join(f.split(".")[:-1]) + "_correction.gbk", "w") as f_out:
            for line in f_in:

                # Get the relevant info
                # gene range
                if line.startswith("     CDS"):
                    # If we have all info, write the gene line
                    if range != "":
                        if chunk != "":
                            f_out.write(chunk)
                            chunk = ""
                        f_out.write(f"     gene            {range}\n")
                        if gene != "":
                            f_out.write(f"                     {gene}\n")
                        if locus_tag != "":
                            f_out.write(f"                     {locus_tag}\n")
                        range = ""
                        locus_tag = ""
                        gene = ""

                    # If beginning of the file
                    if chunk == "":
                        f_out.write(current_chunk)
                        current_chunk = ""

                    range = line.split()[-1]
                    chunk = current_chunk
                    current_chunk = ""
                    
                # gene name
                elif "/gene=" in line:
                    gene = line.split()[-1]

                # locus tag
                elif "/locus_tag=" in line:
                    locus_tag = line.split()[-1]

                current_chunk += line
            f_out.write(current_chunk)

