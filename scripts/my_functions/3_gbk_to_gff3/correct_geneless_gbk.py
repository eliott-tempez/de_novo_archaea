import os

INPUT_FOLD = "/home/eliott.tempez/Documents/archaea_data/complete_122/annotation/"
ORGANISMS = ["Thermococcus_celer_SH1", "Thermococcus_33N", "Thermococcus_sp_690"]

# Get the files
files = os.listdir(INPUT_FOLD)
files_filtered = []
for arch in ORGANISMS:
    is_found = False
    while not is_found:
        for file in files:
            if arch in file and "correction" not in file:
                files_filtered.append(file)
                is_found = True


# Iterate on the 3 files
for file_name in files_filtered:
    with open(INPUT_FOLD + file_name, "r") as f_in:
        chunk = ""
        current_chunk = ""
        gene_range = ""
        locus_tag = ""
        gene = ""
        # Open the output file (with the gene before each CDS)
        output_file = os.path.join(INPUT_FOLD, ".".join(file_name.split(".")[:-1]) + "_correction.gbk")
        with open(output_file, "w") as f_out:
            for line in f_in:
                # Get the relevant info
                # gene range
                if line.startswith("     CDS"):
                    # If we have all info, write the gene line
                    if gene_range != "":
                        if chunk != "":
                            f_out.write(chunk)
                            chunk = ""
                        f_out.write(f"     gene            {gene_range}\n")
                        if gene != "":
                            f_out.write(f"                     {gene}\n")
                        if locus_tag != "":
                            f_out.write(f"                     {locus_tag}\n")
                        gene_range = ""
                        locus_tag = ""
                        gene = ""

                    # If beginning of the file
                    if chunk == "":
                        f_out.write(current_chunk)
                        current_chunk = ""

                    gene_range = line.split()[-1]
                    chunk = current_chunk
                    current_chunk = ""
                    
                # gene name
                elif "/gene=" in line:
                    gene = line.split()[-1]

                # locus tag
                elif "/locus_tag=" in line:
                    locus_tag = line.split()[-1]

                current_chunk += line


            # Write the last gene
            f_out.write(f"     gene            {gene_range}\n")
            if gene != "":
                f_out.write(f"                     {gene}\n")
            if locus_tag != "":
                f_out.write(f"                     {locus_tag}\n")
            f_out.write(current_chunk)

