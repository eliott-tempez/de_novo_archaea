#!/bin/bash
# This scripts removes the fasta part of all the gff files in a directory
# Using the script https://github.com/jorvis/biocode/blob/master/gff/remove_fasta_from_gff3.py

INPUT_FOLD=/home/eliott.tempez/Documents/archaea_data/complete_122/gff3/
OUTPUT_FOLD=/home/eliott.tempez/Documents/archaea_data/complete_122/gff3_no_fasta/
mkdir -p $OUTPUT_FOLD

for file in $INPUT_FOLD/* ; do
    # Check if the file is a directory
    if [ -d "$file" ]; then
        continue
    fi
    # Get the base name of the file
    base_name=$(basename "$file"  | sed 's/\.[^.]*$//')
    # Create the output file
    output_file="${OUTPUT_FOLD}/${base_name}.gff3"
    # Run the script
    python3 /home/eliott.tempez/Documents/archaea_data/gbk_to_gff/remove_fasta_from_gff3.py -i "$file" -o "$output_file"
done