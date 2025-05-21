#!/bin/bash
# This script converts all the genbank files in the current directory to gff3 files
# Using the script https://github.com/jorvis/biocode/blob/master/gff/convert_genbank_to_gff3.py

log_script="/home/eliott.tempez/Documents/M2_Stage_I2BC/logs/gbk_to_gff3/gbk_to_gff_log_2"
err_script="/home/eliott.tempez/Documents/M2_Stage_I2BC/logs/gbk_to_gff3/gbk_to_gff_err_2"

cd /home/eliott.tempez/Documents/archaea_data/complete_122/annotation

for file in ./*correction.gbk; do
    base_name=$(basename "$file" | sed 's/\.[^.]*$//')
    output_file="../gff3/${base_name}.gff3"
    echo "${base_name}:" >> "$log_script"
    echo "${base_name}:" >> "$err_script"
    python /home/eliott.tempez/Documents/archaea_data/gbk_to_gff/convert_genbank_to_gff3.py -i "$file" -o "$output_file" >> "$log_script" 2>> "$err_script"
done
