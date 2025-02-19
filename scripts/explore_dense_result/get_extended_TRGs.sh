#!/bin/bash
#### This script aims to extract all extended TRG faa files from the tmp files created by dense ####

# Get all file names
FASTA_FILES=$(ls /store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/complete_122/fasta_renamed/*.fa)
cd /datas/ELIOTT/scripts/work/
mkdir extended_CDS_list

for FILE in $FASTA_FILES
do
    # Extract the basename without the extension
    BASENAME=$(basename $FILE .fa)
    # File name to look for
    TARGET_FILE="TRG_multielongated_blastp_$(BASENAME)_CDS_elongated.out"
    # Run ls and get folder name for first match
    TARGET_FILEPATH=$(ls */*/$TARGET_FILE | head -n 1)
    TARGET_FOLDER=$(dirname $TARGET_FILEPATH)
    # Copy the TRH file to the extended_CDS_list folder
    cp $TARGET_FOLDER/TRG_multielongated.faa extended_CDS_list/TRG_multielongated_$(BASENAME).faa
done
