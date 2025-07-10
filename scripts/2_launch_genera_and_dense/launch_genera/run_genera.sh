#!/bin/sh
# /!\ The paths are probably not correct, you will have to change them

########## Parameters ##########

# Get the species from the argument
archaea=$SPECIES

# Data files
NR=/datas/NR/nr_2.0.13
TAXDUMP=/datas/ELIOTT/scripts/work/taxdump/
TAXID_FILE=/datas/ELIOTT/archaea_data/taxid.csv
# Output dir
OUT_DIR=/datas/ELIOTT/archaea_data/genera/
# Log filenames
LOG_OUTPUT=/home/eliott.tempez/genera_output_$archaea.log
LOG_ERROR=/home/eliott.tempez/genera_error_$archaea.log

# Set up environment
source /home/eliott.tempez/miniconda3/bin/activate phylostrat
mkdir -p /datas/ELIOTT/scripts/${archaea}
cd /datas/ELIOTT/scripts/${archaea}


########## Run GenEra ##########

echo $archaea >> $LOG_OUTPUT

# Get taxid
TAXID=$(grep $archaea $TAXID_FILE | cut -f2 -d$'\t')
echo "Taxid: $TAXID" >> $LOG_OUTPUT
# Get filenames
FASTA_CDS=/datas/ELIOTT/archaea_data/reannotated_CDS/${archaea}_CDS.faa
GFF_FILE=/datas/ELIOTT/archaea_data/reannotated_gff_75/${archaea}.gff3
# Create temp taxid file with all other fasta paths
echo "Creating temp taxid file..." >> $LOG_OUTPUT
TAXID_FILE_FASTA=/datas/ELIOTT/archaea_data/genera/taxid_tmp.csv
grep -v $archaea $TAXID_FILE | awk -F'\t' '{print "/datas/ELIOTT/archaea_data/reannotated_CDS/" $1 "_CDS.faa \t " $2}' > $TAXID_FILE_FASTA
echo "Temp taxid file created" >> $LOG_OUTPUT
# Create output dir
mkdir -p $OUT_DIR/out/${archaea}

# Run GenEra
echo "Running GenEra..." >> $LOG_OUTPUT
genEra -t $TAXID -q $FASTA_CDS -a $TAXID_FILE_FASTA -n 32 -b $NR -r /datas/ELIOTT/scripts/ncbi_lineages_2025-02-15.csv -o ${OUT_DIR}out/${archaea}/ >> $LOG_OUTPUT 2>> $LOG_ERROR

# Deactivate conda environment
conda deactivate
rm $TAXID_FILE_FASTA
rm -r /datas/ELIOTT/scripts/${archaea}
echo -e "Job completed successfully\n" >> $LOG_OUTPUT
