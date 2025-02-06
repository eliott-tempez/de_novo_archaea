#!/bin/sh

########## Parameters ##########

#PBS -N genera_archaea_1
#PBS -q bim
#PBS -l ncpus=32 -l host=node04 -l mem=128gb -l walltime=1000:00:00
#PBS -o /home/eliott.tempez/genera_output.log
#PBS -e /home/eliott.tempez/genera_error.log

# Log filenames
LOG_OUTPUT=/home/eliott.tempez/genera_output.log
LOG_ERROR=/home/eliott.tempez/genera_error.log

# Data files
NR=/datas/NR/nr_2.0.13
TAXDUMP=/datas/ELIOTT/scripts/taxdump/
TAXID_FILE=/datas/ELIOTT/archaea_data/taxid.csv

# Output dir
OUT_DIR=/datas/ELIOTT/archaea_data/genera/

# Set up environment
source /home/eliott.tempez/miniconda3/bin/activate phylostrat
cd /datas/ELIOTT/scripts/



########## Run GenEra ##########

# List of archaeas to iterate over
declare -a archaeas=("GCA_028471785@Thermococcus_kodakarensis_TS900")

# Run genera for each of them
for archaea in "${archaeas[@]}"; do

    echo $archaea >> $LOG_OUTPUT

    # Get taxid
    TAXID=$(grep $archaea $TAXID_FILE | cut -f2 -d$'\t')
    echo "Taxid: $TAXID" >> $LOG_OUTPUT
    # Get filenames
    FASTA_CDS=/datas/ELIOTT/archaea_data/genera/CDS/${archaea}_CDS.faa
    GFF_FILE=/datas/ELIOTT/archaea_data/genome/${archaea}.gff3
    # Create temp taxid file with all other fasta paths
    echo "Creating temp taxid file..." >> $LOG_OUTPUT
    TAXID_FILE_FASTA=/datas/ELIOTT/archaea_data/genera/taxid_tmp.csv
    grep -v $archaea $TAXID_FILE | awk -F'\t' '{print "/datas/ELIOTT/archaea_data/genera/CDS/" $1 "_CDS.faa \t " $2}' > $TAXID_FILE_FASTA
    echo "Temp taxid file created" >> $LOG_OUTPUT
    # Create output dir
    mkdir -p $OUT_DIR/out/${archaea}

    # Run GenEra
    echo "Running GenEra..." >> $LOG_OUTPUT
    genEra -t $TAXID -q $FASTA_CDS -a $TAXID_FILE_FASTA -n 32 -b $NR -d $TAXDUMP -o ${OUT_DIR}out/${archaea}/ >> $LOG_OUTPUT 2>> $LOG_ERROR

    # Deactivate conda environment
    conda deactivate
    rm $TAXID_FILE_FASTA
    echo -e "Job completed successfully\n" >> $LOG_OUTPUT

done