#!/bin/sh

########## Parameters ##########

#PBS -N genera_archaea
#PBS -q bim
#PBS -l ncpus=32 -l host=node04 -l mem=128gb -l walltime=200:00:00
#PBS -o /home/eliott.tempez/genera_output.log
#PBS -e /home/eliott.tempez/genera_error.log

# Log filenames
LOG_OUTPUT=/home/eliott.tempez/genera_output.log
LOG_ERROR=/home/eliott.tempez/genera_error.log

# Data files
NR=/datas/NR/nr_2.0.13.dmnd
TAXDUMP=/datas/ELIOTT/work/taxdump
TAXID_FILE=/datas/ELIOTT/archaea_data/taxid.csv
TAXID_FILE_FASTA=/datas/ELIOTT/archaea_data/taxid_fasta.csv
TREE=/datas/ELIOTT/archaea_data/whole_tree.nwk

# Output dir
OUT_DIR=/datas/ELIOTT/archaea_data/genera/

# Set up environment
source /home/eliott.tempez/miniconda3/bin/activate phylostrat
cd /datas/ELIOTT/scripts/



########## Run GenEra ##########
ARCHAEA=GCA_000195935@Pyrococcus_abyssi_GE5


echo $ARCHAEA >> $LOG_OUTPUT
# Get taxid
TAXID=$(grep $ARCHAEA $TAXID_FILE | cut -d, -f2)
# Get filenames
FASTA_CDS=/datas/ELIOTT/archaea_data/genera/CDS/${ARCHAEA}_CDS.fna
GFF_FILE=/datas/ELIOTT/archaea_data/genome/${ARCHAEA}.gff3


echo "Running GenEra" >> $LOG_OUTPUT
genEra -t $TAXID -q $FASTA_CDS -a $TAXID_FILE_FASTA -n 32 -b $NR -d $TAXDUMP -o {$OUT_DIR}out/ >> $LOG_OUTPUT 2>> $LOG_ERROR

# Deactivate conda environment
conda deactivate
echo "Job completed successfully" >> $LOG_OUTPUT