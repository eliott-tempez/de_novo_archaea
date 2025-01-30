#!/bin/sh

########## Parameters ##########

#PBS -N genera_archaea_test
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
TAXID_FILE=/datas/ELIOTT/archaea_data/test/taxid.csv
TREE=/datas/ELIOTT/archaea_data/test/small-tree.nwk

# Output dir
OUT_DIR=/datas/ELIOTT/archaea_data/test/genera/

# Set up environment
source /home/eliott.tempez/miniconda3/bin/activate phylostrat
cd /datas/ELIOTT/



########## Run GenEra ##########
ARCHAEA=GCA_000022545@Thermococcus_sibiricus_MM_739.fa


echo $ARCHAEA >> $LOG_OUTPUT
# Get taxid
TAXID=$(grep $ARCHAEA $TAXID_FILE | cut -d, -f2)
# Get filenames
FASTA_CDS=/datas/ELIOTT/archaea_data/test/genera/CDS/${ARCHAEA}.fa
GFF_FILE=/datas/ELIOTT/archaea_data/test/genome/${ARCHAEA}.gff3













########## Phylostratigraphy ##########
echo "Running GenEra" >> $LOG_OUTPUT
genEra -t $TAXID -q ${OUT_DIR}CDS/${ARCHAEA}_CDS.fna -a $TAXID_FILE -n 32 -b $NR -d $TAXDUMP -o {$OUT_DIR}out/

# Deactivate conda environment
conda deactivate
echo "Job completed successfully" >> $LOG_OUTPUT