#!/bin/sh

########## Parameters ##########

# Get the species from the argument
archaea=$SPECIES

# Data files
GENDIR=/datas/ELIOTT/archaea_data/genome/
NR=/datas/NR/nr_2.0.13.dmnd
TREE=/datas/ELIOTT/archaea_data/whole_tree.nwk
TAXID_FILE=/datas/ELIOTT/archaea_data/taxid.csv
GENERA_OUTFILE=/datas/ELIOTT/archaea_data/genera/out/${archaea}/*_gene_ages.tsv
OUT_DIR=/datas/ELIOTT/archaea_data/dense/${archaea}/
# Log filenames
LOG_OUTPUT=/home/eliott.tempez/dense_output_$archaea.log
LOG_ERROR=/home/eliott.tempez/dense_error_$archaea.log

# Set up environment
source /home/eliott.tempez/miniconda3/bin/activate dense
cd /datas/ELIOTT/scripts/



########## Run Dense ##########

echo $archaea >> $LOG_OUTPUT
mkdir -p $OUT_DIR


# Run Dense
echo "Running Dense..." >> $LOG_OUTPUT
nextflow run /home/eliott.tempez/dense \
    -profile singularity \
    --max_cpus 8 \
    --max_memory 16.GB \
    --max_time 10.h \
    --num_outgroups 2 \
    --gendir $GENDIR \
    --focal $archaea \
    --tree $TREE \
    --taxids $TAXID_FILE \
    --genera_out $GENERA_OUTFILE \
    --trg_node Thermococcaceae \
    --outdir $OUT_DIR >> $LOG_OUTPUT 2>> $LOG_ERROR


# Deactivate conda environment
conda deactivate
echo -e "Job completed successfully\n" >> $LOG_OUTPUT