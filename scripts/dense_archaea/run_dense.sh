#!/bin/sh

########## Parameters ##########

#PBS -N dense_archaea
#PBS -q bim
#PBS -l ncpus=32 -l host=node04 -l mem=128gb -l walltime=300:00:00

# Data files
GENDIR=/datas/ELIOTT/archaea_data/genome/
TREE=/datas/ELIOTT/archaea_data/whole_tree.nwk
TAXID_FILE=/datas/ELIOTT/archaea_data/taxid.csv
NR=/datas/NR/nr_2.0.13.dmnd

# Set up environment
source /home/eliott.tempez/miniconda3/bin/activate dense
cd /datas/ELIOTT/scripts/



########## Run Dense ##########

# List of archaeas to iterate over
declare -a archaeas=("GCA_000195935@Pyrococcus_abyssi_GE5")
# For each of them
for archaea in "${archaeas[@]}"; do

    echo $archaea >> $LOG_OUTPUT

    # Log filenames
    LOG_OUTPUT=/home/eliott.tempez/dense_output_${archaea}.log
    LOG_ERROR=/home/eliott.tempez/dense_error_${archaea}.log

    # Get taxid
    TAXID=$(grep $archaea $TAXID_FILE | cut -f2 -d$'\t')
    echo "Taxid: $TAXID" >> $LOG_OUTPUT
    # Get filenames
    OUT_DIR=/datas/ELIOTT/archaea_data/dense/${archaea}/


    # Run Dense
    echo "Running Dense..." >> $LOG_OUTPUT
    nextflow run /home/eliott.tempez/dense \
        -profile singularity \
        --max_cpus 32 \
        --max_memory 128.GB \
        --max_time 300.h \
        --num_outgroups 2 \
        --gendir $GENDIR \
        --focal $archaea \
        --tree $TREE \
        --taxids $TAXID_FILE \
        --genera_db $NR \
        --trg_node Thermococcaceae \
        --outdir $OUT_DIR >> $LOG_OUTPUT 2>> $LOG_ERROR


    # Deactivate conda environment
    conda deactivate
    echo -e "Job completed successfully\n" >> $LOG_OUTPUT

done