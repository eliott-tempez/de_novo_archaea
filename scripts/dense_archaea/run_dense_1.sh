#!/bin/sh

########## Parameters ##########

#PBS -N dense_archaea_1
#PBS -q bim
#PBS -l ncpus=8 -l host=node04 -l mem=64gb -l walltime=100:00:00
#PBS -o /home/eliott.tempez/genera_output.log
#PBS -e /home/eliott.tempez/genera_error.log

# Log filenames
LOG_OUTPUT=/home/eliott.tempez/dense_output.log
LOG_ERROR=/home/eliott.tempez/dense_error.log

# Data files
GENDIR=/datas/ELIOTT/archaea_data/genome/
TREE=/datas/ELIOTT/archaea_data/whole_tree.nwk
TAXDUMP=/datas/ELIOTT/scripts/taxdump/
TAXID_FILE=/datas/ELIOTT/archaea_data/taxid.csv

# Set up environment
source /home/eliott.tempez/miniconda3/bin/activate dense
cd /datas/ELIOTT/scripts/



########## Run Dense ##########

# List of archaeas to iterate over
declare -a archaeas=("GCA_000195935@Pyrococcus_abyssi_GE5")
# For each of them
for archaea in "${archaeas[@]}"; do

    echo $archaea >> $LOG_OUTPUT

    # Get taxid
    TAXID=$(grep $archaea $TAXID_FILE | cut -f2 -d$'\t')
    echo "Taxid: $TAXID" >> $LOG_OUTPUT
    # Get filenames
    GENERA_OUTFILE=/datas/ELIOTT/archaea_data/genera/out/${archaea}/*_gene_ages.tsv
    OUT_DIR=/datas/ELIOTT/archaea_data/dense/${archaea}/


    # Run Dense
    echo "Running Dense..." >> $LOG_OUTPUT
    nextflow run /home/eliott.tempez/dense \
        -profile singularity \
        --max_cpus 8 \
        --max_memory 64.GB \
        --max_time 100.h \
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

done