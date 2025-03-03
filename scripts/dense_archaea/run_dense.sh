#!/bin/sh

########## Parameters ##########

# Get the species from the argument
archaea=$SPECIES

# Data files
GENDIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/complete_122/gendir_for_dense/
TREE=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/whole_tree.nwk
TAXID_FILE=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/taxid.csv
GENERA_OUTFILE=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/genera/out/${archaea}/*_gene_ages.tsv
OUT_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/dense/${archaea}/
# Log filenames
LOG_OUTPUT=/home/eliott.tempez/dense_output_$archaea.log
LOG_ERROR=/home/eliott.tempez/dense_error_$archaea.log

# Set up environment
source /home/eliott.tempez/miniconda3/bin/activate dense
mkdir -p /scratchlocal/$USER/$SLURM_JOBID
cd /scratchlocal/$USER/$SLURM_JOBID
# Copy all inputs
cp -r $GENDIR .
cp $TREE .
cp $TAXID_FILE .
cp $GENERA_OUTFILE .
GENERA_OUT_TMP=*_gene_ages.tsv
# Create temp output directory
mkdir -p out/



########## Run Dense ##########

echo $archaea >> $LOG_OUTPUT



# Run Dense
echo "Running Dense..." >> $LOG_OUTPUT
nextflow run /home/eliott.tempez/dense \
    -profile singularity \
    --max_cpus 8 \
    --max_memory 32.GB \
    --max_time 10.h \
    --num_outgroups 2 \
    --gendir gendir_for_dense/ \
    --focal $archaea \
    --tree whole_tree.nwk \
    --taxids taxid.csv \
    --genera_out "*gene_ages.tsv" \
    --trg_node Thermococcaceae \
    --outdir out/ >> $LOG_OUTPUT 2>> $LOG_ERROR


# Copy output
mkdir -p $OUT_DIR
cp -r out/* $OUT_DIR

# Deactivate conda environment
rm -r /scratchlocal/$USER/$SLURM_JOBID
conda deactivate
echo -e "Job completed successfully\n" >> $LOG_OUTPUT