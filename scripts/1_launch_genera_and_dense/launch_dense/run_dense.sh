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

# Check if we already have the result
OUT_FILE=$OUT_DIR/denovogenes.tsv
if [ -f $OUT_FILE ]; then
    echo "Output file already exists, skipping job"
    exit 0
fi


# Set up environment
source /home/eliott.tempez/miniconda3/bin/activate dense
module load singularity

# Create a unique scratch directory using PBS job ID
SCRATCH_DIR=/scratchlocal/$USER/$SLURM_JOBID
# Set up Singularity environment
export SINGULARITY_TMPDIR=$SCRATCH_DIR/tmp
export SINGULARITY_CACHEDIR=$SCRATCH_DIR/cache
# Create directories
mkdir -p $SINGULARITY_TMPDIR $SINGULARITY_CACHEDIR $SCRATCH_DIR

cd $SCRATCH_DIR
# Delete old tmp files
find /tmp -maxdepth 1 -type d -name "rootfs-*" -user eliott.tempez -exec rm -rf {} +
# Copy all inputs
cp -r /home/eliott.tempez/dense .
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
nextflow run ./dense \
    -profile singularity \
    --max_cpus 8 \
    --max_memory 64.GB \
    --max_time 10.h \
    --num_outgroups 2 \
    --gendir gendir_for_dense/ \
    --focal $archaea \
    --tree whole_tree.nwk \
    --taxids taxid.csv \
    --genera_out $GENERA_OUT_TMP \
    --trg_node Thermococcaceae \
    --outdir out/ \>> $LOG_OUTPUT 2>> $LOG_ERROR


# Copy output
mkdir -p $OUT_DIR
cp -rL out/* $OUT_DIR

# Deactivate conda environment
cd /home/eliott.tempez
rm -r $SCRATCH_DIR
conda deactivate
echo -e "Job completed successfully\n" >> $LOG_OUTPUT