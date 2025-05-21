#!/bin/sh

########## Parameters ##########

# Get the species from the argument
archaea=$SPECIES

# Data files
OUT_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/complete_122/gff_reannotated/${archaea}/
GENOME=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/complete_122/fasta_renamed/${archaea}.fa
# Log filenames
LOG_OUTPUT=/home/eliott.tempez/prokka_output_$archaea.log
LOG_ERROR=/home/eliott.tempez/prokka_error_$archaea.log


# Set up environment
source /home/eliott.tempez/miniconda3/bin/activate prokka
# Create a unique scratch directory using PBS job ID
SCRATCH_DIR=/scratchlocal/$USER/$PBS_JOBID
mkdir -p $SCRATCH_DIR
cd $SCRATCH_DIR
# Copy all inputs
cp $GENOME .
# Create temp output directory
mkdir -p out/




########## Run Prokka ##########
prokka --outdir out/ --prefix $archaea --cpus 0 --force --kingdom Archaea --addgenes $GENOME >> $LOG_OUTPUT 2>> $LOG_ERROR

# Copy output
mkdir -p $OUT_GFF
cp -r out/* $OUT_GFF

# Clean up
rm -rf $SCRATCH_DIR
conda deactivate