#! /bin/sh

#SBATCH -p common
#SBATCH -J pvals_gb
#SBATCH --time=1000-00
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH -o /home/eliott.tempez/pval_gb_bin_out.log
#SBATCH -e /home/eliott.tempez/pval_gb_bin_err.log

SCRIPTS=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/stage/M2_stage_I2BC/scripts
GOOD_CANDIDATES=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/stage/M2_stage_I2BC/results/14_get_noncoding_match/good_candidates.txt
DESCRIPTORS=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/pvals_descriptors/sequence_features_good_candidates_all.csv
FA_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/complete_122/fasta_renamed
DENSE_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/dense
OUT_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/pvals_descriptors
OUTPUT_LOG=/home/eliott.tempez/pval_gb_bin_out.log
ERROR_LOG=/home/eliott.tempez/pval_gb_bin_err.log


# Set up a trap to clean up on exit or interruption
trap "rm -rf $SCRATCH_DIR; find /tmp -maxdepth 1 -name "*" -user eliott.tempez -exec rm -rf {} +" EXIT SIGTERM SIGINT

# Create environment
source /home/eliott.tempez/miniconda3/bin/activate descriptors
module load singularity

SCRATCH_DIR=/scratchlocal/$USER/$SLURM_JOBID
mkdir -p $SCRATCH_DIR
cd $SCRATCH_DIR
cp -r $SCRIPTS/* .
cp $DESCRIPTORS .
mkdir -p fa/
cp -r $FA_DIR/* fa/
mkdir -p dense/
rsync -a --exclude 'archive*/' --include '*/' --exclude '*/*/*' --include '*.tsv' --exclude '*' $DENSE_DIR/ dense/
cp $GOOD_CANDIDATES .
echo "Environment created" >> $OUTPUT_LOG

# Run script
python 17_compare_denovo_sequences/calculate_pvals_seq_comparison.py >> $OUTPUT_LOG

# Get output
cp pvalues*.csv $OUT_DIR
cp bin_indexes_*.csv $OUT_DIR
# Clean up
rm -rf $SCRATCH_DIR
find /tmp -maxdepth 1 -name "*" -user eliott.tempez -exec rm -rf {} +
conda deactivate
echo "Done" >> $OUTPUT_LOG