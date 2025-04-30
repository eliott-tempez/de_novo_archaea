#! /bin/sh

#SBATCH -p common
#SBATCH -J pval_descriptors
#SBATCH --time=10-00
#SBATCH --cpus-per-task=32
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH -o /home/eliott.tempez/pval_comparison_output.log
#SBATCH -e /home/eliott.tempez/pval_comparison_error.log

SCRIPTS=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/stage/M2_stage_I2BC/scripts
GOOD_CANDIDATES=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/stage/M2_stage_I2BC/results/14_get_noncoding_match/good_candidates.txt
CONTAINER=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/programs/ORFmine/orfmine_latest.sif
IUPRED_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/programs/iupred
TANGO=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/programs/tango/tango_x86_64_release
GENERA_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/genera/out
DENSE_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/dense
FA_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/complete_122/fasta_renamed
CDS_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/dense/GCA_000007305@Pyrococcus_furiosus_DSM_3638/CDS
GFF_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/complete_122/reannotated_gff_75
OUTPUT_LOG=/home/eliott.tempez/pval_comparison_output.log
ERROR_LOG=/home/eliott.tempez/pval_comparison_error.log
OUT_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/pvals_descriptors

# Set up a trap to clean up on exit or interruption
trap "rm -rf $SCRATCH_DIR; find /tmp -maxdepth 1 -name "*" -user eliott.tempez -exec rm -rf {} +" EXIT

# Create environment
source /home/eliott.tempez/miniconda3/bin/activate descriptors
module load singularity

SCRATCH_DIR=/scratchlocal/$USER/$SLURM_JOBID
SOFTWARE_PATH=$SCRATCH_DIR/softwares
mkdir -p $SCRATCH_DIR
mkdir -p $SOFTWARE_PATH
cd $SCRATCH_DIR
cp -r $SCRIPTS/* .
cp $GOOD_CANDIDATES .
cp $CONTAINER .
# Softwares
cp -r $IUPRED_DIR/* $SOFTWARE_PATH
cp $TANGO $SOFTWARE_PATH
mkdir -p genera/
cp -r $GENERA_DIR/* genera/
mkdir -p dense/
rsync -a --exclude 'archive*/' --include '*/' --include 'blast_out/**' --include '*/' --include '*/[^/]*.tsv' --exclude '*' $DENSE_DIR .
mkdir -p cds/
cp -r $CDS_DIR/* cds/
mkdir -p fa/
cp -r $FA_DIR/* fa/
mkdir -p gff/
cp -r $GFF_DIR/* gff/
mkdir -p out/
echo "Environment created" >> $OUTPUT_LOG

# Launch program
python 13_plot_dense_results/compare_sequences_statistically.py >> $OUTPUT_LOG

# Get output 
cp -f out/* $OUT_DIR

# Clean up
rm -rf $SCRATCH_DIR
find /tmp -maxdepth 1 -name "*" -user eliott.tempez -exec rm -rf {} +
conda deactivate
echo "Done" >> $OUTPUT_LOG