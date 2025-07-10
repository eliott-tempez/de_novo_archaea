#! /bin/sh
# /!\ The paths are probably not correct, you will have to change them

#SBATCH -p common
#SBATCH -J run_descriptors
#SBATCH --time=10-00
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --mem=50gb
#SBATCH -o /home/eliott.tempez/comparison_output.log
#SBATCH -e /home/eliott.tempez/comparison_error.log

SCRIPTS=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/stage/M2_stage_I2BC/scripts
GOOD_CANDIDATES=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/stage/M2_stage_I2BC/results/14_get_noncoding_match/good_candidates.txt
IORF_FILE=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/stage/M2_stage_I2BC/results/17_compare_denovo_sequences/iorfs.txt
INTERGENIC_GC=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/stage/M2_stage_I2BC/results/10_analyse_intergenic/intergenic_gc.tsv
CONTAINER=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/programs/ORFmine/orfmine_latest.sif
IUPRED_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/programs/iupred
TANGO=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/programs/tango/tango_x86_64_release
GENERA_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/genera/out
DENSE_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/dense
FA_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/complete_122/fasta_renamed
CDS_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/dense/GCA_000007305@Pyrococcus_furiosus_DSM_3638/CDS
GFF_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/complete_122/reannotated_gff_75
OUTPUT_LOG=/home/eliott.tempez/comparison_output.log
ERROR_LOG=/home/eliott.tempez/comparison_error.log
OUT_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/pvals_descriptors

# Set up a trap to clean up on exit or interruption
trap "rm -rf $SCRATCH_DIR; find /tmp -maxdepth 1 -name "*" -user eliott.tempez -exec rm -rf {} +" EXIT SIGTERM SIGINT

# Create environment
source /home/eliott.tempez/miniconda3/bin/activate descriptors
module load singularity

SCRATCH_DIR=/scratchlocal/$USER/$SLURM_JOBID
mkdir -p $SCRATCH_DIR
cd $SCRATCH_DIR
cp -r $SCRIPTS/* .
cp $GOOD_CANDIDATES .
cp $IORF_FILE .
cp $INTERGENIC_GC .
cp $CONTAINER .
# Softwares
cp -r $IUPRED_DIR/* .
cp $TANGO .
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
python 17_compare_denovo_sequences/compare_sequences.py >> $OUTPUT_LOG

# Get output 
cp -f out/* $OUT_DIR

# Clean up
rm -rf $SCRATCH_DIR
find /tmp -maxdepth 1 -name "*" -user eliott.tempez -exec rm -rf {} +
conda deactivate
echo "Done" >> $OUTPUT_LOG