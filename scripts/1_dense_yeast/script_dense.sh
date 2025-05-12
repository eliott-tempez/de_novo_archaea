#! /bin/sh
#!/bin/bash
#SBATCH --job-name=dense_yeast
#SBATCH --partition=bim
#SBATCH --nodelist=node04
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --time=200:00:00
#SBATCH --output=/home/eliott.tempez/dense_yeast_output.log
#SBATCH --error=/home/eliott.tempez/dense_yeast_error.log

# Log files
LOG_FILE="/home/eliott.tempez/dense_yeast_output.log"
ERROR_FILE="/home/eliott.tempez/dense_yeast_error.log"

# Set up environment
echo "Creating environment" >> $LOG_FILE
source /home/eliott.tempez/miniconda/bin/activate dense


# Copy data to node
DENSE_GENOMES_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/yeast_data/dense_initial/
PARAMS=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/stage/M2_stage_I2BC/scripts/1_dense_yeast/params_yeast.config
TREE=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/yeast_data/dense_initial/Saccharomyces_species.nwk
TAXIDS=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/yeast_data/dense_initial/taxids.tsv
DENSE=/home/eliott.tempez/dense

mkdir -p /datas/ELIOTT/
cd /datas/ELIOTT/
mkdir -p dense_genomes
cp -r $DENSE_GENOMES_DIR/* dense_genomes/
cp $TREE .
cp $TAXIDS .
mkdir -p dense/
cp -r $DENSE/* dense/
cp $PARAMS .


# Run DENSE
echo "Running Nextflow" >> $LOG_FILE
nextflow run ./dense -profile singularity -c params_yeast.config >> $LOG_FILE 2>> $ERROR_FILE
EXIT_CODE=$?

# Check exit code
if [ $EXIT_CODE -ne 0 ]; then
    echo "Nextflow run failed with exit code $EXIT_CODE" >> $ERROR_FILE
    exit $EXIT_CODE
fi


# Deactivate conda environment
conda deactivate
echo "Job completed successfully" >> $LOG_FILE