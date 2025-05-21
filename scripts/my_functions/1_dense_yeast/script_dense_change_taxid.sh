#!/bin/sh
#PBS -N dense_change_taxid
#PBS -q bim
#PBS -l ncpus=32 -l host=node04 -l mem=128gb -l walltime=240:00:00
#PBS -o /home/eliott.tempez/dense_change_taxid_output.log
#PBS -e /home/eliott.tempez/dense_change_taxid_error.log

# Set up environment
echo "Activating conda environment" >> /home/eliott.tempez/dense_change_taxid_output.log
source /home/eliott.tempez/miniconda3/bin/activate dense
cd /datas/ELIOTT/

# Get data 1
echo "---- First way ----" >> /home/eliott.tempez/dense_change_taxid_output.log
mkdir /datas/ELIOTT/result_dense_yeast_3/
rm -rf /datas/ELIOTT/yeast_data/*
cp -r /store/EQUIPES/BIM/MEMBERS/eliott.tempez/yeast_data/dense_change_taxid/* /datas/ELIOTT/yeast_data/
echo "Data imported" >> /home/eliott.tempez/dense_change_taxid_output.log

# Run Dense 1
# Run DENSE
echo "Running Nextflow" >> /home/eliott.tempez/dense_change_taxid_output.log
nextflow run /home/eliott.tempez/dense -profile singularity -c /home/eliott.tempez/params_yeast_change_taxid-1.config >> /home/eliott.tempez/dense_change_taxid_output.log 2>> /home/eliott.tempez/dense_change_taxid_error.log
EXIT_CODE=$?

# Check exit code
if [ $EXIT_CODE -ne 0 ]; then
    echo "Nextflow run failed with exit code $EXIT_CODE" >> /home/eliott.tempez/dense_change_taxid_error.log
    exit $EXIT_CODE
fi




# Get data 2
echo "----- Second way -----" >> /home/eliott.tempez/dense_change_taxid_output.log
mkdir /datas/ELIOTT/result_dense_yeast_4/
rm -rf /datas/ELIOTT/yeast_data/*
cp -r /store/EQUIPES/BIM/MEMBERS/eliott.tempez/yeast_data/dense_change_taxid_2/* /datas/ELIOTT/yeast_data/
echo "Data imported" >> /home/eliott.tempez/dense_change_taxid_output.log

# Run Dense 2
# Run DENSE
echo "Running Nextflow" >> /home/eliott.tempez/dense_change_taxid_output.log
nextflow run /home/eliott.tempez/dense -profile singularity -c /home/eliott.tempez/params_yeast_change_taxid-2.config >> /home/eliott.tempez/dense_change_taxid_output.log 2>> /home/eliott.tempez/dense_change_taxid_error.log
EXIT_CODE=$?

# Check exit code
if [ $EXIT_CODE -ne 0 ]; then
    echo "Nextflow run failed with exit code $EXIT_CODE" >> /home/eliott.tempez/dense_change_taxid_error.log
    exit $EXIT_CODE
fi

# Deactivate conda environment
echo "Deactivating conda environment" >> /home/eliott.tempez/dense_change_taxid_output.log
conda deactivate
echo "Job completed successfully" >> /home/eliott.tempez/dense_change_taxid_output.log


