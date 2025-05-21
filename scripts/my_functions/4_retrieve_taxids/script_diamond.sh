#! /bin/sh
#PBS -N diamond_test
#PBS -q bim
#PBS -l ncpus=32 -l host=node04 -l mem=128gb -l walltime=240:00:00
#PBS -o /home/eliott.tempez/diamond_test_output.log
#PBS -e /home/eliott.tempez/diamond_test_error.log

cd /datas/ELIOTT/archaea_data/retrieve_taxid
cp -f /store/EQUIPES/BIM/MEMBERS/eliott.tempez/stage/M2_stage_I2BC/scripts/retrieve_taxids/retrieve_taxids.py .
source /home/eliott.tempez/miniconda3/bin/activate get_taxids
python3 retrieve_taxids.py >> /home/eliott.tempez/diamond_test_output.log 2>> /home/eliott.tempez/diamond_test_error.log
conda deactivate