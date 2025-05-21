#! /bin/sh
#PBS -N reblast_rank1
#PBS -q bim
#PBS -l ncpus=8 -l host=node04 -l mem=64gb -l walltime=50:00:00
#PBS -o /home/eliott.tempez/reblast_rank1_output.log
#PBS -e /home/eliott.tempez/reblast_rank1_error.log

cd /datas/ELIOTT/scripts
cp cp /store/EQUIPES/BIM/MEMBERS/eliott.tempez/stage/M2_stage_I2BC/scripts/explore_genera_results/re_blast_rank_1.py .
source /home/eliott.tempez/miniconda3/bin/activate run_diamond
python re_blast_rank_1.py >> /home/eliott.tempez/reblast_rank1_output.log 2>> /home/eliott.tempez/reblast_rank1_error.log
source /home/eliott.tempez/miniconda3/bin/deactivate
rm re_blast_rank_1.py