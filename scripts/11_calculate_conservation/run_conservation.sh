#! /bin/sh

#SBATCH -p run2
#SBATCH -J conservation_calculations
#SBATCH --time=3-00
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH -o /home/eliott.tempez/conservation_calculations.log
#SBATCH -e /home/eliott.tempez/conservation_calculations.log

SCRIPTS=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/stage/M2_stage_I2BC/scripts
GENERA_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/genera/out
DENSE_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/dense/dense_files_only
FA_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/complete_122/fasta_renamed
GBK_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/complete_122/reannotated_gbk_75
OUTPUT_LOG=/home/eliott.tempez/conservation_calculations.log
ERROR_LOG=/home/eliott.tempez/conservation_calculations.log
OUT_DIR=/store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/conservation

# Create environment
SCRATCH_DIR=/scratchlocal/$USER/$SLURM_JOBID
mkdir -p $SCRATCH_DIR
cd $SCRATCH_DIR
cp -r $SCRIPTS/* .
mkdir -p genera/
cp -r $GENERA_DIR/* genera/
mkdir -p dense/
cp -r $DENSE_DIR/* dense/
mkdir -p cds/
cp dense/GCA_000007*/CDS/* cds/
mkdir -p fa/
cp -r $FA_DIR/* fa/
mkdir -p gbk/
cp -r $GBK_DIR/* gbk/
mkdir -p out/
echo "Environment created" >> $OUTPUT_LOG

# Set up a trap to clean up on exit or interruption
trap "rm -rf $SCRATCH_DIR; find /tmp -maxdepth 1 -type d -name 'rootfs-*' -user eliott.tempez -exec rm -rf {} +; find /tmp -maxdepth 1 -name 'tmp*' -user eliott.tempez -exec rm -rf {} +; cp -r out/* $OUT_DIR" EXIT




# Launch calculation
source /home/eliott.tempez/miniconda3/bin/activate python_basics
declare -a archaeas=("GCA_000007305@Pyrococcus_furiosus_DSM_3638" "GCA_000009965@Thermococcus_kodakarensis_KOD1" "GCA_000011105@Pyrococcus_horikoshii_OT3" "GCA_000018365@Thermococcus_onnurineus_NA1" "GCA_000022365@Thermococcus_gammatolerans_EJ3" "GCA_000022545@Thermococcus_sibiricus_MM_739" "GCA_000151105@Thermococcus_barophilus_MP" "GCA_000151205@Thermococcus_sp_AM4" "GCA_000195935@Pyrococcus_abyssi_GE5" "GCA_000211475@Pyrococcus_sp_NA2" "GCA_000215995@Pyrococcus_yayanosii_CH1" "GCA_000221185@Thermococcus_sp_4557" "GCA_000246985@Thermococcus_litoralis_DSM_5473" "GCA_000258515@Thermococcus_zilligii_AN1" "GCA_000263735@Pyrococcus_sp_ST04" "GCA_000265525@Thermococcus_cleftensis_CL1" "GCA_000275605@Pyrococcus_furiosus_COM1" "GCA_000517445@Thermococcus_paralvinellae_ES1" "GCA_000585495@Thermococcus_nautili_30_1" "GCA_000725425@Palaeococcus_pacificus_DY20341" "GCA_000769655@Thermococcus_eurythermalis_A501" "GCA_000816105@Thermococcus_guaymasensis_DSM_11113" "GCA_001433455@Thermococcus_barophilus_CH5" "GCA_001484685@Thermococcus_sp_2319x1" "GCA_001577775@Pyrococcus_kukulkanii_NCB100" "GCA_001592435@Thermococcus_peptonophilus_OG_1" "GCA_001647085@Thermococcus_piezophilus_CDGS" "GCA_002197185@Thermococcus_sp_5_4" "GCA_002214365@Thermococcus_celer_Vu_13" "GCA_002214385@Thermococcus_gorgonarius_W_12" "GCA_002214465@Thermococcus_barossii_SHCK_94" "GCA_002214485@Thermococcus_pacificus_P4" "GCA_002214505@Thermococcus_siculi_RG_20" "GCA_002214525@Thermococcus_sp_P6" "GCA_002214545@Thermococcus_thioreducens_OGL_20P" "GCA_002214565@Thermococcus_radiotolerans_EJ2" "GCA_002214585@Thermococcus_profundus_DT_5432" "GCA_002214605@Thermococcus_chitonophagus" "GCA_006274605@Thermococcus_indicus_IOH1" "GCA_008152015@Thermococcus_aciditolerans_SY113" "GCA_008245085@Pyrococcus_furiosus_DSM_3638_Reg" "GCA_020386975@Thermococcus_bergensis_T7324" "GCA_023746595@Thermococcus_argininiproducens_IOH2" "GCA_024022995@Thermococcus_aggregans_TY" "GCA_024707485@Thermococcus_sp_813A4" "GCA_028471785@Thermococcus_kodakarensis_TS900" "GCA_028471865@Thermococcus_kodakarensis_TS901" "GCA_035621495@Thermococcus_sp_SY098" "GCA_900012635@Thermococcus_chitonophagus_isolate_1" "GCA_900198835@Thermococcus_henrietii_EXT12c" "GCA_902813195@Thermococcus_sp_2319x1_Essen" "GCA_904067545@Thermococcus_camini_IRI35c" "Palaeococcus_helgesoni_DSM_15127" "Pyrococcus_31_3" "Pyrococcus_32_1" "Pyrococcus_abyssi_GE2" "Pyrococcus_endeavori_ES4" "Pyrococcus_EXT16" "Pyrococcus_sp_32_3" "Pyrococcus_sp_32_4" "Pyrococcus_sp_EXT_15c" "Pyrococcus_sp_GE_23" "Pyrococcus_sp_IRI_42C" "Thermococcus_15_2" "Thermococcus_23_2" "Thermococcus_29_3" "Thermococcus_31_3" "Thermococcus_aegaeus_DSM_12767" "Thermococcus_aggregans_DSM_12819" "Thermococcus_alcaliphilus_DSM_10322" "Thermococcus_AMTc11" "Thermococcus_AMTc30" "Thermococcus_AMTc51" "Thermococcus_AMTc85" "Thermococcus_atlanticus_DSM_15226" "Thermococcus_barophilus_423_3_23" "Thermococcus_barophilus_CH1" "Thermococcus_barophilus_DT4" "Thermococcus_coalescens_TS1" "Thermococcus_E15P35" "Thermococcus_EXT11c" "Thermococcus_fumicolans_ST557" "Thermococcus_IRI_06c" "Thermococcus_IRI15c" "Thermococcus_marinus_DSM_15227" "Thermococcus_prieuri_Col3" "Thermococcus_profundus_DSM_9503" "Thermococcus_sp_23_4" "Thermococcus_sp_26_2" "Thermococcus_sp_33_3" "Thermococcus_sp_9_3" "Thermococcus_sp_AMTc09" "Thermococcus_sp_AMTc102" "Thermococcus_sp_AMTc71" "Thermococcus_sp_AMTc72" "Thermococcus_sp_AMTc79" "Thermococcus_sp_AMTc94" "Thermococcus_sp_CIR_03a" "Thermococcus_sp_CIR_10a" "Thermococcus_sp_DS1" "Thermococcus_sp_E10P11" "Thermococcus_sp_E10P8" "Thermococcus_sp_E15p33" "Thermococcus_sp_EXT09c" "Thermococcus_sp_IRI_05c" "Thermococcus_sp_IRI10c" "Thermococcus_sp_IRI14c" "Thermococcus_sp_IRI24c" "Thermococcus_sp_IRI25c" "Thermococcus_sp_IRI27c2" "Thermococcus_sp_IRI29c" "Thermococcus_sp_IRI_33c" "Thermococcus_sp_TPV" "Thermococcus_sp_unknown" "Thermococcus_stetteri_DSM_5262" "Thermococcus_waiotapuensis_WT1")

for archaea in "${archaeas[@]}"; do
    python 11_calculate_conservation/calculate_conservation.py --focal_species $archaea >> $OUTPUT_LOG 2>> $ERROR_LOG
done

# Copy output
cp -r out/* $OUT_DIR
# Clean up
rm -rf $SCRATCH_DIR
find /tmp -maxdepth 1 -type d -name "rootfs-*" -user eliott.tempez -exec rm -rf {} +
find /tmp -maxdepth 1 -name "tmp*" -user eliott.tempez -exec rm -rf {} +
conda deactivate
echo "Done" >> $OUTPUT_LOG
