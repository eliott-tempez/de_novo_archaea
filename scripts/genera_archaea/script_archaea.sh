#!/bin/sh


########## Parameters ##########
ARCHAEA=GCA_000007305@Pyrococcus_furiosus_DSM_3638
TAXID=186497

#PBS -N genera_archaea_test
#PBS -q bim
#PBS -l ncpus=32 -l host=node04 -l mem=128gb -l walltime=200:00:00
#PBS -o /home/eliott.tempez/logs_genera/${ARCHAEA}_output.log
#PBS -e /home/eliott.tempez/logs_genera/${ARCHAEA}_error.log

# Log filenames
LOG_OUTPUT=/home/eliott.tempez/logs_genera/${ARCHAEA}_output.log
LOG_ERROR=/home/eliott.tempez/logs_genera/${ARCHAEA}_error.log

# Data files
NR=/datas/NR/nr_2.0.13.dmnd
TAXDUMP=/datas/ELIOTT/work/taxdump
TAXID_FILE=/datas/ELIOTT/archaea_data/taxid_fastapath.csv
FASTA_FILE=/datas/ELIOTT/archaea_data/genome/${ARCHAEA}.fa
GFF_FILE=/datas/ELIOTT/archaea_data/genome/${ARCHAEA}.gff3

# Output dir
OUT_DIR=/datas/ELIOTT/archaea_data/genera/

# Set up environment
echo $ARCHAEA >> $LOG_OUTPUT
source /home/eliott.tempez/miniconda3/bin/activate phylostrat
cd /datas/ELIOTT/



########## Extract CDS ##########
echo "Extracting CDS" >> $LOG_OUTPUT
# fai index
samtools faidx $FASTA_FILE
# Keep GFF sequences that are also in the FASTA file.
awk 'BEGIN{FS=OFS="\t"} {print "^"\$1,""} END {print "^#"}' ${FASTA_FILE}.fai | grep -f - $GFF_FILE > gff_filterA
# Also remove mRNA with undefined strand and features with abnormal end
awk -F"\t" '
        FNR==NR {max[\$1]=\$2}
        FNR!= NR && ( /^#/ || (\$7 !~ /?/ && \$5 <= max[\$1]) )
        ' ${FASTA_FILE}.fai gff_filterA > gff_filterB
# Extract the genomic CDS.
# -V discard any mRNAs with CDS having in-frame stop codons
# -x : write a fasta file with spliced CDS for each GFF transcript
gffread -V -g $FASTA_FILE -x ${OUT_DIR}CDS/${ARCHAEA}_CDS.fna gff_filterB >> $LOG_OUTPUT 2>> $LOG_ERROR
# Remove CDS missing a terminal stop codon and get a translated version of the FASTA.
./discard_CDS_missing_terminal_stop_codon.sh ${OUT_DIR}CDS/${ARCHAEA}_CDS.fna >> $LOG_OUTPUT 2>> $LOG_ERROR
rm gff_filterA gff_filterB
echo "CDS extracted" >> $LOG_OUTPUT



########## Phylostratigraphy ##########
echo "Running GenEra" >> $LOG_OUTPUT
genEra -t $TAXID -q ${OUT_DIR}${ARCHAEA}_CDS.fna -a $TAXID_FILE -n 32 -b $NR -d $TAXDUMP -o $OUT_DIR

# Deactivate conda environment
conda deactivate
echo "Job completed successfully" >> $LOG_OUTPUT