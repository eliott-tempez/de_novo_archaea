#!/bin/sh

########## Parameters ##########

#PBS -N CDS_archaea_test
#PBS -q bim
#PBS -l ncpus=4 -l host=node04 -l mem=128gb -l walltime=200:00:00
#PBS -o /home/eliott.tempez/extract_cds_output.log
#PBS -e /home/eliott.tempez/extract_cds_error.log

# Log files
LOG_OUTPUT=/home/eliott.tempez/extract_cds_output.log
LOG_ERROR=/home/eliott.tempez/extract_cds_error.log

# Data files
NR=/datas/NR/nr_2.0.13.dmnd
TAXDUMP=/datas/ELIOTT/work/taxdump
TAXID_FILE=/datas/ELIOTT/archaea_data/test/taxid.csv
TREE=/datas/ELIOTT/archaea_data/test/small-tree.nwk

# Output dir
OUT_DIR=/datas/ELIOTT/archaea_data/test/genera/

# Set up environment
source /home/eliott.tempez/miniconda3/bin/activate phylostrat
cd /datas/ELIOTT/



########## Extract CDS ##########

for file in ./archaea_data/test/genome/*.fa; do
    ARCHAEA=$(basename $file .fa)
    echo $ARCHAEA >> $LOG_OUTPUT
    # Get taxid
    TAXID=$(grep $ARCHAEA $TAXID_FILE | cut -d, -f2)

    # Get filenames
    FASTA_FILE=/datas/ELIOTT/archaea_data/test/genome/${ARCHAEA}.fa
    GFF_FILE=/datas/ELIOTT/archaea_data/test/genome/${ARCHAEA}.gff3

    # fai index
    echo "Indexing $FASTA_FILE..." >> $LOG_OUTPUT
    samtools faidx $FASTA_FILE
    # Keep GFF sequences that are also in the FASTA file.
    awk 'BEGIN{FS=OFS="\t"} {print "^" $1,""} END {print "^#"}' ${FASTA_FILE}.fai | grep -f - $GFF_FILE > gff_filterA
    # Also remove mRNA with undefined strand and features with abnormal end
    awk -F"\t" '
            FNR==NR {max[\$1]=\$2}
            FNR!= NR && ( /^#/ || (\$7 !~ /?/ && \$5 <= max[\$1]) )
            ' ${FASTA_FILE}.fai gff_filterA > gff_filterB
    echo "Extracting CDS..." >> $LOG_OUTPUT
    # Extract the genomic CDS.
    # -V discard any mRNAs with CDS having in-frame stop codons
    # -x : write a fasta file with spliced CDS for each GFF transcript
    gffread -V -g $FASTA_FILE -x ${OUT_DIR}CDS/${ARCHAEA}_CDS.fna gff_filterB >> $LOG_OUTPUT 2>> $LOG_ERROR
    rm gff_filterA gff_filterB
    # Remove CDS missing a terminal stop codon and get a translated version of the FASTA.
    echo "Discarding CDS missing terminal stop codon..." >> $LOG_OUTPUT
    /datas/ELIOTT/discard_CDS_missing_terminal_stop_codon.sh ${OUT_DIR}CDS/${ARCHAEA}_CDS.fna >> $LOG_OUTPUT 2>> $LOG_ERROR
    echo "done!" >> $LOG_OUTPUT
done

