#!/bin/sh

########## Parameters ##########

#PBS -N extract_CDS
#PBS -q bim
#PBS -l ncpus=4 -l host=node04 -l mem=16gb -l walltime=10:00:00
#PBS -o /home/eliott.tempez/extract_cds_output.log
#PBS -e /home/eliott.tempez/extract_cds_error.log

# Log files
LOG_OUTPUT=/home/eliott.tempez/extract_cds_output.log
LOG_ERROR=/home/eliott.tempez/extract_cds_error.log

# Data files
TAXID_FILE=/datas/ELIOTT/archaea_data/taxid.csv
GENOMES_FOLDER=/datas/ELIOTT/archaea_data/fasta_renamed/
GFF_FOLDER=/datas/ELIOTT/archaea_data/reannotated_gff_75/

# Output dir
OUT_DIR=/datas/ELIOTT/archaea_data/reannotated_CDS/

# Set up environment
source /home/eliott.tempez/miniconda3/bin/activate phylostrat
cd /datas/ELIOTT/scripts/
# copy scripts
cp /home/eliott.tempez/dense/bin/discard_CDS_missing_terminal_stop_codon.sh .


########## Extract CDS ##########

for file in $GENOMES_FOLDER*.fa
do
    ARCHAEA=$(basename $file .fa)
    echo $ARCHAEA >> $LOG_OUTPUT
    # Get taxid
    TAXID=$(grep $ARCHAEA $TAXID_FILE | cut -d, -f2)

    # Get filenames
    FASTA_FILE=$file
    GFF_FILE=$GFF_FOLDER${ARCHAEA}.gff3
    # Create output dir
    mkdir -p ${OUT_DIR}

    # fai index
    echo "Indexing $FASTA_FILE..." >> $LOG_OUTPUT
    samtools faidx $FASTA_FILE
    # Keep GFF sequences that are also in the FASTA file.
    awk 'BEGIN{FS=OFS="\t"} {print "^" $1,""} END {print "^#"}' ${FASTA_FILE}.fai | grep -f - $GFF_FILE > gff_filterA
    # Also remove mRNA with undefined strand and features with abnormal end
    awk -F"\t" '
            FNR==NR {max[$1]=$2}
            FNR!= NR && ( /^#/ || ($7 !~ /?/ && $5 <= max[$1]) )
            ' ${FASTA_FILE}.fai gff_filterA > gff_filterB
    echo "Extracting CDS..." >> $LOG_OUTPUT
    # Extract the genomic CDS.
    # -V discard any mRNAs with CDS having in-frame stop codons
    # -x : write a fasta file with spliced CDS for each GFF transcript
    gffread -V -g $FASTA_FILE -x ${OUT_DIR}CDS/${ARCHAEA}_CDS.fna gff_filterB >> $LOG_OUTPUT 2>> $LOG_ERROR
    rm gff_filterA gff_filterB
    # Remove CDS missing a terminal stop codon and get a translated version of the FASTA.
    echo "Discarding CDS missing terminal stop codon..." >> $LOG_OUTPUT
    /datas/ELIOTT/scripts/discard_CDS_missing_terminal_stop_codon.sh ${OUT_DIR}CDS/${ARCHAEA}_CDS.fna >> $LOG_OUTPUT 2>> $LOG_ERROR
    echo "done!" >> $LOG_OUTPUT
done

# Copy output 
cp ${OUT_DIR}/* /store/EQUIPES/BIM/MEMBERS/eliott.tempez/archaea_data/complete_122/reannotated_CDS/
