#!/usr/bin/env bash
#author:    :Gregory Wickham
#date:      :20241119
#version    :1.0.0
#desc       :Align reads against reference
#usage		:bash bowtie2_alignment.sh
#===========================================================================================================
# Exit script on error
set -exuo pipefail

source "$(sudo find ~ -maxdepth 4 -name conda.sh)" #find path to conda base environment
conda activate bowtie2

# Set variables 
READS_DIR=$1        # Directory containing read files
BOWTIE2_REF=$2      # Path to ref file
OUTPUT_DIR=$3       # Directory for output BAM files

# Check if OUTPUT_DIR exists, if not create it
if [ ! -e "$OUTPUT_DIR" ]; then
    mkdir -p $(realpath $OUTPUT_DIR)
    echo "$OUTPUT_DIR not found, creating new directory."
fi

# Create the Bowtie2 ref index if not already created
if [ -f "${BOWTIE2_REF}.1.bt2" ]; then
    echo "Bowtie2 index $(basename $BOWTIE2_REF) found"
elif [ -f "${BOWTIE2_REF%.f*}".f* ]; then
    echo "Creating new index $BOWTIE2_REF"
    bowtie2-build $BOWTIE2_REF "${BOWTIE2_REF%.f*}"
    BOWTIE2_REF="${BOWTIE2_REF%.f*}"
else 
    echo "ERROR: File $BOWTIE2_REF not found. Exiting"
    exit
fi

#loop through READS_DIR
for R1 in $READS_DIR/*1.fastq*; do
    # Check R2 pairmate exsits, if not write to errors.txt/.,
    R2=${R1/1.fastq/2.fastq}  # Get R2 file name from R1
    if [ ! -f "$R2" ]; then
        echo "ERROR: No R2 file for $(basename $R1). Skipping read."
        > $OUTPUT_DIR/errors.txt
        echo "$R1 missing R2 pairmate" >> $OUTPUT_DIR/errors.txt
        exit 1
    fi

    # Extract sample name (e.g., sample from sample_R1.fastq)
    for filename in $R1; do
        if [[ $filename == *R1.fastq.gz ]]; then
            # Remove the *_R1.fastq.gz part
            SAMPLE_NAME=$(basename ${filename%_R1.fastq*})
        else
            # Remove the *_1.fastq.gz part
            SAMPLE_NAME=$(basename ${filename%_1.fastq*})
        fi
    done

    # Align reads and generate BAM file
    echo "Aligning $SAMPLE_NAME with bowtie2"
    bowtie2 \
        -x $BOWTIE2_REF \
        -1 $R1 -2 $R2 \
        | samtools view -bS - > "$OUTPUT_DIR/$SAMPLE_NAME.bam"  # pipe to samtools for sam to bam conversion

    echo "Alignment complete. Output written to $OUTPUT_DIR/$SAMPLE_NAME.bam"
done
