#!/bin/bash

SAMPLE_NAME="$1"

INPUT_DIR="2_exp_fastq"
OUTPUT_DIR="4_stats"
OUTPUT_SAMPLE_DIR="$OUTPUT_DIR/$SAMPLE_NAME/"

# add directory to PATH environment
export PATH="utils/quack/:$PATH"

# Ensure the output directory exists
mkdir -p "$OUTPUT_SAMPLE_DIR"

# Run Quack on the FastQ files
quack -1 "$INPUT_DIR/${SAMPLE_NAME}_R1.fastq.gz" -2 "$INPUT_DIR/${SAMPLE_NAME}_R2.fastq.gz" -n $SAMPLE_NAME > $OUTPUT_SAMPLE_DIR/${SAMPLE_NAME}_quack_QC.svg

echo -e "FASTQ QC for ${SAMPLE_NAME} done..."

