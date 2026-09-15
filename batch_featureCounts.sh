#!/bin/bash

set -euo pipefail

SECONDS=0

# step 0: set working directories
BAM_DIR="HISAT2/out/bam"
OUT_DIR="counts/out"
GTF="annotation/gencode.gtf"

echo "Starting featureCounts..."

# step 1: loop through all .bam files
for BAM in "$BAM_DIR"/*.sorted.bam
do
    # step 2: extract file names
    SAMPLE=$(basename "$BAM" .sorted.bam)

    # step 3: run featureCounts
    featureCounts \
        -a "$GTF" \
        -o "$OUT_DIR/${SAMPLE}.counts.tsv" \
        -T 8 \
        -p \
        -B \
        --format TSV \
        "$BAM"
done

echo "Finished. $(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds elapsed."