#!/bin/bash

SUBJ_ID=$1
POOL_DIR=$2

echo "Running ${SUBJ_ID}."

# 1. samtools sort -n -o out.bam in.bam (sort by query name)
samtools sort -n -o ${POOL_DIR}/${SUBJ_ID}.pooled_reads.sorted.bam ${POOL_DIR}/${SUBJ_ID}.pooled_reads.bam

# 2. samtools fixmate -O bam in.bam out.bam
samtools fixmate -O bam ${POOL_DIR}/${SUBJ_ID}.pooled_reads.sorted.bam ${POOL_DIR}/${SUBJ_ID}.pooled_reads.fix.bam

# 3. samtools sort -o out.bam in.bam (sort by coordinates)
samtools sort -o ${POOL_DIR}/${SUBJ_ID}.pooled_reads.final.bam ${POOL_DIR}/${SUBJ_ID}.pooled_reads.fix.bam

# 4. samtools index
samtools index ${POOL_DIR}/${SUBJ_ID}.pooled_reads.final.bam

# 5. cleanup
rm ${POOL_DIR}/${SUBJ_ID}.pooled_reads.sorted.bam
rm ${POOL_DIR}/${SUBJ_ID}.pooled_reads.fix.bam

