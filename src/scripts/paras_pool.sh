#!/bin/bash

#OUTPUT_DIR='/media/GenomeData/parascopy-sabyun/Projects/Parascopy-ExomeDepth/paras-pool/pooled-reads-bam-nomo1'
#HOMOLOGY_DIR='/home/sabyun/parascopy_test/data/homology_table'
#REF_FILEPATH='/media/SHARED1/ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna'

# Parse arguments
FILE_NAME=$1
SRR_NAME=$(echo ${FILE_NAME} | awk -F '::' '{print $2}')
LOCI_NAME=$2
LOCI_REGION=$3
HOM_FP=$4
REF_FP=$5
POOL_DIR=pooled-reads-bam-${LOCI_NAME}

# Run Parascopy pool
echo -e "\nRunning parascopy pool for ${SRR_NAME}.\n"

parascopy pool \
-i ${FILE_NAME} \
-t ${HOM_FP} \
-f ${REF_FP} \
-r ${LOCI_REGION} \
-o ${POOL_DIR}/${SRR_NAME}.pooled_reads.bam

