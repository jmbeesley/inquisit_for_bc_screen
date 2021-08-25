#!/bin/bash

# run with bash INQUISIT.sh <name_of_run> <variant_bed_file> [ e.g. bash inquisit.sh "2021-08-22" BCACFM.CCVs.bed ]

# assign variables
RUN="$1"
SNPs="$2"
INQDIR_RUN=output/${RUN}
DATADIR=reference/breast_annotations/
TEMPDIR=${INQDIR_RUN}/temp

# HPC modules
module load bedtools/2.26.0
module load R/3.4.1

echo
echo "Running INQUISIT for $( basename ${SNPs} ) SNP set using data in $( basename ${DATADIR} ) for ${RUN}"

# set up directories for run: Each run requires VARIANTS and ANNOTATIONS which might differ, therefore make a standard directory structure for each run
mkdir -p ${INQDIR_RUN}/results/{detail,summary} ${TEMPDIR}

# process non-interval annotation data
source code/process_annotations.sh

# intersect SNPs with annotations
source code/intersect_annotations.sh

# scores
echo "Calculating scores..."
echo

# calculate gene scores
Rscript code/compute_scores.R ${INQDIR_RUN} ${RUN}

# list Level 1 genes
cut -f1-3,6 ${INQDIR_RUN}/results/summary/inquisit.${RUN}.summary.* | 
sort |
uniq |
awk -F'\t' 'BEGIN{ print "region""\t""signal""\t""gene" } $NF == "1" { print $1"\t"$2"\t"$3 }' > ${INQDIR_RUN}/results/inquisit_level1_genes.txt

echo
echo "Finished."
