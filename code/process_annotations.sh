#!/bin/bash

# Process non-interval based data

echo
echo "Processing annotation data:"

echo
echo "...Variants..."

# SNPs
echo -e "varChr\tvarStart\tvarEnd\tvariant\tsnp_icogs\trsID\tregion\tsignal\tphenotype" |
cat - ${SNPs} > ${TEMPDIR}/variants.txt

echo "...Allele-specific expression..."

# process ASE data
cut -f1-3 ${DATADIR}/ase/ase_data.txt > ${TEMPDIR}/ase_processed.tmp

echo -e "gene\tregion\tsignal" |
cat - ${TEMPDIR}/ase_processed.tmp > ${TEMPDIR}/ase_processed.txt

# echo "...Capture Hi-C..."
# 
# # process CHi-C data
# cp ${DATADIR}/CHiC/CHiC_data.txt ${TEMPDIR}/CHiC_processed.tmp
# 
# echo -e "oeID\tpromoterFrag\tgene\tbiotype\tvariant\tregion\tsignal" |
# cat - ${TEMPDIR}/CHiC_processed.tmp > ${TEMPDIR}/CHiC_processed.txt

echo "...Coding variants..."

# process coding data
for i in ${DATADIR}/coding/*.txt ; do

variant_type=$( basename ${i%.txt} )
awk -F'\t' -v A=${variant_type} '{ print $0"\t""CODING""\t"A }' ${i} > ${TEMPDIR}/coding_${variant_type}.tmp
done

echo -e "variant\tregion\tsignal\tgene\tmethod\ttype" |
cat - ${TEMPDIR}/coding_*.tmp > ${TEMPDIR}/coding_variants.txt

echo "...eQTLs..."

# process eQTL data
cp ${DATADIR}/eQTL/eQTL_data.txt ${TEMPDIR}/eQTL_data.tmp

echo -e "variant\tregion\tsignal\tgene" |
cat - ${TEMPDIR}/eQTL_data.tmp > ${TEMPDIR}/eQTL_data.txt

echo "...Gene expression..."

# process expression data
cp ${DATADIR}/expression/expression.txt ${TEMPDIR}/expression.tmp

echo -e "gene\tvalue" |
cat - ${TEMPDIR}/expression.tmp > ${TEMPDIR}/expression.txt

echo "...Post-translational modifications..."

# process post-translational modification data
cp ${DATADIR}/ptm/ptm_data.txt ${TEMPDIR}/ptm_data.tmp

echo -e "variant\tregion\tsignal\tphenotype\tgene\tptm_cat\tptm_score" |
cat - ${TEMPDIR}/ptm_data.tmp > ${TEMPDIR}/ptm_data.txt

# driver genes
echo -e "gene" | cat - ${DATADIR}/driver_genes/bc.drivers.txt > ${TEMPDIR}/bc_drivers.txt

# gene annotations
echo -e "gene\tbiotype" | cat - ${DATADIR}/gene_annotations/gencode.v19.biotype.txt > ${TEMPDIR}/gencode_v19_biotype.txt

# remove temp files
rm -rf ${TEMPDIR}/*.tmp

echo
echo "Done."

