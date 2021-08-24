#!/bin/bash

# Intersects SNPs with enriched feature annotations and genomic interval annotation data

echo
echo "Intersecting SNPs:"
echo
echo "...in enriched biofeatures..."

# general
TSS=data/breast_annotations/gene_annotations/gencode.v19.basic.tss.1.1kb.bed.gz
VARIANTS=${SNPs}

# intersect enriched feature annotations
for j in ${DATADIR}/enriched_features/annotation_files/*.bed.gz ; do

    ANNOTATION=$( basename ${j} )
   
    bedtools intersect -a ${VARIANTS} -b ${j} -loj |
    awk -F'\t' -v A=${ANNOTATION} '$10 != "." { print $0"\t"A }'

done > ${TEMPDIR}/enriched_features_intersect.tmp

# add header
echo -e "varChr\tvarStart\tvarEnd\tvariant\tsnp_icogs\trsID\tregion\tsignal\tphenotype\tfeatChr\tfeatStart\tfeatEnd\tannotationFile" |
cat - ${TEMPDIR}/enriched_features_intersect.tmp > ${TEMPDIR}/enriched_features_intersect.txt

echo "...in annotated enhancers and interactions..."

# intersect target prediction annotations
for j in ${DATADIR}/target_gene_links/annotation_files/*.bed.gz ; do

    ANNOTATION=$( basename ${j} )

    bedtools intersect -a ${VARIANTS} -b ${j} -loj |
    awk -F'\t' -v A=${ANNOTATION} ' { print $0"\t"A }'

done > ${TEMPDIR}/target_gene_links_intersect.tmp

# add header
echo -e "varChr\tvarStart\tvarEnd\tvariant\tsnp_icogs\trsID\tregion\tsignal\tphenotype\toeChr\toeStart\toeEnd\tgene\tannotation\toeID\tfile" |
cat - ${TEMPDIR}/target_gene_links_intersect.tmp > ${TEMPDIR}/target_gene_links_intersect.txt

# intersect promoter annotations
echo "...in annotated promoters..."

# intersect promoter annotations
for j in ${DATADIR}/promoter/annotation_files/*.bed.gz ; do

    ANNOTATION=$( basename ${j} )

    bedtools intersect -a ${VARIANTS} -b ${TSS} -wa -wb |
    bedtools intersect -a - -b ${j} -loj |
    awk -F'\t' -v A=${ANNOTATION} ' { print $0"\t"A }'

done > ${TEMPDIR}/promoter_intersect.tmp

# add header
echo -e "varChr\tvarStart\tvarEnd\tvariant\tsnp_icogs\trsID\tregion\tsignal\tphenotype\tpChr\tpStart\tpEnd\tpGene\tp2chr\tp2start\tp2end\tpAnnotation\tfile" |
cat - ${TEMPDIR}/promoter_intersect.tmp > ${TEMPDIR}/promoter_intersect.txt

# TADs
echo "...variant and TSS TAD intersection..."

# intersect TAD annotations
for j in ${DATADIR}/TAD/*.bed ; do

    ANNOTATION=$( basename ${j} )

    bedtools intersect -a ${VARIANTS} -b ${j} -wa -wb > ${TEMPDIR}/variant_tad_intersect.tmp
    bedtools intersect -a ${TSS} -b ${j} -wa -wb > ${TEMPDIR}/tss_tad_intersect.tmp

done 

# add header
echo -e "varChr\tvarStart\tvarEnd\tvariant\tsnp_icogs\trsID\tregion\tsignal\tphenotype\ttChr\ttStart\ttEnd\ttad" |
cat - ${TEMPDIR}/variant_tad_intersect.tmp > ${TEMPDIR}/variant_tad_intersect.txt

# add header
echo -e "tssChr\ttssStart\ttssEnd\tgene\ttChr\ttStart\ttEnd\ttad" |
cat - ${TEMPDIR}/tss_tad_intersect.tmp > ${TEMPDIR}/tss_tad_intersect.txt

# remove temp files
rm -rf ${TEMPDIR}/*.tmp

echo
echo "Done."
echo
