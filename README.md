# inquisit_bc_screen_prototype

## Background

Candidate breast cancer risk genes analysed in the CRISPR screens in XXXX were selected primarily using the INQUISIT heuristic.

**INQUISIT** (*Integrated eQTL and In Silico Prediction of GWAS Targets*) was originally developed to prioritise candidate target genes underlying breast cancer risk signals detected by GWAS. This was presented in Michailidou et al., Nature, 2017 and subsequently in an updated version to prioritise candidate target genes at fine-mapped breast cancer risk loci (Fachal et al., Nature Genetics, 2019). This (now out-of-date) is the version presented in this workflow. 

It is a composite method which considers SNPs act on effector genes by affecting a gene's coding potential (protein sequence or through splicing), via modulaing the promoter, or through distal gene regulation. The method combines evidence from a range of experimental functional genomics assays and computational approach that variously aim to link regulatory elements such as enhancers with target genes. The hypothesis is that candidate variants within a GWAS locus (e.g. correlated with the index SNP) will predominantly act through one or more of these categories on one or more target genes.

## Data

  - **Gene annotations**: GENCODE v19 basic
  - **Promoter regions** - GENCODE v19 basic transcription start sites, upstream 1kb, and downstream 100bp.
  - **Coding variants** - Fine-mapped breast cancer risk variants pre-processed with VEP, alamut, 
  - **Allele-specific expression** - [obsolete] Allelic-imbalance analysis, probably erroneous
  - **Enriched biofetures** - Epigenomic annotations enriched with BCAC FM CCVs
  - **eQTL** - eQTLs from breast cancer studies (METABRIC, TCGA, NHS), filtered for study wide P < 1e-4, then for p value ratio with best CCV < 100.
  - **Expression** - average gene expression from a range a breast cell lines
  - **Promoter** - genomic annotations from Roadmap ChromHMM analysis from breast samples. Intervals with TSS characteristic signatures
  - **Post translational modifications** - Mutations in protein PTM sites. From the AWESOMEdb.
  - **Topological-associating domains** - ENCODE T47D TAD domains
  - **Target gene links**: 
      - Experimental chromatin interaction data from breast samples: 
          - ENCODE (ChIA-PET), 
          - Beesley et al., Genome Biol 2019 (Capture Hi-C), 
          - Rao et al., Cell 2013 (HiC); 
      - Computation correlative methods using epigenomic data from breast cells: 
          - PreSTIGE, Hnisz, FANTOM5, IM-PET
  - **Breast cancer driver genes** - 

## Running INQUISIT

This code simply copies the pre-porcessed annotation data, performs intersections between input variants and these annoations, then computes scores for each category of potential SNP effect (coding, promoter, distal).

bash inquisit.sh <run_name> snp_list.bed

## Interpreting results

Genes prioritised in each of the mechanistic categories are ranked by these score, where the score reflects the degree of spporting evidence from all the candidate causal SNPs in the signal. Scores are converted to levels (1-3), where level 1 genes have the most evidence of potentially functional SNPs acting upon them. Level 2 genes shoould be viewed as possibile targets, but demand less focus. Level 3 are very unlikely targets.


