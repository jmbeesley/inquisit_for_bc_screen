# INQUISIT GWAS target gene prioritisation

## Background

Candidate breast cancer risk genes included in the CRISPR screens in Tuano et al., were selected primarily using the INQUISIT heuristic.

**INQUISIT** (*Integrated eQTL and In Silico Prediction of GWAS Targets*) was originally developed to prioritise candidate target genes underlying breast cancer risk signals detected by GWAS. This was presented in Michailidou et al., Nature, 2017 and subsequently in an updated version to prioritise candidate target genes at fine-mapped breast cancer risk loci (Fachal et al., Nature Genetics, 2020). The latter (now out-of-date) version is presented in this workflow. 

INQUISIT is a composite method which hypothesises that trait-associated SNPs impact effector genes by affecting coding potential (protein sequence or through splicing), through modulaing the promoter, or through distal gene regulation. The method combines evidence from a range of functional genomics assays and computational approaches that variously aim to link regulatory elements such as promoters and enhancers with target genes. The hypothesis is that candidate variants within a GWAS locus (e.g. correlated with the index SNP) will predominantly act through one or more of these categories on one or more target genes.

## Data

  - **Gene annotations**: GENCODE v19 basic
  - **Promoter regions** - GENCODE v19 basic transcription start sites, upstream 1kb, and downstream 100bp.
  - **Coding variants** - Fine-mapped breast cancer risk variants pre-processed with VEP and alamut, 
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
          - Beesley et al., Genome Biol 2020 (Capture Hi-C), 
          - Rao et al., Cell 2014 (HiC); 
      - Computational correlative methods using epigenomic data from breast cells: 
          - PreSTIGE, Hnisz, FANTOM5, IM-PET
  - **Breast cancer driver genes** - 

## Running INQUISIT

This code simply copies the pre-processed annotation data, performs intersections between input variants and these annotations, then computes scores for each category of potential SNP effects (coding, promoter, distal). A bed file containing BCAC candidate casual variants from the fine-mapping analysis (Fachal et al., 2019) is provided in the reference directory.

```
bash code/inquisit.sh <run_name> <snp_list>
```

## Interpreting results

Genes prioritised in each of the mechanistic categories are ranked by these scores which reflect the degree of supporting evidence from all the candidate causal variants contained within the association signal. Because the score ranges differ between categories, scores are converted to levels (1-3), where level 1 genes are supported by the strongest evidence, generally from multiple sources or methods. Level 2 genes should be viewed as possibile targets, but demand less focus. Level 3 are very unlikely to be true targets.

## References

[Michailidou et al., Nature, 2017](https://www.nature.com/articles/nature24284)
[Fachal et al., Nature Genetics 2020](https://www.nature.com/articles/s41588-019-0537-1)
[Beesley et al., Genome Biol 2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1877-y)
[Rao et al., Cell 2014](https://www.cell.com/fulltext/S0092-8674(14)01497-4)
