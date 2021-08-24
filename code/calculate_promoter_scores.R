#!/usr/bin/Rscript

# compile and calculate promoter based scores
if (dim(promoterData)[1] > 0) {
  
  # filter for variants of interest
  promoterData <- promoterData %>%
    filter(variant %in% variants$variant)
  
  # Add eQTL and ASE data
  promoter_eQTL_ASE <- left_join(promoterData, eQTLs, by = c("variant", "pGene" = "gene")) %>%
    left_join(ASE2, by = c("region.x" = "region", "signal.x" = "signal", "pGene" = "gene")) %>%
    mutate(score_eQTL = ifelse(method == "eQTL" | N_ASE_links == 1, 1, 0)) %>%
    dplyr::select(-region.y, -signal.y, -method, -N_ASE_links) %>%
    distinct() %>%
    replace_na(list(score_eQTL = 0))
  
  # link to biofeatures
  promoter_enriched_feature <- left_join(promoter_eQTL_ASE, enriched_features, by = "variant") %>%
    mutate(enriched_feature = ifelse(!is.na(annotationFile), 1, 0)) %>%
    dplyr::select(pGene, region.x, signal.x, variant, score_eQTL, enriched_feature) %>%
    distinct()
  
  # annotate proximal TSS
  promoter_enrF_eQTL <- promoter_enriched_feature %>%
    mutate(proximal = 1) %>%
    rename(gene = pGene, region = region.x, signal = signal.x)
  
  # link to expression table
  promoter_enrF_eQTL_Exp <- left_join(promoter_enrF_eQTL, expression_score, by = "gene") %>% 
    replace_na(list(gene_expression = 0))
  
  # calculate scores
  PROM_score_step1 <- promoter_enrF_eQTL_Exp %>%
    mutate(score_compiled = (score_eQTL + enriched_feature + proximal) * gene_expression)
  
  # add driver status
  PROM_driver <- left_join(PROM_score_step1, biotype, by = "gene") %>%
    left_join(drivers, by = "gene") %>%
    left_join(SNPinfo, by = c("region", "signal")) %>%
    replace_na(list(driver_status = 0))
  
  # find max score gene
  PROM_score_summary <- PROM_driver %>%
    dplyr::select(-variant) %>%
    filter(gene != ".") %>%
    distinct() %>%
    group_by(region, signal, gene) %>%
    mutate(gene_score = max(score_compiled) + driver_status) %>%
    mutate(gene_inquisit_level = case_when(
      gene_score >= 3 ~ 1,
      gene_score == 1 | gene_score == 2 ~ 2,
      gene_score < 1 & gene_score > 0 ~ 3,
      gene_score == 0 ~ 0
    )) %>%
    arrange(region, signal, desc(gene_score)) %>%
    ungroup() %>% 
    dplyr::select(region,signal,gene,N_SNPs_signal,gene_score,gene_inquisit_level) %>% 
    distinct() 
  
  # find max score gene, include variant
  PROM_score_detail <- PROM_driver %>%
    filter(gene != ".") %>%
    distinct() %>%
    group_by(region, signal, gene) %>%
    mutate(gene_score = max(score_compiled) + driver_status) %>%
    mutate(gene_inquisit_level = case_when(
      gene_score >= 3 ~ 1,
      gene_score == 1 | gene_score == 2 ~ 2,
      gene_score < 1 & gene_score > 0 ~ 3,
      gene_score == 0 ~ 0
    )) %>%
    arrange(region, signal, desc(gene_score)) %>%
    ungroup() %>% 
    group_by(region, signal, gene) %>%
    ungroup() %>% 
    dplyr::select(region, signal, variant, gene, biotype, N_SNPs_signal, gene_expression, everything()) 
  
  message("Writing promoter target results...")
  
  # write results
  PROM_score_summary %>% 
    write.table(paste0(results_dir_summary, "inquisit.", run, ".summary.promoter.tsv"), quote = F, sep = "\t", row.names = F)
  
  # details
  PROM_score_detail %>% 
    write.table(paste0(results_dir_detail, "inquisit.", run, ".detail.promoter.tsv"), quote = F, sep = "\t", row.names = F) 
  
} else {
  x <- data.frame(matrix(nrow = 1, ncol = 6))
  names(x) <- c("region", "signal", "gene", "N_SNPs_signal", "gene_score", "gene_inquisit_level")
  x %>% write.table(paste0(results_dir_summary, "inquisit.", run, ".summary.promoter.tsv"), quote = F, sep = "\t", row.names = F)
  y <- data.frame(matrix(nrow = 1, ncol = 6))
  names(y) <- c("region", "signal", "gene", "N_SNPs_signal", "gene_score", "gene_inquisit_level")
  y %>% write.table(paste0(results_dir_detail, "inquisit.", run, ".detail.promoter.tsv"), quote = F, sep = "\t", row.names = F)
}
