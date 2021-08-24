#!/usr/bin/Rscript

# Calculate scores for distal candidate target genes based on groups of variants
if(dim(targetAnnot)[1] > 0){
  
  # get all SNP-gene combinations
  snp_gene_all <- bind_rows(targetAnnot, eQTLs, ASE) %>%
    dplyr::select(variant, region, signal, gene) %>%
    distinct()
  
  # link SNP to enriched biofeatures
  snp_gene_all_enriched_feature <- left_join(snp_gene_all, enriched_features, by = "variant") %>%
    mutate(enriched_feature = ifelse(!is.na(annotationFile), 1, 0)) %>%
    dplyr::select(variant:region.x, signal.x, gene, enriched_feature) %>%
    distinct() %>%
    rename(region = region.x, signal = signal.x)
  
  # categorise methods
  computationmethods <- c("PreSTIGE", "hnisz", "FANTOM5", "IMPET")
  experimentalmethods <- c("CHIAPET","CHiC","isHiC")
  
  targetAnnot2 <- targetAnnot %>%
    filter(gene != ".", method != ".") %>%
    mutate(method2 = case_when(method %in% computationmethods ~ "compMethod", 
                               method %in% experimentalmethods ~ "expMethod")) %>%
    dplyr::select(variant, region, signal, gene, method2, oeID) %>%
    rename(method = method2) %>%
    distinct()
  
  # count up distinct enhancers for comp and exp -- join this up to main table
  nDistinctEnhancers <- targetAnnot2 %>%
    group_by(variant, region, signal, gene, method) %>%
    summarise(A = n_distinct(oeID)) %>%
    ungroup() %>% 
    spread(method, A, fill = 0) %>%
    rename(N_uniq_compMth_Enh_per_sig = compMethod, N_uniq_expMth_Enh_per_sig = expMethod) 
  
  # get compMethods for TADs
  TAD_variant_gene_check <- left_join(snp_gene_all_enriched_feature, TAD_variant, by = "variant") %>%
    left_join(TAD_tss, by = "gene") %>%
    mutate(TAD_check = (tad.x == tad.y)) %>%
    dplyr::select(variant, region, signal, gene, enriched_feature, TAD_check) %>%
    replace_na(list(TAD_check = 0))
  
  # group by variant...gene, count methods
  distal_targets <- targetAnnot2 %>%
    group_by(variant, region, signal, gene, method) %>%
    count() %>%
    spread(method, n, fill = 0) %>%
    rename(N_compMethod_links = compMethod, N_expMethod_links = expMethod)
  
  # count eQTL combinations
  eQTL_links <- eQTLs %>%
    group_by(variant, region, signal, gene, method) %>%
    count() %>%
    spread(method, n, fill = 0) %>%
    rename(N_eQTL_links = eQTL)
  
  # format ASE data - ASE links only contribute when there hasn't been an eQTL point awarded for the ASE gene (rationale is that the methods are similar, therefore maximum 1 point contribution)
  if (dim(ASE)[1] == 0) {
    ASE2 <- ASE %>% mutate(N_ASE_links = 0)
  } else {
    ASE2 <- anti_join(ASE, eQTLs, by = c("region", "signal", "gene")) %>%
      group_by(region, signal, gene, method) %>%
      count() %>%
      spread(method, n, fill = 0) %>%
      rename(N_ASE_links = ASE)
  }
  
  # merge variant-gene combinations for each method
  join_tables_1 <- left_join(snp_gene_all_enriched_feature, distal_targets, by = c("variant", "region", "signal", "gene")) %>%
    left_join(eQTL_links, by = c("variant", "region", "signal", "gene")) %>%
    left_join(nDistinctEnhancers, by = c("variant", "region", "signal", "gene")) %>%
    left_join(ASE2, by = c("region", "signal", "gene")) %>%
    replace_na(list(N_compMethod_links = 0, N_expMethod_links = 0, N_eQTL_links = 0, N_uniq_compMth_Enh_per_sig = 0, N_uniq_expMth_Enh_per_sig = 0, N_ASE_links = 0))
  
  # merge all tables
  join_tables_2 <- left_join(snp_gene_all_enriched_feature, join_tables_1, by = c("variant", "region", "signal", "gene")) %>%
    left_join(TAD_variant_gene_check, by = c("variant", "region", "signal", "gene")) %>%
    left_join(expression_score, by = "gene") %>%
    replace_na(list(gene_expression = 0)) %>%
    dplyr::select(-enriched_feature.x, -enriched_feature.y)
  
  # convert to score levels:  eg. compMeth 0, expMeth - convert to strata: 0 = 0, 1 = 1, >1 = 2
  DIST_score_step1 <- join_tables_2 %>%
    # scores for prediction methods
    mutate(score_compMeths = ifelse(N_compMethod_links >= 1, 1, 0)) %>%
    mutate(score_expMeths = ifelse(N_expMethod_links >= 1, 2, 0)) %>%
    # bonus for overlap with enriched feature
    mutate(score_compMeths_multi = ifelse((N_compMethod_links >= 1 & N_uniq_compMth_Enh_per_sig >= 1), 1, 0)) %>%
    mutate(score_expMeths_multi = ifelse(N_expMethod_links >= 1 & N_uniq_expMth_Enh_per_sig > 1, 1, 
                                              ifelse(N_expMethod_links >= 1 & N_uniq_expMth_Enh_per_sig == 1, 1, 0))) %>%
    # bonus for overlap with enriched feature
    mutate(score_compMeths_enriched_feature = ifelse((N_compMethod_links >= 1 & enriched_feature >= 1), 1, 0)) %>%
    mutate(score_expMeths_enriched_feature = ifelse((N_expMethod_links >= 1 & enriched_feature >= 1), 1, 0)) %>% 
    # eQTL score
    mutate(score_eQTL = ifelse(N_eQTL_links == 1, 1, 0)) %>%
    # ASE score
    mutate(score_ASE = ifelse(N_ASE_links == 1, 1, 0)) %>%
    # check that computational prediction is within TAD and penalise if not (not applied to experimental links)
    mutate(score_TAD = ifelse((N_eQTL_links > 0 | N_compMethod_links > 0) & N_expMethod_links == 0 & TAD_check == 0, 0.2, 1)) %>%
    dplyr::select(-enriched_feature, -N_compMethod_links, -N_expMethod_links, -N_uniq_compMth_Enh_per_sig, -N_uniq_expMth_Enh_per_sig, -N_eQTL_links, -N_ASE_links, -TAD_check) %>%
    mutate(score_compiled = ((score_compMeths + score_expMeths + score_compMeths_multi + score_expMeths_multi + score_compMeths_enriched_feature + score_expMeths_enriched_feature + score_eQTL + score_ASE) * score_TAD) * gene_expression)
  
  # add driver status
  DIST_driver <- left_join(DIST_score_step1, biotype, by = 'gene') %>%
    left_join(drivers, by = 'gene') %>%
    left_join(SNPinfo, by = c("region","signal")) %>% 
    replace_na(list(driver_status=0)) 
  
  # find max score gene
  DIST_score_summary <- DIST_driver %>%
    dplyr::select(-variant) %>%
    filter(gene != ".") %>%
    distinct() %>%
    group_by(region, signal, gene) %>%
    mutate(gene_score = max(score_compiled) + driver_status ) %>%
    mutate(gene_inquisit_level = case_when(gene_score >4 ~ 1,
                                       gene_score <=4 & gene_score >= 1 ~ 2,
                                       gene_score <1 & gene_score > 0 ~ 3,
                                       gene_score == 0 ~ 0)) %>% 
    arrange(region, signal, desc(gene_score)) %>%
    ungroup() %>% 
    dplyr::select(region,signal,gene,N_SNPs_signal,gene_score,gene_inquisit_level) %>% 
    distinct() 
  
  # find max score gene, include variant
  DIST_score_detail <- DIST_driver %>%
    filter(gene != ".") %>%
    distinct() %>%
    group_by(region, signal, gene) %>%
    mutate(gene_score = max(score_compiled) + driver_status ) %>%
    mutate(gene_inquisit_level = case_when(gene_score >4 ~ 1,
                                           gene_score <=4 & gene_score >= 1 ~ 2,
                                           gene_score <1 & gene_score > 0 ~ 3,
                                           gene_score == 0 ~ 0)) %>% 
    arrange(region, signal, desc(gene_score)) %>%
    ungroup() %>% 
    group_by(region, signal, gene) %>%
    ungroup() %>% 
    dplyr::select(region, signal, variant, gene, biotype, N_SNPs_signal, everything()) 
  
  message("Writing distal target results...")
  
  # write results
  DIST_score_summary %>% 
    write.table(paste0(results_dir_summary, "inquisit.", run, ".summary.distal.tsv"), quote = F, sep = "\t", row.names = F)
  
  # details
  DIST_score_detail %>% 
    write.table(paste0(results_dir_detail, "inquisit.", run, ".detail.distal.tsv"), quote = F, sep = "\t", row.names = F) 
  
} else {
  
  x <- data.frame(matrix(nrow = 1, ncol = 6))
  names(x) <- c("region", "signal", "gene", "N_SNPs_signal", "gene_score", "gene_inquisit_level")
  x %>% write.table(paste0(results_dir_summary, "inquisit.", run, ".summary.distal.tsv"), quote = F, sep = "\t", row.names = F)
  y <- data.frame(matrix(nrow = 1, ncol = 6))
  names(y) <- c("region", "signal", "gene", "N_SNPs_signal", "gene_score", "gene_inquisit_level")
  y %>% write.table(paste0(results_dir_detail, "inquisit.", run, ".detail.distal.tsv"), quote = F, sep = "\t", row.names = F)
  
}

