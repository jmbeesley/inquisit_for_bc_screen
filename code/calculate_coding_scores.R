#!/usr/bin/Rscript

# coding variants
if (dim(codingData)[1] > 0) {

  # import
  codingData <- codingData %>%
    filter(variant %in% variants$variant)

  # PTM data
  ptmData <- ptmData %>%
    filter(variant %in% variants$variant) %>%
    dplyr::select(variant:signal, gene) %>%
    distinct() %>%
    mutate(method = "CODING", type = "PTM")

  # PTM variants only add a point if variant is not counted in missense set
  ptmData_novel <- anti_join(ptmData, codingData, by = c("region", "signal", "gene"))

  # concat all coding changes
  coding_2 <- bind_rows(codingData, ptmData_novel)

  #
  coding_count <- coding_2 %>%
    dplyr::select(region, signal, gene, type) %>%
    distinct() %>%
    filter(!is.na(gene)) %>%
    mutate(score = 1) %>%
    spread(type, score, fill = 0)

  # link to gene expression
  codingScore_preExp <- left_join(coding_count, expression_score, by = "gene")

  # calculate scores
  CODING_score_step1 <- codingScore_preExp %>%
    mutate(score_compiled = (MISSENSE + PTM + SPLICING + STOP) * gene_expression)

  # add driver status
  COD_driver <- left_join(CODING_score_step1, biotype, by = "gene") %>%
    left_join(drivers, by = "gene") %>%
    left_join(SNPinfo, by = c("region", "signal")) %>%
    replace_na(list(driver_status = 0))

  # find max score gene
  COD_score_summary <- COD_driver %>%
    filter(gene != ".") %>%
    distinct() %>%
    group_by(region, signal, gene) %>%
    mutate(gene_score = max(score_compiled) + driver_status) %>%
    mutate(gene_inquisit_level = case_when(
      gene_score > 1 ~ 1,
      gene_score == 1 ~ 2,
      gene_score < 1 & gene_score > 0 ~ 3,
      gene_score == 0 ~ 0
    )) %>%
    arrange(region, signal, desc(gene_score)) %>%
    ungroup() %>% 
    dplyr::select(region,signal,gene,N_SNPs_signal,gene_score,gene_inquisit_level) %>% 
    distinct() 
  
  # coding score is slightly different to distal and promoter 
  COD_score_detail <- coding_2 %>%
    dplyr::select(region, signal, variant, gene, type) %>%
    distinct() %>%
    filter(!is.na(gene)) %>%
    mutate(score = 1) %>%
    spread(type, score, fill = 0) %>% 
    left_join(COD_score_summary, by = c("region", "signal","gene")) %>% 
    left_join(expression_score, by = "gene") %>% 
    left_join(biotype, by = "gene") %>%
    left_join(drivers, by = "gene") %>% 
    replace_na(list(driver_status = 0)) %>% 
    dplyr::select(region, signal, variant, gene, biotype, N_SNPs_signal, gene_expression, MISSENSE:STOP, driver_status, gene_score, gene_inquisit_level) 

  message("Writing coding target results...")

  # write results
  COD_score_summary %>%
    write.table(paste0(results_dir_summary, "inquisit.", run, ".summary.coding.tsv"), quote = F, sep = "\t", row.names = F)
  
  # details
  COD_score_detail %>%
    write.table(paste0(results_dir_detail, "inquisit.", run, ".detail.coding.tsv"), quote = F, sep = "\t", row.names = F)

  } else {
  
  x <- data.frame(matrix(nrow = 1, ncol = 6))
  names(x) <- c("region", "signal", "gene", "N_SNPs_signal", "gene_score", "gene_inquisit_level")
  x %>% write.table(paste0(results_dir_summary, "inquisit.", run, ".summary.coding.tsv"), quote = F, sep = "\t", row.names = F)
  y <- data.frame(matrix(nrow = 1, ncol = 6))
  names(y) <- c("region", "signal", "gene", "gene_expression", "gene_score", "gene_inquisit_level")
  y %>% write.table(paste0(results_dir_detail, "inquisit.", run, ".detail.coding.tsv"), quote = F, sep = "\t", row.names = F)

}
