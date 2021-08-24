#!/usr/bin/Rscript

# Import and process data

# load packages
suppressMessages(library(dplyr, warn.conflicts = FALSE))
suppressMessages(library(reshape2))
suppressMessages(library(tidyr))

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

# processed data
datadir <- paste0(run_dir, "/temp/")

# results destinations
results_dir_summary <- paste0(run_dir, "/results/summary/")
results_dir_detail <- paste0(run_dir, "/results/detail/")

# functions
import_data <- function(x){
  a <- read.delim(x, sep = "\t", stringsAsFactors = F)
  return(a)
}

message("Importing and formatting data...")
message("")

# import data
processed_data_files <- paste0(datadir, dir(datadir))
processed_data_list <- list()
processed_data_list <- lapply(processed_data_files, import_data)
names(processed_data_list) <- gsub("\\.txt", "", dir(datadir))

# variants
variants <- processed_data_list[["variants"]] 

# locus and SNP info
SNPinfo <- variants %>%
  group_by(region, signal) %>%
  summarise(N_SNPs_signal = n())

# format gene link data 
targetAnnot <- processed_data_list[["target_gene_links_intersect"]] %>% 
  separate(annotation, into = c("celltype", "annotInfo", "method"), sep = "-", fill = "left")

# format eQTL data
eQTLs <- processed_data_list[["eQTL_data"]] %>% 
  dplyr::select("variant", "region", "signal", "gene") %>% 
  mutate(method = "eQTL") %>% 
  distinct() %>%
  filter(variant %in% variants$variant)

# add row
if (dim(eQTLs)[1] == 0) {
  eQTLs[nrow(eQTLs) + 1, ] <- 0
  eQTLs <- eQTLs %>% mutate(method = "eQTL")
}

# get relevant ASE genes
ASE <- inner_join(processed_data_list[["ase_processed"]], variants, by = c("region", "signal")) %>%
  dplyr::select("gene", "region", "signal") %>% mutate(method = "ASE") %>% 
  distinct()

# add row
if (dim(ASE)[1] == 0) {
  ASE[nrow(ASE) + 1, ] <- 0
}

# format TAD data
TAD_variant <- processed_data_list[["variant_tad_intersect"]] %>% 
  dplyr::select("variant", "tad") %>% distinct()
TAD_tss <- processed_data_list[["tss_tad_intersect"]] %>% 
  dplyr::select("gene", "tad") %>% distinct()

# gene expression
expression_score <- processed_data_list[["expression"]] %>%
  mutate(gene_expression = ifelse(value > 0, 1, 0.1)) %>%
  dplyr::select(-value) %>%
  distinct()

enriched_features <- processed_data_list[["enriched_features_intersect"]] 
promoterData <- processed_data_list[["promoter_intersect"]]
codingData <- processed_data_list[["coding_variants"]] 
ptmData <- processed_data_list[["ptm_data"]] 
drivers <- processed_data_list[["bc_drivers"]] %>% mutate(driver_status = 1)
biotype <- processed_data_list[["gencode_v19_biotype"]]

message("Done.")
message("")


