#!/usr/bin/Rscript

# Calculate scores for each candidate target gene based on groups of variants

# load packages
suppressMessages(library(dplyr, warn.conflicts = FALSE))
suppressMessages(library(reshape2))
suppressMessages(library(tidyr))

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

# args
args <- commandArgs(TRUE)
args[1] -> run_dir
args[2] -> run

# pre-process
source("code/import_annotations.R")

# compute distal scores
source("code/calculate_distal_scores.R")
# compute promoter scores
source("code/calculate_promoter_scores.R")
# compute coding scores
source("code/calculate_coding_scores.R")

