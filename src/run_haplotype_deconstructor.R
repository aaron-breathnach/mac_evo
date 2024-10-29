#!/usr/bin/env Rscript

library(HaplotypeDeconstructor)
library(NMF)
library(tidyverse)

source("src/get_vcf.R")

.run_haplotype_deconstructor <- function(mat, threshold) {
  
  mat <- mat[rowSums(mat) != 0,]
  
  set.seed(123)
  mat <- mat[sample.int(nrow(mat), 1000),]
  mat <- mat[, colSums(mat) > 0]
  
  res <- assessNumberHaplotyes(mat, 2:10, nReplicates = 3)
  
  exp_v_hap <- tibble(num_hap = res$NumberHaplotyes,
                      explained = 100*res$ExplainedVariance) %>%
    distinct() %>%
    group_by(num_hap) %>%
    summarise(explained = mean(explained)) %>%
    ungroup()
  
  number_of_haplotypes <- exp_v_hap %>%
    filter(explained >= threshold) %>%
    filter(explained == min(explained)) %>%
    pull(num_hap)
  
  decomposed <- findHaplotypes(mat, number_of_haplotypes)
  
  haplotypes <- decomposed$samples %>%
    as.data.frame() %>%
    rownames_to_column("isolate") %>%
    pivot_longer(!isolate, names_to = "haplotype", values_to = "contribution") %>%
    mutate(haplotype = haplotype %>%
             str_replace("H", "") %>%
             str_pad(2, "left", 0)) %>%
    mutate(haplotype = paste0("H", haplotype)) %>%
    group_by(isolate) %>%
    mutate(rel_contribution = 100 * contribution / sum(contribution)) %>%
    ungroup()
  
  list(exp_v_hap, haplotypes)
  
}

run_haplotype_deconstructor <- function(threshold = 80, present_study = FALSE) {
  
  if (present_study) {
    rds <- "data/haplotype_deconstructor.irish_cohort.RDS"
  } else {
    rds <- "data/haplotype_deconstructor.RDS"
  }
  
  if (!file.exists(rds)) {
    
    files <- list.files("data/lofreq/common_reference", full.names = TRUE)
    
    # patients <- readLines("data/reinfected_patients.txt")
    
    metadata <- read_delim("data/metadata.tsv") %>%
      filter(bracken_pass)
    
    if (present_study) {
      
      metadata <- metadata %>%
        filter(study == "Present")
      
    }
    
    vcfs <- metadata %>%
      filter(bracken_pass) %>%
      mutate(vcf = sprintf("data/lofreq/common_reference/%s.lofreq.vcf.gz", isolate)) %>%
      pull(vcf)
    
    dat <- purrr::map(vcfs, function(x) {
      get_vcf(x) %>%
        dplyr::rename(CHROM = 1) %>%
        mutate(CHROM = str_replace(CHROM, ".*\\|", "")) %>%
        mutate(GENOME = basename(x) %>% str_replace("\\..*", "")) %>%
        mutate(SNP = sprintf("%s.%s.%s_%s", CHROM, POS, REF, ALT)) %>%
        mutate(AF = INFO %>%
                 str_replace(".*AF=", "") %>%
                 str_replace(".SB=.*", "") %>%
                 as.numeric()) %>%
        select(SNP, GENOME, AF)
    }) %>%
      bind_rows() %>%
      pivot_wider(names_from = "GENOME", values_from = "AF", values_fill = 0)
    
    mat <- dat %>%
      column_to_rownames("SNP") %>%
      as.matrix()
    
    res <- .run_haplotype_deconstructor(MAT, threshold)
    
    saveRDS(res, rds)
    
  }
}

if (sys.nframe() == 0) {
  run_haplotype_deconstructor()
}
