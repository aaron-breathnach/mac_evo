#!/usr/bin/env Rscript

library(tidyverse)
library(IRanges)

source("src/get_vcf.R")

"Usage:
   filter_vcf.R [--patients <patients> --metadata <metadata> --gen_len <genome_lengths> --out_dir <out_dir>]

Options:
   --patients List of patients to include [default: data/patients.txt]
   --metadata Genome metadata [default: data/metadata.tsv]
   --gen_len Genome lengths [default: data/genome_lengths.tsv]
   --out_dir Output directory [default: data/]
   
" -> doc

opts <- docopt::docopt(doc)

## this function is taken from https://github.com/gtonkinhill/pneumo_withinhost_manuscript/blob/main/variant_calling.Rmd
filter_regions <- function(minor_vafs, genome_length, expected_null_count=1, p.threshold=0.05) {
  
  #determine null rate of polymorphisms
  nvar <- nrow(minor_vafs)
  null_rate <- nvar/genome_length
  w <- floor(expected_null_count/null_rate)
  
  nvar_adj <- sum(map_dbl(split(minor_vafs, minor_vafs$`#CHROM`), function(chrom_var,i){
    sw_vec <- map_dbl(chrom_var$POS, ~{
      sw <- sum((chrom_var$POS<=(.x+w/2)) & (chrom_var$POS>=(.x-w/2)))
      binom.test(sw, p = null_rate,  n = w, alternative = 'greater')$p.value
    })
    sum(p.adjust(sw_vec, method = 'BH') > p.threshold)
  }))
  
  null_rate <- nvar_adj/genome_length
  w <- max(1e4, floor(expected_null_count/null_rate))
  
  all_sig_vafs <- map_dfr(split(minor_vafs, minor_vafs$`#CHROM`), function(chrom_var,i){
    # Test initial regions
    sw_vec <- map_dbl(chrom_var$POS, ~{
      sw <- sum((chrom_var$POS<=(.x+w/2)) & (chrom_var$POS>=(.x-w/2)))
      binom.test(sw, p = null_rate,  n = w, alternative = 'greater')$p.value
    })
    chrom_var$start <- ceiling(pmax(0, chrom_var$POS - w/2))
    chrom_var$end <- floor(chrom_var$POS + w/2)
    
    # Merge overlapping regions
    sig_vafs <- chrom_var[p.adjust(sw_vec, method = 'BH') < p.threshold, , drop = FALSE] %>% arrange(POS)
    ir <- IRanges(sig_vafs$start, sig_vafs$end)
    sig_vafs$group <- subjectHits(findOverlaps(ir, reduce(ir)))
    regions <- sig_vafs %>% group_by(group) %>%
      summarise(
        start = min(start),
        end = max(end)
      )
    chrom_var$region <- map_int(chrom_var$POS, ~{
      r <- regions$group[(regions$start<=.x) & (regions$end>=.x)]
      if (length(r)<1){
        return(NA)
      } else {
        return(r)
      }
    })
    sig_vafs <- chrom_var %>% filter(!is.na(region))
    
    # refine regions
    sig_vafs <- map_dfr(split(sig_vafs, sig_vafs$region), function(regiondf){
      should_trim <- TRUE
      while (should_trim && (nrow(regiondf) >= 3)) {
        should_trim <- FALSE
        sr <- nrow(regiondf)
        lr <- max(regiondf$POS[[sr]] - regiondf$POS[[1]] + 1, sr)
        p1 <- sr/lr
        llkr <- dbinom(sr, lr, p1, log = TRUE) - dbinom(sr, lr, null_rate, log = TRUE)
        remove_vars <- c()
        
        # attempt trim at start
        st <- sr - 1 
        lt <- max(regiondf$POS[[sr]] - regiondf$POS[[2]] + 1, st)
        pt <- st/lt
        llkt_start <- dbinom(st, lt, pt, log = TRUE) - dbinom(st, lt, null_rate, log = TRUE)
        
        # attempt trim at end
        st <- sr - 1 
        lt <- max(regiondf$POS[[sr - 1]] - regiondf$POS[[1]] + 1, st)
        pt <- st/lt
        llkt_end <- dbinom(st, lt, pt, log = TRUE) - dbinom(st, lt, null_rate, log = TRUE)
        
        if ((llkt_start > llkr) && (llkt_start > llkt_end)) {
          regiondf <- regiondf[-1,,drop = FALSE]
          should_trim <- TRUE
        } else if (llkt_end > llkr) {
          regiondf <- regiondf[-sr,,drop = FALSE]
          should_trim <- TRUE
        }
      }
      
      # Final test to see if region passes threshold
      sr <- nrow(regiondf)
      lr <- max(regiondf$POS[[sr]] - regiondf$POS[[1]] + 1, sr)
      p1 <- sr/lr
      if (binom.test(x = sr, n = lr, p = null_rate, alternative = 'greater')$p.value >= p.threshold/(genome_length/lr)) {
        # region does not pass threshold
        regiondf <- regiondf[0,]
      }
      
      return(regiondf)
    })
    return(sig_vafs)
  })
  
  all_sig_vafs$flag <- NULL
  all_sig_vafs$start <- NULL
  all_sig_vafs$end <- NULL
  minor_vafs$flag <- NULL
  
  keep_vafs <- minor_vafs[!paste(minor_vafs$`#CHROM`, minor_vafs$POS) %in% paste(all_sig_vafs$`#CHROM`, all_sig_vafs$POS),]
  keep_vafs$region <- NA
  return(rbind(keep_vafs, all_sig_vafs))
}

get_genome_size <- function(x) {
  
  lines <- readLines(x)
  
  isolate <- x %>%
    basename() %>%
    str_replace(".txt", "")
  
  bases <- lines[3] %>%
    str_replace(".* ", "") %>%
    as.numeric()
  
  tibble(isolate = isolate, bases = bases)
  
}

.run_filter_regions <- function(vcf, genome_length) {
  
  prefix <- vcf %>%
    basename() %>%
    str_replace("\\..*", "")
  
  minor_vafs <- get_vcf(vcf)
  
  keep <- suppressWarnings(filter_regions(minor_vafs, genome_length)) %>%
    select(1, 2, region) %>%
    filter(is.na(region)) %>%
    mutate(chr_pos = sprintf("%s-%s", `#CHROM`, POS)) %>%
    pull(chr_pos)
  
  vcfR::read.vcfR(vcf, verbose = FALSE) %>%
    slot("fix") %>%
    as_tibble() %>%
    mutate(chr_pos = sprintf("%s-%s", `CHROM`, POS)) %>%
    filter(chr_pos %in% keep) %>%
    mutate(GENOME = prefix) %>%
    select(10, 1:8)
    
}

run_filter_regions <- function(pid, metadata, genome_lengths) {
  
  met <- metadata %>%
    filter(patient == pid)
  
  reference <- met %>%
    filter(time_from_diagnosis == 0) %>%
    pull(isolate)
  
  genome_length <- genome_lengths %>%
    filter(isolate == reference) %>%
    pull(bases)
  
  vcfs <- met %>%
    filter(time_from_diagnosis > 0) %>%
    mutate(vcf = sprintf("data/lofreq/within_host/%s.lofreq.snpeff.vcf.gz", isolate)) %>%
    pull(vcf)
  
  purrr::map(vcfs, function(x) .run_filter_regions(x, genome_length))
  
}

filter_vcf <- function(PATIENTS, METADATA, GEN_LEN, OUT_DIR) {
  
  if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
  
  patients <- readLines(PATIENTS)
  
  metadata <- read_delim(METADATA) %>%
    filter(patient %in% patients)
  
  gen_len <- read_delim(GEN_LEN)
  
  filtered <- patients %>%
    purrr::map(function(x) run_filter_regions(x, metadata, gen_len)) %>%
    bind_rows()
  
  filename <- sprintf("%s/filtered_vcfs.tsv", OUT_DIR)
  
  write_tsv(filtered, filename)
  
}

if (sys.nframe() == 0) {
  filter_vcf(opts$patients, opts$metadata, opts$gen_len, opts$out_dir)
}
