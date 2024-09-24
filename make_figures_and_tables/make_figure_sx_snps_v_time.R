library(tidyverse)

metadata <- read_delim("data/metadata.tsv") %>%
  filter(!multiple_carriage)

dat <- read_delim("data/minor_vafs/filtered_vcfs.txt") %>%
  filter(GENOME %in% metadata$isolate)

get_af <- function(x) {
  x %>%
    str_replace(".*AF=", "") %>%
    str_replace("\\;.*", "") %>%
    as.numeric()
}

parse_info <- function(x, n) {
  x %>%
    str_split("\\|") %>%
    unlist() %>%
    nth(n)
}

dat$allele_frequency <- purrr::map(dat$INFO, function(x) get_af(x)) %>%
  unlist()

dat$gene <- purrr::map(dat$INFO, function(x) parse_info(x, 4)) %>%
  unlist()

dat$variant <- purrr::map(dat$INFO, function(x) parse_info(x, 2)) %>%
  unlist()

df <- dat %>%
  mutate(snp = paste0(CHROM, "|", POS)) %>%
  select(GENOME, gene, snp, allele_frequency, variant) %>%
  filter(variant %in% "missense_variant") %>%
  rename(isolate = 1) %>%
  inner_join(metadata, by = "isolate") %>%
  filter(allele_frequency > 0.95)

count_numb_snps_per_timepoint <- function(pid, metadata, dat) {
  
  meta <- metadata %>%
    filter(patient == pid & time_from_diagnosis > 0) %>%
    arrange(time_from_diagnosis)
  
  snps <- c()
  DF <- c()
  
  for (i in 1:nrow(meta)) {
    
    ISOLATE <- meta[[i, "isolate"]]
    
    snp_ids <- df %>%
      filter(isolate == ISOLATE) %>%
      pull(snp)
    
    n <- length(snp_ids)
    
    for (snp_id in snp_ids) {
      if (!snp_id %in% snps) {
        snps <- c(snps, snp_id)
      }
    }
    
    tmp <- tibble(isolate = ISOLATE, numb_snps = n)
    DF <- rbind(DF, tmp)
    
  }
  
  return(DF)
  
}

patients <- metadata %>%
  filter(time_from_diagnosis > 0) %>%
  group_by(patient) %>%
  tally() %>%
  ungroup() %>%
  filter(n > 3) %>%
  pull(patient)

res <- purrr::map(patients, function(x) count_numb_snps_per_timepoint(x, metadata, dat)) %>%
  bind_rows() %>%
  inner_join(metadata, by = "isolate") %>%
  select(numb_snps, time_from_diagnosis, patient) %>%
  group_by(patient) %>%
  mutate(n = sum(numb_snps)) %>%
  ungroup() %>%
  filter(n > 0)

ggplot(res, aes(x = time_from_diagnosis, y = numb_snps)) +
  facet_wrap(~patient, scales = "free") +
  geom_point() +
  geom_smooth(method = "lm")
