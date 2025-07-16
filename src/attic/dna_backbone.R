library(tidyverse)

sort_pids <- function(x, y) {
  values <- sort(c(x, y))
  paste0(values[1], "___", values[2])
}

list_acquired_mutations <- function(pid, metadata, variants) {
  
  final_isolate <- metadata %>%
    filter(patient == pid) %>%
    pull(isolate)
  
  variants %>%
    filter(GENOME %in% final_isolate) %>%
    pull(gene_name) %>%
    unique()
  
}

metadata <- read_delim("data/metadata.tsv")

patients <- metadata %>%
  filter(!multiple_carriage) %>%
  pull(patient) %>%
  unique()

cols <- c("pair", "patient", "genome", "time_from_diagnosis")

snp_dists <- read_delim("data/snp_dists.tsv")

meta <- read_delim("data/metadata.tsv") %>%
  select(isolate, patient, country, time_from_diagnosis) %>%
  dplyr::rename(name = 1)

dat <- snp_dists %>%
  dplyr::rename(genome_1 = 1) %>%
  filter(genome_1 != "Reference") %>%
  select(-any_of("Reference")) %>%
  pivot_longer(!genome_1, names_to = "genome_2", values_to = "dist") %>%
  inner_join(meta[,c(1, 2, 4)], by = c("genome_1" = "name")) %>%
  inner_join(meta[,c(1, 2, 4)], by = c("genome_2" = "name")) %>%
  filter(patient.x %in% patients & patient.y %in% patients) %>%
  filter(patient.x != patient.y) %>%
  filter(dist <= 12) %>%
  mutate(row_num = row_number()) %>%
  group_by(row_num) %>%
  mutate(pair = sort_pids(patient.x, patient.y)) %>%
  ungroup() %>%
  select(pair, patient.x, genome_1, time_from_diagnosis.x) %>%
  distinct() %>%
  group_by(patient.x, pair) %>%
  filter(time_from_diagnosis.x == min(time_from_diagnosis.x)) %>%
  ungroup() %>%
  arrange(pair) %>%
  setNames(cols)

variants <- read_delim("https://raw.githubusercontent.com/aaron-breathnach/mac_evo/refs/heads/main/data/filtered_variants.tsv")
gene_data <- read_delim("data/panaroo/gene_data.csv")

variants$annotation_id <- variants$INFO %>%
  purrr::map(function(x) {
    x %>%
      str_split("\\|") %>%
      unlist() %>%
      nth(5) 
  }) %>%
  unlist()

variants$variant <- variants$INFO %>%
  purrr::map(function(x) {
    x %>%
      str_split("\\|") %>%
      unlist() %>%
      nth(2) 
  }) %>%
  unlist()

variants <- variants %>%
  filter(variant == "missense_variant") %>%
  inner_join(gene_info, by = "annotation_id") %>%
  select(GENOME, gene_name) %>%
  drop_na()

list_shared_mutations <- function(pair, metadata, variants) {
  
  pids <- pair %>%
    str_split("___") %>%
    unlist()
  
  mut_a <- list_acquired_mutations(pids[1], metadata, variants)
  mut_b <- list_acquired_mutations(pids[2], metadata, variants)
  intersect(mut_a, mut_b)
  
}

pairs <- unique(dat$pair) 

shared_mutations <- pairs %>%
  purrr::map(function(x) list_shared_mutations(x, metadata, variants))
shared_mutations
