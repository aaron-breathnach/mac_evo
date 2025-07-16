library(tidyverse)

sort_pids <- function(x, y) {
  values <- sort(c(x, y))
  paste0(values[1], "___", values[2])
}

list_acquired_mutations <- function(pid, metadata, variants) {
  
  final_isolate <- metadata %>%
    filter(patient == pid) %>%
    filter(time_from_diagnosis == max(time_from_diagnosis)) %>%
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

get_pairs <- function(cid, dat) {
  
  x <- dat %>%
    filter(clade == cid) %>%
    pull(patient) %>%
    unique()
  
  apply(combn(x, 2), 2, paste, collapse = "___")
  
}

pids <- readLines("data/patients.txt")

dat <- read_delim("data/snp_dists.tsv") %>%
  filter(patient.x != patient.y & dist <= 12)

pairs <- purrr::map(clades, function(x) get_pairs(x, dat)) %>%
  unlist()

gene_data <- read_delim("data/panaroo/gene_data.csv")

variants <- read_delim("data/filtered_variants.tsv") %>%
  get_annotation_ids() %>%
  filter(grepl("missense_variant", INFO)) %>%
  inner_join(gene_data, by = "annotation_id") %>%
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
