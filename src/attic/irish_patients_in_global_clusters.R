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
  filter(dist <= 12)

strains <- snp_dists %>%
  dplyr::rename(genome_1 = 1) %>%
  pivot_longer(!genome_1, names_to = "genome_2", values_to = "dist") %>%
  filter(genome_1 != genome_2) %>%
  filter(dist <= 3)

clusters <- strains %>%
  tidygraph::as_tbl_graph(directed = FALSE) %>%
  igraph::cluster_leiden()

clusters <- tibble(isolate = clusters$names, cluster = clusters$membership) %>%
  inner_join(metadata, by = "isolate")

present <- clusters %>%
  group_by(cluster, study) %>%
  tally() %>%
  ungroup() %>%
  pivot_wider(names_from = "study", values_from = n, values_fill = 0) %>%
  filter(Present > 0 & (`Van Tonder` + Wetzstein) > 0) %>%
  pull(cluster)

patient_info <- read_delim("data/patient_info.tsv")

patients <- clusters %>%
  filter(cluster %in% present & study == "Present") %>%
  select(patient, cluster) %>%
  mutate(cluster = sprintf("cluster_%s", str_pad(cluster, 2, "left", "0"))) %>%
  distinct() %>%
  arrange(patient) %>%
  arrange(cluster) %>%
  inner_join(patient_info, by = "patient")

write_tsv(patients, "data/patients_trans_clust.tsv")
