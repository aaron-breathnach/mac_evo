library(tidyverse)

mutation_rate <- as.numeric(readLines("data/mutation_rate.txt"))

metadata <- read_delim("data/metadata.tsv") %>%
  filter(!multiple_carriage & study == "Present") %>%
  select(isolate, patient, time_from_diagnosis)

snp_dists <- read_delim("data/snp_dists.tsv") %>%
  filter(genome_1 != genome_2)

dists <- snp_dists %>%
  filter(patient.x == patient.y) %>%
  select(patient.x, dist) %>%
  distinct() %>%
  group_by(patient.x) %>%
  filter(dist == max(dist)) %>%
  ungroup() %>%
  mutate(
    IQR = IQR(dist),
    O_upper = quantile(dist, probs = 0.75, na.rm = FALSE) + 3 * IQR,
    O_lower = quantile(dist, probs = 0.25, na.rm = FALSE) - 3 * IQR
  ) %>%
  filter(O_lower <= dist & dist <= O_upper) %>%
  pull(dist) 

x <- dists %>%
  quantile(probs = 0.95) %>%
  as.numeric()

threshold <- floor(x + mutation_rate)

writeLines(as.character(threshold), "data/snp_threshold.txt")
