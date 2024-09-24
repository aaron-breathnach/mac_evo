library(tidyverse)

metadata <- read_delim("data/metadata.tsv") %>%
  filter(study == "Present") %>%
  arrange(time_from_diagnosis) %>%
  group_by(patient) %>%
  mutate(timepoint = row_number()) %>%
  ungroup() %>%
  select(patient, isolate, timepoint, multiple_carriage) %>%
  arrange(patient)

dat <- read_delim("data/snp_dists.tsv")

t0 <- metadata %>%
  group_by(patient) %>%
  filter(timepoint == min(timepoint)) %>%
  pull(isolate)

df <- dat %>%
  rename(genome_1 = 1) %>%
  filter(genome_1 %in% t0) %>%
  select(1, all_of(metadata$isolate)) %>%
  pivot_longer(!genome_1, names_to = "genome_2", values_to = "dist") %>%
  filter(genome_1 != genome_2) %>%
  inner_join(metadata, by = c("genome_1" = "isolate")) %>%
  inner_join(metadata, by = c("genome_2" = "isolate")) %>%
  filter(patient.x == patient.y) %>%
  select(patient.x, timepoint.y, dist, multiple_carriage.x) %>%
  setNames(c("patient", "timepoint", "dist", "multiple_carriage")) %>%
  mutate(timepoint = paste0("T", timepoint)) %>%
  pivot_wider(names_from = "timepoint", values_from = "dist") %>%
  arrange(patient)

write_tsv(df, "data/table_s3.tsv")

info <- read_delim("data/patient_info.tsv") %>%
  inner_join(df, by = "patient")

write_tsv(info, "meeting.tsv")
