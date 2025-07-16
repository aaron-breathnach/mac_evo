library(tidyverse)

metadata <- read_delim("data/metadata.tsv")

read_delim("data/fastbaps.tsv") %>%
  inner_join(metadata, by = "isolate") %>%
  select(level_1, patient) %>%
  distinct() %>%
  group_by(level_1) %>%
  tally() %>%
  ungroup() %>%
  filter(n > 1) %>%
  pull(level_1) %>%
  as.character() %>%
  writeLines("data/fastbaps_clusters.txt")

