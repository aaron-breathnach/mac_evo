#!/usr/bin/env Rscript

library(tidyverse)

"Usage:
   list_reinfected_patients.R [--fastbaps <fastbaps> --metadata <metadata>]

Options:
   --fastbaps fastbaps output [default: data/fastbaps.tsv]
   --metadata Genome metadata [default: data/metadata.tsv]

" -> doc 

opts <- docopt::docopt(doc)

list_reinfected_patients <- function(METADATA, FASTBAPS) {
  
  metadata <- read_delim(METADATA)
  
  reinfected_patients <- read_delim("data/fastbaps.tsv") %>%
    inner_join(metadata, by = "isolate") %>%
    select(patient, level_1) %>%
    distinct() %>%
    group_by(patient) %>%
    tally() %>%
    filter(n > 1) %>%
    pull(patient)
  
  writeLines(reinfected_patients, "data/reinfected_patients.txt")
  
  patients <- metadata %>%
    filter(!patient %in% reinfected_patients) %>%
    pull(patient) %>%
    unique() %>%
    sort()
  
  writeLines(patients, "data/patients.txt")
  
}

if (sys.nframe() == 0) {
  list_reinfected_patients(opts$metadata, opts$fastbaps)
}
