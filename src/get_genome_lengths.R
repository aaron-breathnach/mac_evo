library(tidyverse)

get_genome_length <- function(isolate) {
  
  bases <- sprintf("prokka/%s/%s.txt", isolate, isolate) %>%
    readLines() %>%
    nth(3) %>%
    str_replace(".* ", "") %>%
    as.numeric()
  
  tibble(isolate = isolate, bases = bases)
  
}

get_genome_lengths <- function() {
  
  isolates <- list.files("prokka")
  
  genome_lengths <- purrr::map(isolates, function(x) get_genome_length(x)) %>%
    bind_rows()
  
  write_tsv(genome_lengths, "data/genome_lengths.tsv")
  
}

if (sys.nframe() == 0) {
  get_genome_lengths()
}
