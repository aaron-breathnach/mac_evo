#!/usr/bin/env Rscript

library(tidyverse)

"Usage:
   list_mycobacterial_args.R [--card <card> --out_dir <out_dir>]

Options:
   --card Path to the CARD data directory
   --out_dir Output directory [default: data]
   
" -> doc

opts <- docopt::docopt(doc)

list_card_mycobacterial_args <- function(card, out_dir) {
  
  dir.create(out_dir, FALSE, TRUE)
  
  filenames <- paste0("card-data/", c("shortname_pathogens.tsv", "shortname_antibiotics.tsv", "aro_index.tsv"))
  
  myc <- read_delim(filenames[1]) %>%
    filter(grepl("Mycobacterium", Pathogen)) %>%
    pull(Abbreviation)
  
  abx <- read_delim(filenames[2]) %>%
    dplyr::rename(abbreviation = 1, antibiotic = 2)
  
  aro <- read_delim(filenames[3]) %>%
    setNames(str_to_lower(gsub(" ", "_", names(.))))
  
  aro <- purrr::map(myc, function(x) filter(aro, grepl(x, card_short_name))) %>%
    bind_rows() %>%
    select(card_short_name) %>%
    separate_wider_delim(card_short_name,
                         "_",
                         names = c("species", "gene", "abbreviation")) %>%
    inner_join(abx, by = "abbreviation", relationship = "many-to-many") %>%
    select(-abbreviation) %>%
    group_by(gene) %>%
    summarise(antibiotic = antibiotic %>%
                paste0(collapse = ";") %>%
                str_split(";") %>%
                unlist() %>%
                unique() %>%
                paste0(collapse = ";"))
  
  write_tsv(aro, sprintf("%s/card_mycobacterial_args.tsv", out_dir))
  
}

if (sys.nframe() == 0) {
  list_card_mycobacterial_args(opts$card, opts$out_dir)
}
