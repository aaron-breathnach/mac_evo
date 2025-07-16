library("tidyverse")

parse_dists <- function(f, metadata) {
  
  fastbaps_cluster <- f %>%
    str_split("/") %>%
    unlist() %>%
    nth(2)
  
  meta <- metadata %>%
    select(isolate, patient, country) %>%
    rename(name = 1)
  
  dat <- read_delim(f)
  
  if (nrow(dat) > 0) {
    
    dat %>%
      dplyr::rename("genome_1" = 1) %>%
      filter(genome_1 != "Reference") %>%
      select(-any_of("Reference")) %>%
      pivot_longer(!genome_1, names_to = "genome_2", values_to = "dist") %>%
      inner_join(meta, by = c("genome_1" = "name")) %>%
      inner_join(meta, by = c("genome_2" = "name")) %>%
      mutate(fastbaps_cluster = fastbaps_cluster) %>%
      select(ncol(.), 1:(ncol(.) - 1))
    
  }
}

metadata <- read_delim("data/metadata.tsv")

files <- list.files("ska_x_fastbaps",
                    pattern = "snp_dists.txt",
                    full.names = TRUE,
                    recursive = TRUE)

snp_dists <- purrr::map(files, function(x) parse_dists(x, metadata)) %>%
  bind_rows()

write_tsv(snp_dists, "data/snp_dists.tsv")
