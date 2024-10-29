#!/usr/bin/env Rscript

library(tidyverse)

"Usage:
   run_string_analysis.R [--sig_gen <sig_gen> --score <score> --limit <limit> --prefix <prefix> --out_dir <out_dir>]

Options:
   --sig_gen Table containing a list of genes in with enriched mutations [default: data/sig_gen.tsv]
   --score Minimum required interaction score [default: 400]
   --limit Max number of interactors to show [default: 30]
   --prefix Prefix for output filenames [default: network]
   --out_dir Output directory [default: data]
   
" -> doc

opts <- docopt::docopt(doc)

get_component_membership <- function(x) {
  
  igraph::components(x)$membership %>%
    data.frame() %>%
    rownames_to_column("name") %>%
    dplyr::rename(membership = 2)
  
}

query_string_db <- function(dat, score, limit, prefix, out_dir) {
  
  dir.create(out_dir, FALSE, TRUE)
  
  df <- dat %>%
    setNames(c("name", "node_size"))
  
  gois <- unique(df$name)
  
  ids <- paste0(gois, collapse = "%0d")
  
  part_1 <- "https://string-db.org/api/tsv/network"
  part_2 <- paste0("?identifiers=", ids)
  part_3 <- paste0("&limit=", limit)
  part_4 <- "&species=83332"
  part_5 <- paste0("&required_score=", score)
  string <- paste0(part_1, part_2, part_3, part_4, part_5)
  
  out <- sprintf("%s/%s.limit_%s.score_%s.tsv", out_dir, prefix, limit, score)
  
  utils::download.file(string, out)
  
}

get_enrichment <- function(gois, out_dir) {
  
  ids <- paste0(gois, collapse = "%0d")
  
  enrichment <- sprintf(
    "https://string-db.org/api/tsv/enrichment?identifiers=%s",
    ids
  )
  
  filename <- sprintf("%s/enrichment.tsv", out_dir)
  
  download.file(enrichment, filename)
  
}

get_functional_annotations <- function(gene, out_dir) {
  
  functional_annotation <- sprintf(
    "https://string-db.org/api/tsv/functional_annotation?identifiers=%s&species=83332",
    gene
  )
  
  dir.create(out_dir, FALSE, TRUE)
  
  filename <- sprintf("%s/%s.tsv", out_dir, gene)
  
  try(download.file(functional_annotation, filename))
  
}

run_string_analysis <- function(gen_pos, score, limit, prefix, out_dir) {
  
  dir.create(out_dir, FALSE, TRUE)
  
  dat <- read_delim(gen_pos) %>%
    filter(bonferroni < 0.05) %>%
    select(gene_id_prokka, size) %>%
    setNames(c("gene", "size"))
  
  query_string_db(dat, score, limit, prefix, out_dir)
  
  filename <- sprintf("%s/%s.limit_%s.score_%s.tsv", out_dir, prefix, limit, score)
  
  genes <- read_delim(filename) %>%
    select(3, 4) %>%
    pivot_longer(cols = everything()) %>%
    filter(value %in% dat$gene) %>%
    pull(value) %>%
    unique()
  
  get_enrichment(genes, out_dir)
  
}

if (sys.nframe() == 0) {
  run_string_analysis(opts$sig_gen, opts$score, opts$limit, opts$prefix, opts$out_dir)
}
