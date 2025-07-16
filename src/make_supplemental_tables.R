library(tidyverse)

get_alias <- function(x) {
  x %>%
    str_replace("P", "P0") %>%
    str_replace("Wetzstein_", "W") %>%
    str_replace("VT_", "V")
}

make_table_s4 <- function() {
  
  mat <- readRDS("data/haplotype_deconstructor.RDS")$decomposed$signatures
  
  gff <- "prokka/reference/reference.gff"
  
  snps <- tibble(snp = rownames(mat)) %>%
    separate_wider_delim(
      snp,
      delim = ".",
      names = c("chrom", "pos", "mut"),
      cols_remove = FALSE
    ) %>%
    mutate(pos = as.numeric(pos))
  
  ref <- rtracklayer::readGFF(gff) %>%
    as.data.frame() %>%
    filter(type == "gene") %>%
    mutate(chrom = gsub(".*\\|", "", seqid)) %>%
    select(chrom, start, end, gene) %>%
    as_tibble()
  
  res <- fuzzyjoin::fuzzy_inner_join(
    snps,
    ref,
    by = c("chrom" = "chrom", "pos" = "start", "pos" = "end"),
    match_fun = list(`==`, `>=`, `<=`)
  ) %>%
    select(-chrom.y) %>%
    rename(chrom = chrom.x)
  
  dup_snp <- res %>%
    group_by(snp) %>%
    tally() %>%
    filter(n > 1) %>%
    pull(snp)
  
  part_1 <- res %>%
    filter(!snp %in% dup_snp)
  
  part_2 <- res %>%
    filter(snp %in% dup_snp) %>%
    drop_na() %>%
    group_by(snp) %>%
    sample_n(1) %>%
    ungroup()
  
  ann <- rbind(part_1, part_2) %>%
    separate_wider_delim(mut, "_", names = c("ref", "alt")) %>%
    select(snp, chrom, pos, ref, alt, gene)
  
  dat <- mat %>%
    as.data.frame() %>%
    rownames_to_column("snp") %>%
    as_tibble()
  
  table_s4 <- inner_join(ann, dat, by = "snp")
  
}

if (!dir.exists("supplemental_tables")) dir.create("supplemental_tables")

metadata <- read_delim("data/metadata.tsv") %>%
  mutate(patient_alias = patient %>%
           str_replace("P", "P0") %>%
           str_replace("Wetzstein_", "W") %>%
           str_replace("VT_", "V")) %>%
  select(1, ncol(.), 2:(ncol(.) - 1))

## table s1
table_s1 <- metadata %>%
  filter(study == "Present")

## table s2
table_s2 <- metadata %>%
  filter(study != "Present")

## table s3
table_s3 <- read_delim("data/patient_info.tsv")

## table s4
table_s4 <- make_table_s4()

## write the tables
tables   <- list(table_s1, table_s2, table_s3, table_s4)
files    <- sprintf("supplemental_tables/table_s%s.tsv", 1:4)
captions <- readLines("supplemental_tables/supplemental_table_captions.txt")

for (n in 1:4) {
  write_tsv(tables[[n]], files[n])
  lines <- c(captions[n], readLines(files[n]))
  writeLines(lines, files[n])
}
