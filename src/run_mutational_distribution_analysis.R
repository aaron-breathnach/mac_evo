#!/usr/bin/env Rscript

library(tidyverse)
source("src/run_string_analysis.R")

"Usage:
   run_mutational_distribution_analysis.R [--panaroo <panaroo> --metadata <metadata> --vcf <vcf> --ref_cds <ref_cds> --out_dir <out_dir>]

Options:
   --panaroo Path to the directory containing output from Panaroo [default: data/panaroo]
   --metadata Isolate metadata [default: data/metadata.tsv]
   --vcf File containing recombination-filtered variant calls [default: data/filtered_variants.tsv]
   --ref_cds dndscv-formatted CDS table [default: data/ref_cds.tsv]
   --out_dir Output directory [default: data]
   
" -> doc

opts <- docopt::docopt(doc)

panaroo_to_vcf <- function(pid, dat, metadata, pres_abs) {
  
  meta <- filter(metadata, patient == pid)
  
  ref_genome <- meta %>%
    filter(time_from_diagnosis == 0) %>%
    pull(isolate)
  
  genes <- get_gene_to_anno(pres_abs)$gene_id
  
  ## loci that are shared with GCF_009741445.1 + the patient's reference genome
  gen_pre_abs <- pres_abs %>%
    dplyr::rename(gene_id = 1) %>%
    select(gene_id, all_of(ref_genome)) %>%
    dplyr::rename(annotation_id = 2) %>%
    filter(gene_id %in% genes) %>%
    separate_longer_delim(annotation_id, ";")
  
  inner_join(gen_pre_abs, dat, by = "annotation_id")
  
}

get_gene_to_anno <- function(pres_abs) {
  
  pres_abs %>%
    select(Gene, ncol(.)) %>%
    drop_na() %>%
    dplyr::rename("gene_id" = 1, "annotation_id" = 2)
  
}

get_uniq_vars <- function(dat, metadata) {
  
  inner_join(dat, metadata[,c("patient", "isolate")], by = c("GENOME" = "isolate")) %>%
    select(patient, gene_id, variant, POS, REF, ALT) %>%
    distinct()
  
}

calc_dn_ds <- function(dat, metadata) {
  
  unique_variants <- inner_join(dat, metadata[,c("patient", "isolate")], by = c("GENOME" = "isolate")) %>%
    select(patient, gene_id, variant, POS, REF, ALT) %>%
    distinct()
  
  unique_variants %>%
    filter(variant %in% c("missense_variant", "synonymous_variant")) %>%
    select(gene_id, variant) %>%
    group_by(gene_id, variant) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    pivot_wider(names_from = variant, values_from = n, values_fill = 0) %>%
    setNames(c("gene_id", "dn", "ds")) %>%
    mutate(dn_ds = dn / ds)
  
}

run_poisson_test <- function(dat, metadata, num_mut, total_length, gene_data) {
  
  uniq_vars <- inner_join(dat, metadata[,c("patient", "isolate")], by = c("GENOME" = "isolate")) %>%
    select(patient, gene_id, variant, POS, REF, ALT) %>%
    distinct()
  
  uniq_vars %>%
    filter(variant == "missense_variant") %>%
    group_by(gene_id) %>%
    tally() %>%
    ungroup() %>%
    inner_join(gene_data, by = "gene_id") %>%
    group_by(gene_id) %>%
    mutate(p = poisson.test(
      n,
      num_mut * (gene_length / total_length),
      alternative = "greater")$p.value
    ) %>%
    ungroup()
  
}

get_gen_pos <- function(dndscv, sel, sizes, dn_ds) {
  
  left_join(sel, sizes, by = "gene_id") %>%
    inner_join(dn_ds, by = "gene_id") %>%
    inner_join(dndscv, by = "annotation_id") %>%
    mutate(p = replace_na(p, 1)) %>%
    mutate(pos = pos / 1e6) %>%
    filter(!grepl("~~~", gene_id_prokka)) %>%
    mutate(
      bonferroni = p.adjust(p, method = "bonferroni"),
      fdr = p.adjust(p, method = "fdr"),
      holm = p.adjust(p, method = "holm")
    )
  
}

run_mutational_distribution_analysis <- function(PANAROO,
                                                 METADATA,
                                                 VARIANTS,
                                                 REF,
                                                 OUT_DIR = "data") {
  
  dir.create(OUT_DIR, FALSE, TRUE)
  
  ## panaroo presence/absence
  gene_presence_absence <- sprintf("%s/gene_presence_absence.csv", PANAROO) %>%
    read_delim()
  
  gene_to_anno <- get_gene_to_anno(gene_presence_absence)
  
  ## panaroo gene data
  gene_data <- sprintf("%s/gene_data.csv", PANAROO) %>%
    read_delim(col_select = c(4, 6)) %>%
    inner_join(gene_to_anno, by = "annotation_id") %>%
    mutate(gene_length = nchar(dna_sequence)) %>%
    select(1, 3, 4)
  
  ## genome metadata
  metadata <- read_delim(METADATA)
  
  ## variants
  dat <- read_delim(VARIANTS)
  
  dat$annotation_id <- dat$INFO %>%
    purrr::map(function(x) str_split(x, "\\|") %>% unlist() %>% nth(5)) %>%
    unlist()
  
  dat$variant <- dat$INFO %>%
    purrr::map(function(x) str_split(x, "\\|") %>% unlist() %>% nth(2)) %>%
    unlist()
  
  patients <- metadata %>%
    filter(bracken_pass & !multiple_carriage) %>%
    pull(patient) %>%
    unique()
  
  dat <- patients %>%
    purrr::map(function(x) panaroo_to_vcf(x, dat, metadata, gene_presence_absence)) %>%
    bind_rows() %>%
    filter(variant %in% c("missense_variant", "synonymous_variant")) %>%
    filter(!grepl("group_", gene_id))
  
  no_mutations <- dat %>%
    filter(variant == "missense_variant") %>%
    nrow()
  
  total_length <- sum(gene_data$gene_length)
  sel <- run_poisson_test(dat, metadata, no_mutations, total_length, gene_data)
  
  ## count the number of patients in whom variants were detected
  sizes <- dat %>%
    inner_join(metadata, by = c("GENOME" = "isolate")) %>%
    select(patient, gene_id) %>%
    distinct() %>%
    group_by(gene_id) %>%
    summarise(size = n())
  
  ## map Prokka gene names to annotation identifiers
  prok_to_anno <- read_delim(REF, col_select = c(1, 5)) %>%
    dplyr::rename(locus_tag = 1, pos = 2) %>%
    mutate(gene_id_prokka = gsub(".*\\:", "", locus_tag)) %>%
    filter(!grepl("IMFOIEKC_", gene_id_prokka)) %>%
    mutate(gene_id_prokka = gsub("_.*", "", gene_id_prokka)) %>%
    mutate(annotation_id = gsub("\\:.*", "", locus_tag))
  
  dn_ds <- calc_dn_ds(dat, metadata)
  
  gen_pos <- get_gen_pos(prok_to_anno, sel, sizes, dn_ds)
  
  write_tsv(gen_pos, sprintf("%s/gen_pos.tsv", OUT_DIR))
  
}

PANAROO <- opts$panaroo
METADATA <- opts$metadata
PATIENTS <- opts$patients
VARIANTS <- opts$vcf
REF <- opts$ref_cds
ARG <- opts$myco_args
OUT_DIR <- opts$out_dir

if (sys.nframe() == 0) {
  run_mutational_distribution_analysis(
    opts$panaroo,
    opts$metadata,
    opts$vcf,
    opts$ref_cds,
    opts$out_dir
  )
}
