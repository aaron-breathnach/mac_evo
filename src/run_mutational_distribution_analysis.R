#!/usr/bin/env Rscript

library(tidyverse)
source("src/run_string_analysis.R")

"Usage:
   run_mutational_distribution_analysis.R [--panaroo <panaroo> --metadata <metadata> --vcf <vcf> --ref_cds <ref_cds> --myco_args <myco_args> --out_dir <out_dir>]

Options:
   --panaroo Path to the directory containing output from Panaroo [default: data/panaroo]
   --metadata Isolate metadata [default: data/metadata.tsv]
   --vcf File containing recombination-filtered variant calls [default: data/filtered_variants.tsv]
   --ref_cds dndscv-formatted CDS table [default: data/ref_cds.tsv]
   --myco_args Mycobacterial antibiotic resistance genes [default: data/card_mycobacterial_args.tsv]
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

run_poisson_test <- function(dat, no_mutations, total_length, gene_data) {
  
  dat %>%
    filter(variant %in% c("missense_variant", "synonymous_variant")) %>%
    select(gene_id, variant) %>%
    group_by(gene_id, variant) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    pivot_wider(names_from = variant, values_from = n, values_fill = 0) %>%
    setNames(c("gene_id", "dn", "ds")) %>%
    inner_join(gene_data, by = "gene_id") %>%
    group_by(gene_id) %>%
    mutate(p = poisson.test(dn,
                            r = no_mutations * (gene_length / total_length),
                            alternative = "greater")$p.value) %>%
    ungroup() %>%
    mutate(fdr = p.adjust(p, method = "fdr"))
  
}

get_gen_pos <- function(dndscv, sel, sizes, amr_gen) {
  
  left_join(sel, sizes, by = "gene_id") %>%
    inner_join(dndscv, by = "annotation_id") %>%
    mutate(p = replace_na(p, 1)) %>%
    mutate(fdr = replace_na(fdr, 1)) %>%
    mutate(pos = pos / 1e6) %>%
    filter(!grepl("~~~", gene_id_prokka)) %>%
    mutate(SIZE = ifelse(fdr <= 0.05 | gene_id_prokka %in% amr_gen, size, 0)) %>%
    mutate(sig = ifelse(fdr > 0.05, "n", "y"))
  
}

export_results <- function(gen_pos, out_dir) {
  
  sig_gen <- gen_pos %>%
    filter(fdr <= 0.05) %>%
    select(gene_id_prokka, size)
  
  filename <- sprintf("%s/sig_gen.tsv", out_dir)
  
  write_tsv(sig_gen, filename)
  
}

run_mutational_distribution_analysis <- function(PANAROO,
                                                 METADATA,
                                                 VARIANTS,
                                                 REF,
                                                 ARG,
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
  
  ## list of patients to include in the analysis
  patients <- metadata %>%
    filter(bracken_pass & !multiple_carriage) %>%
    pull(patient) %>%
    unique()
  
  ## variants
  dat <- read_delim(VARIANTS)
  
  dat$annotation_id <- dat$INFO %>%
    purrr::map(function(x) str_split(x, "\\|") %>% unlist() %>% nth(5)) %>%
    unlist()
  
  dat$variant <- dat$INFO %>%
    purrr::map(function(x) str_split(x, "\\|") %>% unlist() %>% nth(2)) %>%
    unlist()
  
  dat <- patients %>%
    purrr::map(function(x) panaroo_to_vcf(x, dat, metadata, gene_presence_absence)) %>%
    bind_rows() %>%
    filter(variant %in% c("missense_variant", "synonymous_variant")) %>%
    filter(!grepl("group_", gene_id))
  
  ## adapted from Tonkin-Hill et al. 2022
  too_close <- unlist(
    imap(
      split(dat$POS, paste(dat$CHROM, dat$GENOME)), 
      function(pos, mut){
        if (length(pos) > 1){
          cmb <- combinat::combn(pos, 2, simplify = FALSE)
          cmb <- unique(unlist(cmb[map_dbl(cmb, diff) <= 150]))
          return(paste(mut, cmb))
        } else{
          return(NULL)
        }
      }
    )
  )
  
  dat <- dat %>%
    filter(!paste(CHROM, GENOME, POS) %in% too_close)
  
  no_mutations <- nrow(dat)
  total_length <- sum(gene_data$gene_length)
  sel <- run_poisson_test(dat, no_mutations, total_length, gene_data)
  
  ## count the number of patients in whom variants were detected
  sizes <- dat %>%
    inner_join(metadata, by = c("GENOME" = "isolate")) %>%
    select(patient, gene_id) %>%
    distinct() %>%
    group_by(gene_id) %>%
    summarise(size = n())
  
  ## list of mycobacterial antibiotic resistance genes
  amr_gen <- read_delim(ARG) %>%
    pull(gene)
  
  ## map Prokka gene names to annotation identifiers
  prok_to_anno <- read_delim(REF, col_select = c(1, 5)) %>%
    dplyr::rename(locus_tag = 1, pos = 2) %>%
    mutate(gene_id_prokka = gsub(".*\\:", "", locus_tag)) %>%
    filter(!grepl("IMFOIEKC_", gene_id_prokka)) %>%
    mutate(gene_id_prokka = gsub("_.*", "", gene_id_prokka)) %>%
    mutate(annotation_id = gsub("\\:.*", "", locus_tag))
  
  gen_pos <- get_gen_pos(prok_to_anno, sel, sizes, amr_gen)
  
  sig_gen <- gen_pos %>%
    filter(fdr <= 0.05) %>%
    select(gene_id, gene_id_prokka, size)
  
  missense_variants_per_patient <- dat %>%
    filter(variant == "missense_variant") %>%
    inner_join(sig_gen, by = "gene_id") %>%
    inner_join(metadata, by = c("GENOME" = "isolate")) %>%
    select(patient, gene_id, gene_id_prokka, GENOME, variant)
  
  write_tsv(missense_variants_per_patient, "data/missense_variants_per_patient.tsv")
  
  ## export the results
  write_tsv(gen_pos, sprintf("%s/gen_pos.tsv", OUT_DIR))
  write_tsv(sig_gen, sprintf("%s/sig_gen.tsv", OUT_DIR))
  
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
    opts$myco_args,
    opts$out_dir
  )
}
