#!/usr/bin/env Rscript

library(tidyverse)

"Usage:
   make_dndscv_fmt_ref_tab.R [--gff <gff> --out_dir <out_dir>]

Options:
   --gff Reference genome GFF file
   --out_dir Output directory
   
" -> doc

opts <- docopt::docopt(doc)

make_dndscv_fmt_ref_tab <- function(gff, out) {
  
  cols <- c("gene.id",
            "gene.name",
            "cds.id",
            "chr",
            "chr.coding.start",
            "chr.coding.end",
            "cds.start",
            "cds.end",
            "length",
            "strand")
  
  gff <- rtracklayer::readGFF(gff) %>%
    as.data.frame() %>%
    as_tibble() %>%
    select(seqid, type, start, end, strand, gene, locus_tag) %>%
    mutate(strand = ifelse(strand == "+", 1, -1)) %>%
    mutate(gene = ifelse(is.na(gene), locus_tag, gene))
  
  part_1 <- gff %>%
    filter(type == "gene") %>%
    select(gene, locus_tag) %>%
    mutate(gene.id = ifelse(is.na(gene),
                            sprintf("%s:%s", locus_tag, locus_tag),
                            sprintf("%s:%s", locus_tag, gene))) %>%
    mutate(gene.name = gene.id) %>%
    select(locus_tag, gene.id, gene.name)
  
  part_2 <- gff %>%
    filter(type == "gene") %>%
    select(locus_tag, seqid, start, end) %>%
    setNames(c("locus_tag", "chr", "chr.coding.start", "chr.coding.end"))
  
  cols_p3 <- c("locus_tag", "cds.start", "cds.end", "strand", "length")
  
  part_3 <- gff %>%
    filter(type == "CDS") %>%
    select(locus_tag, start, end, strand) %>%
    mutate(cds.start = 1) %>%
    mutate(cds.end = abs(end - (start - 1))) %>%
    mutate(length = cds.end) %>%
    select(all_of(cols_p3)) %>%
    setNames(c("locus_tag", "cds.start", "cds.end", "strand", "length"))
  
  ref_table <- inner_join(part_1, part_2, by = "locus_tag") %>%
    inner_join(part_3, by = "locus_tag") %>%
    mutate(cds.id = gene.id) %>%
    select(all_of(cols)) %>%
    arrange(cds.start)
  
  dir.create(out_dir, TRUE, FALSE)
  filename <- sprintf("%s/ref_cds.tsv")
  write_tsv(ref_table, filename)
  
}

if (sys.nframe() == 0) {
  make_dndscv_fmt_ref_tab(opts$gff, opts$out_dir)
}
