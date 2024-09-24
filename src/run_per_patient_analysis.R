#!/usr/bin/env Rscript

library(tidyverse)
library(ggnewscale)
library(ggtree)
library(phytools)

source("src/compare_first_vs_last.R")

plot_tree <- function(pid, tree, metadata, show_legend = FALSE) {
  
  ann <- metadata %>%
    mutate(NODE = ifelse(patient == pid, 1, NA)) %>%
    select(isolate, patient, NODE)
  
  ggtree::ggtree(tree) %<+%
    ann +
    geom_tippoint(aes(subset = !is.na(NODE)),
                  size = 4,
                  colour = "red") +
    ggtitle(pid) +
    theme(legend.position = "none")
  
}

make_umap <- function(p_inp,
                      colour,
                      title = NULL,
                      show_legend = TRUE,
                      legend_title = NULL,
                      pal = NULL) {
  
  p <- ggplot(p_inp, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes_string(colour = colour),
               size = 4,
               show.legend = show_legend) +
    theme_classic() +
    theme(axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          legend.text = ggtext::element_markdown()) +
    labs(title = title, colour = legend_title)
  
  if (!is.null(pal)) {
    p <- p +
      scale_colour_manual(values = pal)
  }
  
  return(p)
  
}

umap_by_patient <- function(pid, p_inp) {
  
  df <- p_inp %>%
    mutate(colour = ifelse(patient == pid, "a", "b")) %>%
    arrange(desc(colour))
  
  make_umap(
    df,
    colour = "colour",
    show_legend = FALSE,
    pal = c("red", "lightgrey")
  )
  
}

run_umap <- function(snp_dists, metadata) {
  
  snp_dists <- snp_dists[rownames(snp_dists) %in% metadata$isolate,
                         colnames(snp_dists) %in% metadata$isolate]
  
  set.seed(321)
  snps_umap <- umap::umap(snp_dists)
  
  points <- snps_umap$layout %>%
    as.data.frame() %>%
    rownames_to_column("isolate") %>%
    rename(UMAP1 = 2, UMAP2 = 3) %>%
    filter(isolate != "Reference") %>%
    inner_join(metadata, by = "isolate")
  
  return(points)
  
}

plot_haplotype <- function(pid, haplo_relab, metadata, pal) {
  
  p_inp <- haplo_relab %>%
    inner_join(metadata, by = "isolate") %>%
    filter(patient == pid & rel_contribution > 0) %>%
    mutate(time_from_diagnosis = as.character(round(time_from_diagnosis, 2)))
  
  ggplot(p_inp, aes(x = time_from_diagnosis, y = contribution)) +
    geom_bar(
      aes(fill = haplotype),
      stat = "identity",
      position = "fill",
      colour = "black"
    ) +
    theme_bw(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold")) +
    scale_fill_manual(values = pal) +
    labs(x = "Time from diagnosis (years)", y = "Abundance", fill = "Haplotype")
    
}

wrapper <- function(pid, tree, umap, hapl, comp, metadata, pal) {
  
  p_tree <- plot_tree(pid, tree, metadata)
  
  p_umap <- umap_by_patient(pid, umap)
  
  p_hapl <- plot_haplotype(pid, haplo_relab, metadata, pal)
  
  p_comp <- first_vs_last_bar(pid, comp)
  
  patchwork::wrap_plots(p_tree, p_umap, p_hapl, p_comp, widths = c(1, 1, 2/3, 2/3))
  
}

per_patient_analysis <- function(TREE = "data/mavium.final_tree.tre",
                                 SNP_DISTS = "data/snp_dists.tsv",
                                 HAPLOTYPE = "data/haplotype_deconstructor.80_percent.RDS",
                                 METADATA = "data/metadata.tsv") {
  
  metadata <- read_delim(METADATA) %>%
    filter(study == "Present") %>%
    select(2, 1, 3:ncol(.))
  
  isolates <- metadata$isolate
  
  tree <- ggtree::read.tree(TREE) %>%
    ape::keep.tip(isolates)
  
  snp_dists <- read_delim(SNP_DISTS) %>%
    rename(isolate = 1) %>%
    column_to_rownames("isolate")
  
  umap <- run_umap(snp_dists, metadata)
  
  haplo_relab <- readRDS(HAPLOTYPE)[[2]]
  
  haplotypes <- sort(unique(haplo_relab$haplotype))
  pal <- RColorBrewer::brewer.pal(length(haplotypes), "Set3")
  names(pal) <- haplotypes
  
  comp <- compare_first_vs_last(haplo_relab, metadata)
  
  patients <- metadata %>%
    filter(study == "Present") %>%
    pull(patient) %>%
    unique()
  
  plot_list <- purrr::map(patients, function(x) wrapper(x, tree, umap, hapl, comp, metadata, pal))
  
  pdf("plots/per_patient_plots.pdf", width = 17.5, height = 5)
  print(plot_list)
  dev.off()
  
}

per_patient_analysis()
