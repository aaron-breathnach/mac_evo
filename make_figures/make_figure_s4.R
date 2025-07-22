library(tidyverse)

map_snp2gen <- function(mat, gff = "prokka/reference/reference.gff") {
  
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
  
  rbind(part_1, part_2) %>%
    separate_wider_delim(mut, "_", names = c("ref", "alt")) %>%
    select(snp, chrom, pos, ref, alt, gene)
  
}

mat2dat <- function(x) {
  x %>%
    as.data.frame() %>%
    rownames_to_column("isolate") %>%
    as_tibble()
}

## haplotype frequency x snp frequency plots
.viz_hf_x_maf <- function(hap, metadata, hap_x_snp, mat_hap_sub, mat_snp) {
  
  snps <- hap_x_snp %>%
    filter(haplotype == hap) %>%
    pull(snp)
  
  tmp <- mat_hap_sub %>%
    filter(haplotype == hap)
  
  haps <- tmp %>%
    dplyr::rename(id = 2, frequency = 3) %>%
    mutate(feature = "HF")
  
  isolates <- mat_hap_sub$isolate
  
  p_inp <- mat2dat(mat_snp) %>%
    pivot_longer(!isolate, names_to = "id", values_to = "frequency") %>%
    filter(isolate %in% isolates & id %in% snps) %>%
    mutate(haplotype = hap) %>%
    mutate(feature = "HF * MAF") %>%
    select(isolate, id, frequency, feature) %>%
    inner_join(tmp[,c(1,3)], by = "isolate") %>%
    mutate(frequency = frequency * abundance) %>%
    select(1:4) %>%
    bind_rows(haps) %>%
    inner_join(metadata, by = "isolate")
  
  p_inp_a <- p_inp %>%
    filter(feature != "HF")
  
  p_inp_b <- p_inp %>%
    filter(feature == "HF")
  
  ggplot(p_inp_a, aes(x = time_from_diagnosis, y = frequency)) +
    geom_line(aes(group = id),
              show.legend = FALSE,
              linewidth = 0.1,
              colour = "steelblue",
              alpha = 0.25) +
    geom_line(data = p_inp_b, aes(group = id),
              show.legend = FALSE,
              linewidth = 1,
              colour = "firebrick",
              linetype = "dashed") +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold"),
          legend.title = element_blank(),
          plot.title = element_text(face = "bold")) +
    labs(x = "Time from diagnosis (years)",
         y = "Frequency",
         title = hap)
  
}

viz_hf_x_maf <- function(pid, metadata, mat_hap, mat_snp, hap_x_snp) {
  
  if (!dir.exists("output")) dir.create("output")
  
  meta <- metadata %>%
    filter(patient == pid)
  
  mat_hap_sub <- mat2dat(mat_hap) %>%
    filter(isolate %in% meta$isolate) %>%
    pivot_longer(!isolate, names_to = "haplotype", values_to = "abundance") %>%
    filter(abundance > 0)
  
  haplotypes <- mat_hap_sub %>%
    filter(abundance > 0) %>%
    group_by(haplotype) %>%
    tally() %>%
    ungroup() %>%
    filter(n >= 3) %>%
    pull(haplotype) %>%
    sort()
  
  mat_hap_sub <- mat_hap_sub %>%
    filter(haplotype %in% haplotypes)
  
  plot_list <- haplotypes %>%
    purrr::map(function(x) .viz_hf_x_maf(x, meta, hap_x_snp, mat_hap_sub, mat_snp))
  
  num_hap <- length(haplotypes)
  
  num_row <- ceiling(num_hap / 3)
  
  if (num_row == 1) {
    h <- 3.75
    w <- num_hap * 3.75
  } else {
    h <- num_row * 3.75
    w <- 11.25
  }
  
  num_col <- ifelse(num_hap > 3, 3, num_hap)
  
  cowplot::plot_grid(plotlist = plot_list, nrow = num_row, ncol = num_col)
  
}

make_figure_s4 <- function() {
  
  haplo <- readRDS("data/haplotype_deconstructor.RDS")
  
  hap_x_snp <- haplo$decomposed$signatures %>%
    as.data.frame() %>%
    rownames_to_column("snp") %>%
    pivot_longer(!snp, names_to = "haplotype", values_to = "loading") %>%
    filter(loading > 0)
  
  ## haplotype abundances
  mat_hap <- haplo$decomposed$samples
  ## minor allele frequencies
  mat_snp <- t(haplo$decomposed$observed)
  
  ## correlation analysis
  # corr <- Hmisc::rcorr(x = mat_hap, y = mat_snp)
  # 
  # rval <- corr$r[1:ncol(mat_hap), (ncol(mat_hap) + 1):ncol(corr$r)] %>%
  #   as.data.frame() %>%
  #   rownames_to_column("haplotype") %>%
  #   pivot_longer(!haplotype, names_to = "snp", values_to = "r")
  # 
  # order <- corr$r[1:ncol(mat_hap), (ncol(mat_hap) + 1):ncol(corr$r)] %>%
  #   t() %>%
  #   dist() %>%
  #   hclust() %>%
  #   as.dendrogram() %>%
  #   labels()
  # 
  # rval$snp <- factor(rval$snp, levels = order)
  # 
  # p1 <- ggplot(rval, aes(x = haplotype, y = snp)) +
  #   geom_tile(aes(fill = r)) +
  #   theme_minimal(base_size = 12.5) +
  #   theme(axis.text.y = element_blank(),
  #         axis.ticks.y = element_blank(),
  #         axis.title = element_text(face = "bold"),
  #         legend.title = element_text(face = "bold.italic")) +
  #   scale_fill_gradient2(
  #     low = "blue",
  #     mid = "white",
  #     high = "red"
  #   ) +
  #   scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  #   scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  #   labs(x = "Haplotype", y = "SNP", fill = "r")
  
  metadata <- read_delim("data/metadata.tsv")
  
  pids <- metadata %>%
    filter(multiple_carriage) %>%
    group_by(patient) %>%
    tally() %>%
    filter(n > 3) %>%
    pull(patient)
  
  p2 <- viz_hf_x_maf("Wetzstein_033", metadata, mat_hap, mat_snp, hap_x_snp)
  
  # p <- cowplot::plot_grid(p1, p2,
  #                         rel_widths = c(1, 2.5),
  #                         scale = 0.95,
  #                         labels = "AUTO")
  
  ggsave("plots/figure_s4.png", p2, height = 3.75, width = 10, bg = "white")
  
}
