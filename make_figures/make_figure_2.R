library(tidyverse)

`%<+%` <- ggtree::`%<+%`

list_pids <- function(metadata) {
  
  patients <- metadata %>%
    filter(multiple_carriage) %>%
    pull(patient) %>%
    unique()
  
  meta <- metadata %>%
    filter(study == "Present" & patient %in% patients)
  
  unique(meta$patient)
  
}

get_jumps <- function(pid, metadata, fastbaps) {
  
  meta <- metadata %>%
    filter(patient == pid) %>%
    select(isolate, time_from_diagnosis)
  
  n <- nrow(meta)
  
  pairs <- 2:n %>%
    purrr::map(function(x) tibble(genome_1 = meta[[x - 1, 1]], genome_2 = meta[[x, 1]])) %>%
    bind_rows()
  
  pairs %>%
    inner_join(fastbaps, by = c("genome_1" = "isolate")) %>%
    inner_join(fastbaps, by = c("genome_2" = "isolate")) %>%
    filter(level_1.x != level_1.y) %>%
    mutate(patient = pid) %>%
    select(ncol(.), 1:(ncol(.) - 1))
  
}

make_figure_2a <- function(fastbaps, metadata) {
  
  p_inp <- inner_join(fastbaps, metadata, by = "isolate") %>%
    select(patient, level_1) %>%
    distinct() %>%
    group_by(patient) %>%
    tally() %>%
    ungroup() %>%
    as.data.frame() %>%
    rename(lineages_per_patient = 2) %>%
    mutate(lineages_per_patient = lineages_per_patient %>%
             english::english() %>%
             str_to_sentence()) %>%
    group_by(lineages_per_patient) %>%
    tally()
  
  ggplot(p_inp, aes(x = lineages_per_patient, y = n)) +
    geom_bar(stat = "identity", colour = "black", fill = "steelblue") +
    geom_text(aes(label = paste0("n=", n)), vjust = -0.25) +
    theme_classic() +
    theme(axis.title = element_text(face = "bold", size = 12.5),
          axis.text = element_text(size = 12.5)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Num. lineages", y = "Num. patients")
  
}

make_figure_2b <- function(min_snp_dis, pal) {
  
  ggplot(min_snp_dis, aes(x = patient, y = dist)) +
    geom_bar(aes(fill = patient),
             stat = "identity",
             colour = "black",
             show.legend = FALSE) +
    geom_text(aes(label = paste0("n=", dist)), vjust = -0.25) +
    theme_classic() +
    theme(axis.title = element_text(face = "bold", size = 12.5),
          axis.text = element_text(size = 12.5)) +
    scale_fill_manual(values = pal) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Patient", y = "Num. SNPs")
  
}

make_figure_2c <- function(pal, tree, metadata, taxalinks) {
  
  pids <- list_pids(metadata)
  
  tippoints <- metadata %>%
    mutate(patient = ifelse(patient %in% pids, patient, NA)) %>%
    select(isolate, patient)
  
  ggtree::ggtree(tree) %<+%
    tippoints +
    ggtree::geom_taxalink(
      data = taxalinks,
      aes(taxa1 = genome_1, taxa2 = genome_2),
      colour = "black",
      size = 1,
    ) +
    ggtree::geom_taxalink(
      data = taxalinks,
      aes(taxa1 = genome_1, taxa2 = genome_2, colour = patient),
      show.legend = FALSE
    ) +
    ggtree::geom_tippoint(
      aes(fill = patient, subset = !is.na(patient)),
      pch = 21,
      size = 3
    ) +
    scale_size_continuous(range = c(0, 2)) +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = pal) +
    labs(fill = "Patient") +
    theme(legend.text = element_text(size = 12.5),
          legend.title = element_text(face = "bold", size = 12.5)) +
    ggnewscale::new_scale_colour() +
    ggnewscale::new_scale_fill()
  
}

make_figure_2d <- function(pal, metadata, snp_dists, threshold) {
  
  metadata <- metadata %>%
    select(isolate, patient) %>%
    rename(name = 1)
  
  snp_mat <- snp_dists %>%
    select(2, 3, 4) %>%
    filter(genome_1 %in% metadata$name & genome_2 %in% metadata$name)
  
  duplicates <- snp_mat %>%
    filter(dist == 0 & genome_1 != genome_2) %>%
    tidygraph::as_tbl_graph(directed = FALSE) %>%
    igraph::cluster_leiden()
  
  non_duplicates <- metadata %>%
    filter(!name %in% duplicates$names) %>%
    pull(name)
  
  genome_groups <- tibble(
    name = duplicates$names,
    cluster = duplicates$membership
  )
  
  g <- snp_mat %>%
    filter(dist <= threshold) %>%
    tidygraph::as_tbl_graph()
  
  components <- igraph::components(g)$membership %>%
    as.data.frame() %>%
    rename(component = 1) %>%
    rownames_to_column("name")
  
  keep <- components %>%
    inner_join(metadata, by = "name") %>%
    select(component, patient) %>%
    distinct() %>%
    group_by(component) %>%
    tally() %>%
    ungroup() %>%
    filter(n > 1) %>%
    pull(component)
  
  isolates_a <- components %>%
    filter(component %in% keep) %>%
    pull(name)
  
  cids <- genome_groups %>%
    inner_join(metadata, by = "name") %>%
    select(cluster, patient) %>%
    distinct() %>%
    group_by(cluster) %>%
    tally() %>%
    filter(n > 1) %>%
    pull(cluster)
  
  isolates_b <- genome_groups %>%
    filter(cluster %in% cids) %>%
    pull(name)
  
  isolates <- unique(c(isolates_a, isolates_b))
  
  iden_geno_reps <- genome_groups %>%
    group_by(cluster) %>%
    sample_n(1) %>%
    ungroup() %>%
    rename(representative = 1) %>%
    select(representative, cluster)
  
  dat_a <- tibble(representative = non_duplicates, name = non_duplicates)
  
  dat_b <- inner_join(iden_geno_reps, genome_groups, by = "cluster") %>%
    select(representative, name)
  
  dat <- bind_rows(dat_a, dat_b)
  
  df_1 <- snp_mat %>%
    filter(dist <= threshold & genome_1 != genome_2) %>%
    tidygraph::as_tbl_graph() %>%
    filter(name %in% isolates) %>%
    inner_join(dat, by = "name") %>%
    inner_join(metadata, by = "name")
  
  df_a <- as_tibble(df_1) %>%
    select(representative, patient) %>%
    distinct() %>%
    as_tibble() %>%
    group_by(representative) %>%
    summarise(size = n())
  
  df_b <- as_tibble(df_1) %>%
    select(representative, patient) %>%
    distinct() %>%
    mutate(n = 1) %>%
    arrange(patient) %>%
    pivot_wider(names_from = patient, values_from = n, values_fill = 0)
  
  df <- inner_join(df_b, df_a, by = "representative")
  
  g <- df_1 %>%
    filter(name %in% df$representative) %>%
    inner_join(df, by = "representative") %>%
    igraph::mst() %>%
    tidygraph::as_tbl_graph()
  
  set.seed(123)
  layout <- ggraph::create_layout(g, "fr")
  igraph::V(g)$x <- layout$x
  igraph::V(g)$y <- layout$y
  
  clades <- components %>%
    filter(name %in% pull(g, name)) %>%
    select(component) %>%
    distinct() %>%
    arrange(component)
  
  ann <- components %>%
    filter(name %in% pull(g, name)) %>%
    inner_join(layout[, c("name", "x", "y")], by = "name") %>%
    group_by(component) %>%
    summarise(x = mean(x), y = 1.05 * max(y)) %>%
    ungroup() %>%
    inner_join(clades, by = "component") %>%
    arrange(desc(y)) %>%
    mutate(clade = paste0("Cluster ", row_number()))
  
  ggraph::ggraph(g, "manual", x = igraph::V(g)$x, y = igraph::V(g)$y) +
    ggraph::geom_edge_fan(aes(label = dist), colour = "lightgrey") +
    scatterpie::geom_scatterpie(data = layout,
                                aes(x = x, y = y, r = 0.2),
                                cols = colnames(df_b)[-1],
                                colour = "black") +
    scatterpie::geom_scatterpie(data = layout,
                                aes(x = x, y = y, r = 0.2),
                                cols = colnames(df_b)[-1],
                                colour = NA) +
    geom_text(data = ann,
              aes(x = x, y = y, label = clade),
              vjust = -2.5) +
    coord_fixed() +
    scale_fill_brewer(palette = "Set3") +
    scale_fill_manual(values = pal) +
    theme(legend.text = element_text(size = 12.5),
          legend.title = element_text(face = "bold", size = 12.5),
          panel.background = element_blank(),
          title = element_text(face = "bold")) +
    labs(fill = "Patient") +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.125)))
  
}

make_figure_2 <- function() {
  
  set_1 <- c("P02", "P04", "P06", "P07")
  set_2 <- c("P02", "P04", "P06", "P10", "P16", "P19")
  set <- sort(unique(c(set_1, set_2)))
  pal <- RColorBrewer::brewer.pal(length(set), "Set3")
  names(pal) <- set
  
  ## genome metadata
  metadata <- read_delim("data/metadata.tsv") %>%
    filter(study == "Present")
  
  ## species-wide SNP distances
  snp_dists_species <- read_delim("data/snp_dists.species_wide.txt") %>%
    rename(genome_1 = 1) %>%
    pivot_longer(!genome_1, names_to = "genome_2", values_to = "dist")
  
  ## cluster-wide SNP distances
  snp_dists_cluster <- read_delim("data/snp_dists.cluster_wide.tsv")
  
  ## fastbaps clusters
  fastbaps <- read_delim("data/fastbaps.tsv", col_select = c(1, 2))
  
  ## tree from IQ-TREE
  tree <- ggtree::read.tree("data/iqtree.treefile") %>%
    ape::keep.tip(metadata$isolate)
  
  ## snp threshold
  threshold <- as.numeric(readLines("data/threshold.txt"))
  
  pids <- unique(metadata$patient)
  
  taxalinks <- purrr::map(pids, function(x) get_jumps(x, metadata, fastbaps)) %>%
    bind_rows()
  
  min_snp_dis <- inner_join(taxalinks, snp_dists_species, by = c("genome_1", "genome_2")) %>%
    group_by(patient) %>%
    filter(dist == min(dist)) %>%
    ungroup() %>%
    arrange(-dist) %>%
    mutate(patient = factor(patient, levels = patient))
  
  p_a <- make_figure_2a(fastbaps, metadata)
  p_b <- make_figure_2b(min_snp_dis, pal)
  p_c <- make_figure_2c(pal, tree, metadata, taxalinks)
  p_d <- make_figure_2d(pal, metadata, snp_dists_cluster, threshold)
  
  plot_list <- list(p_a, p_b, p_c, p_d)
  
  design <- "
ABC
ABC
DDC
DDC
DDC
"
  
  p <- patchwork::wrap_plots(plot_list) +
    patchwork::plot_layout(design = design, widths = c(1/3, 1, 1.5)) +
    patchwork::plot_annotation(tag_levels = "A")
  
  ggsave("plots/figure_2.png", p, width = 12.5, height = 7.5)
  
}
