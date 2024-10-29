library(tidyverse)

`%<+%` <- ggtree::`%<+%`

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

make_umap <- function(snp_dists, pids, metadata, pal) {
  
  fastbaps <- read_delim("data/fastbaps.tsv") %>%
    mutate(lineage = str_pad(level_1, 2, "left", "0"))
  
  inp <- snp_dists %>%
    rename(isolate = 1) %>%
    column_to_rownames("isolate")
  
  p_inp <- run_umap(inp, metadata) %>%
    inner_join(fastbaps, by = "isolate")
  
  ggplot(p_inp, aes(x = UMAP1, y = UMAP2)) +
    geom_point(
      aes(fill = lineage),
      pch = 21,
      colour = "black",
      size = 3
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    ) +
    labs(fill = "fastbaps lineage") +
    scale_fill_brewer(palette = "Set3")
  
}

get_jumps <- function(pid, metadata, snp_dists) {
  
  meta <- metadata %>%
    filter(patient == pid) %>%
    select(isolate, time_from_diagnosis)
  
  n <- nrow(meta)
  
  pairs <- 2:n %>%
    purrr::map(function(x) tibble(gen_1 = meta[[x - 1, 1]], gen_2 = meta[[x, 1]])) %>%
    bind_rows()
  
  if (nrow(pairs) > 1) {
    
    out <- snp_dists %>%
      filter(genome_1 %in% meta$isolate) %>%
      select(1, all_of(meta$isolate)) %>%
      pivot_longer(!genome_1, names_to = "genome_2", values_to = "dist") %>%
      filter(genome_1 != genome_2) %>%
      inner_join(meta, by = c("genome_1" = "isolate")) %>%
      inner_join(meta, by = c("genome_2" = "isolate")) %>%
      mutate(gen_1 = case_when(
        time_from_diagnosis.x < time_from_diagnosis.y ~ genome_1,
        time_from_diagnosis.y < time_from_diagnosis.x ~ genome_2
      )) %>%
      mutate(gen_2 = case_when(
        time_from_diagnosis.x > time_from_diagnosis.y ~ genome_1,
        time_from_diagnosis.y > time_from_diagnosis.x ~ genome_2
      )) %>%
      select(gen_1, gen_2, dist) %>%
      distinct() %>%
      inner_join(pairs, by = c("gen_1", "gen_2")) %>%
      filter(dist > 12) %>%
      mutate(patient = pid) %>%
      select(4, 1, 2)
    
  } else {
    
    out <- pairs %>%
      mutate(patient = pid) %>%
      select(3, 1, 2)
    
  }
  
  return(out)
  
}

make_fastbaps_bar_chart <- function() {
  
  fastbaps <- read_delim("data/fastbaps.tsv")
  
  metadata <- read_delim("data/metadata.tsv") %>%
    filter(study == "Present")
  
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
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Num. lineages", y = "Num. patients")
  
}

make_tree <- function(pal) {
  
  metadata <- read_delim("data/metadata.tsv") %>%
    filter(study == "Present" & bracken_pass)
  
  snp_dists <- read_delim("data/snp_dists.tsv")
  
  tree <- ggtree::read.tree("data/mavium.final_tree.tre")
  
  isolates <- unlist(lapply(tree$tip.label, function(x) if (grepl("^2", x)) x))
  
  tree <- ape::keep.tip(tree, isolates)
  
  patients <- metadata %>%
    filter(multiple_carriage) %>%
    pull(patient) %>%
    unique()
  
  meta <- metadata %>%
    filter(study == "Present" & patient %in% patients & isolate %in% isolates)
  
  snp_dists <- read_delim("data/snp_dists.tsv") %>%
    rename(genome_1 = 1)
  
  pids <- unique(meta$patient)
  
  df_1 <- metadata %>%
    mutate(patient = ifelse(patient %in% pids, patient, NA)) %>%
    select(isolate, patient)
  
  df_2 <- purrr::map(pids, function(x) get_jumps(x, meta, snp_dists)) %>%
    bind_rows()
  
  ggtree::ggtree(tree) %<+%
    df_1 +
    ggtree::geom_taxalink(
      data = df_2,
      aes(taxa1 = gen_1, taxa2 = gen_2),
      colour = "black",
      size = 1,
    ) +
    ggtree::geom_taxalink(
      data = df_2,
      aes(taxa1 = gen_1, taxa2 = gen_2, colour = patient),
      show.legend = FALSE
    ) +
    ggtree::geom_tippoint(
      aes(fill = patient, subset = !is.na(patient)),
      pch = 21,
      size = 2
    ) +
    scale_size_continuous(range = c(0, 2)) +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = pal) +
    labs(fill = "Patient") +
    theme(legend.title = element_text(face = "bold"))
  
}

make_subspecies_bar_chart <- function() {
  
  metadata <- read_delim("data/metadata.tsv") %>%
    filter(bracken_perc > 70 & study == "Present")
  
  ref_gen <- read_delim("data/reference_genomes.tsv") %>%
    select(1, 3) %>%
    rename(subspecies = 1, reference = 2) %>%
    mutate(subspecies = recode(
      subspecies,
      "Mycobacterium avium subsp. avium" = "*M. avium* subsp. *avium*",
      "Mycobacterium avium subsp. hominissuis" = "*M. avium* subsp. *hominissuis*"
    )) %>%
    mutate(reference = reference %>%
             basename() %>%
             str_replace(".fna.*", ""))
  
  cols <- c("query", metadata$isolate, ref_gen$reference)
  
  mash_dists <- read_delim("data/mash.dist.txt") %>%
    rename(query = 1) %>%
    setNames(str_replace(basename(names(.)), ".fna.*", "")) %>%
    mutate(query = str_replace(basename(query), ".fna.*", "")) %>%
    filter(query %in% metadata$isolate) %>%
    select(query, all_of(cols))
  
  subspecies <- mash_dists %>%
    rename(query = 1) %>%
    filter(grepl("^2", query)) %>%
    pivot_longer(!query, names_to = "reference", values_to = "dist") %>%
    inner_join(ref_gen, by = "reference") %>%
    group_by(query, subspecies) %>%
    summarise(dist = mean(dist)) %>%
    ungroup() %>%
    group_by(query) %>%
    filter(dist == min(dist)) %>%
    ungroup()
  
  p_inp_1 <- subspecies %>%
    group_by(subspecies) %>%
    tally() %>%
    ungroup() %>%
    mutate(abbreviation = ifelse(grepl("hominissuis", subspecies), "MAH", "MAA"))
  
  ggplot(p_inp_1, aes(x = abbreviation, y = n)) +
    geom_bar(stat = "identity", colour = "black", fill = "steelblue") +
    geom_text(aes(label = paste0("n=", n)), vjust = -0.25) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold")) +
    labs(x = "Subspecies", y = "Num. isolates") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
}

viz_putative_transmission_clusters <- function(pal) {
  
  metadata <- read_delim("data/metadata.tsv") %>%
    filter(study == "Present" & bracken_pass) %>%
    select(isolate, patient) %>%
    rename(name = 1)
  
  snp_mat <- read_delim("data/snp_dists.tsv") %>%
    rename(genome_1 = 1) %>%
    filter(genome_1 != "Reference") %>%
    select(-any_of("Reference")) %>%
    pivot_longer(!genome_1, names_to = "genome_2", values_to = "dist") %>%
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
  
  # dup_rep <- genome_groups %>%
  #   group_by(cluster) %>%
  #   sample_n(1) %>%
  #   pull(name)
  
  g <- snp_mat %>%
    filter(dist <= 12) %>%
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
    filter(dist <= 12 & genome_1 != genome_2) %>%
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
  
  layout <- ggraph::create_layout(g, "fr")
  igraph::V(g)$x <- layout[, 1]
  igraph::V(g)$y <- layout[, 2]
  
  clades <- components %>%
    filter(name %in% pull(g, name)) %>%
    select(component) %>%
    distinct() %>%
    arrange(component)
  
  ann <- components %>%
    filter(name %in% pull(g, name)) %>%
    inner_join(layout[, c("name", "x", "y")], by = "name") %>%
    group_by(component) %>%
    summarise(x = mean(x), y = max(y)) %>%
    ungroup() %>%
    inner_join(clades, by = "component") %>%
    arrange(desc(y)) %>%
    mutate(clade = paste0("Cluster ", row_number()))
  
  ggraph::ggraph(g, "manual", x = igraph::V(g)$x, y = igraph::V(g)$y) +
    ggraph::geom_edge_fan(aes(label = dist), colour = "lightgrey") +
    geom_text(data = ann,
              aes(x = x, y = y, label = clade),
              vjust = -1.25) +
    scatterpie::geom_scatterpie(data = layout,
                                aes(x = x, y = y, r = 0.251 * size / max(layout$size)),
                                cols = colnames(df_b)[-1],
                                colour = "black") +
    scatterpie::geom_scatterpie(data = layout,
                                aes(x = x, y = y, r = 0.25 * size / max(layout$size)),
                                cols = colnames(df_b)[-1],
                                colour = NA) +
    coord_fixed() +
    scale_fill_manual(values = pal) +
    theme(panel.background = element_blank(),
          title = element_text(face = "bold")) +
    labs(fill = "Patient") +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.125)))
  
}

make_figure_1 <- function() {
  
  set_1 <- c("P02", "P04", "P06", "P07", "P08", "P12", "P16", "P19")
  set_2 <- c("P02", "P04", "P06", "P10", "P19")
  set <- sort(unique(c(set_1, set_2)))
  pal <- RColorBrewer::brewer.pal(length(set), "Set3")
  names(pal) <- set
  
  p0 <- make_subspecies_bar_chart()
  p1 <- make_fastbaps_bar_chart()
  p2 <- make_tree(pal)
  p3 <- viz_putative_transmission_clusters(pal)
  
  plot_list <- list(p0, p1, p2, p3)
  
  design <- "
ABC
ABC
DDC
DDC
DDC
"
  
  p <- patchwork::wrap_plots(plot_list) +
    patchwork::plot_layout(design = design, widths = c(0.5, 0.5, 1.5)) +
    patchwork::plot_annotation(tag_levels = "A")
  
  ggsave("plots/figure_1.png", p, width = 10, height = 5)
  
}

make_figure_1()
