library(tidyverse)

`%<+%` <- ggtree::`%<+%`

find_clades_with_shared_strains <- function(snp_dists, metadata, threshold) {
  
  strains <- snp_dists %>%
    filter(dist <= threshold) %>%
    select(genome_1, genome_2, dist)
  
  g <- strains %>%
    tidygraph::as_tbl_graph(directed = FALSE)
  
  components <- igraph::components(g)$membership %>%
    as.data.frame() %>%
    rownames_to_column("isolate") %>%
    rename(component = 2)
  
  keep <- components %>%
    inner_join(metadata, by = "isolate") %>%
    select(component, country) %>%
    distinct() %>%
    group_by(component) %>%
    tally() %>%
    ungroup() %>%
    filter(n > 1) %>%
    pull(component)
  
  comp2clus <- tibble(component = keep, cluster = paste("Cluster", 1:length(keep)))
  
  components %>%
    filter(component %in% keep) %>%
    inner_join(comp2clus, by = "component") %>%
    inner_join(metadata, by = "isolate") %>%
    arrange(cluster)
  
}

make_figure_3a <- function(fastbaps, metadata, pal) {
  
  dat <- inner_join(fastbaps, metadata, by = "isolate") %>%
    mutate(lineage = str_pad(level_1, 2, "left", "0"))
  
  p_inp <- dat %>%
    group_by(lineage, country) %>%
    tally() %>%
    ungroup()
  
  x_ord <- p_inp %>%
    group_by(lineage) %>%
    summarise(n = sum(n)) %>%
    arrange(desc(n)) %>%
    pull(lineage)
  
  p_inp$lineage <- factor(p_inp$lineage, levels = x_ord)
  
  p_inp$country <- factor(p_inp$country, levels = c("Ireland", "UK", "Germany"))
  
  ggplot(p_inp, aes(x = lineage, y = n)) +
    geom_bar(aes(fill = country), stat = "identity", colour = "black", show.legend = FALSE) +
    scale_fill_manual(values =  pal) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold")) +
    labs(x = "fastbaps lineage", y = "Number of genomes", fill = "Country")
  
}

make_figure_3b <- function(metadata, snp_dists, pal, threshold) {
  
  meta <- metadata %>%
    select(isolate, patient, country) %>%
    dplyr::rename(name = 1)
  
  snp_mat <- snp_dists %>%
    select(-1)
  
  identical_genomes <- snp_mat %>%
    filter(dist == 0 & genome_1 != genome_2) %>%
    tidygraph::as_tbl_graph(directed = FALSE) %>%
    igraph::cluster_leiden()
  
  non_identical_genomes <- meta %>%
    filter(!name %in% identical_genomes$names) %>%
    pull(name)
  
  genome_groups <- tibble(
    name = identical_genomes$names,
    cluster = identical_genomes$membership
  )
  
  identical_genome_representatives <- genome_groups %>%
    group_by(cluster) %>%
    sample_n(1) %>%
    pull(name)
  
  genomes <- c(non_identical_genomes, identical_genome_representatives)
  
  g <- snp_mat %>%
    filter(dist <= threshold) %>%
    tidygraph::as_tbl_graph()
  
  components <- igraph::components(g)$membership %>%
    as.data.frame() %>%
    dplyr::rename(component = 1) %>%
    rownames_to_column("name")
  
  keep <- components %>%
    inner_join(meta, by = "name") %>%
    select(component, country) %>%
    distinct() %>%
    group_by(component) %>%
    tally() %>%
    filter(n > 1) %>%
    pull(component)
  
  isolates <- components %>%
    filter(component %in% keep) %>%
    pull(name)
  
  g <- snp_mat %>%
    filter(genome_1 %in% genomes & genome_2 %in% genomes) %>%
    filter(dist <= threshold & genome_1 != genome_2) %>%
    tidygraph::as_tbl_graph() %>%
    igraph::mst() %>%
    tidygraph::as_tbl_graph() %>%
    filter(name %in% isolates) %>%
    inner_join(meta, by = "name")
  
  sizes <- inner_join(genome_groups, meta, by = "name") %>%
    select(cluster, patient, country) %>%
    distinct() %>%
    group_by(cluster, country) %>%
    tally() %>%
    ungroup()

  part_1 <- genome_groups %>%
    filter(name %in% pull(g, name)) %>%
    inner_join(sizes, by = "cluster") %>%
    select(name, country, n)
  
  part_2 <- g %>%
    as_tibble() %>%
    mutate(n = 1) %>%
    select(name, country, n) %>%
    filter(!name %in% part_1$name)

  part_3 <- bind_rows(part_1, part_2) %>%
    arrange(country) %>%
    pivot_wider(names_from = country, values_from = n, values_fill = 0) %>%
    mutate(size = Germany + Ireland + UK)
  
  g <- g %>%
    inner_join(part_3, by = "name")
  
  layout <- ggraph::create_layout(g, "fr")
  igraph::V(g)$x <- layout[, 1]$x
  igraph::V(g)$y <- layout[, 2]$y
  
  clades <- components %>%
    filter(name %in% pull(g, name)) %>%
    select(component) %>%
    distinct() %>%
    arrange(component) %>%
    mutate(clade = sprintf("Cluster %s", row_number()))
  
  components %>%
    filter(component %in% keep) %>%
    inner_join(clades, by = "component") %>%
    rename(isolate = 1) %>%
    write_tsv("data/international_clusters.tsv")
  
  ann <- components %>%
    filter(name %in% pull(g, name)) %>%
    inner_join(layout[, c("name", "x", "y")], by = "name") %>%
    group_by(component) %>%
    summarise(x = mean(x), y = max(y)) %>%
    ungroup() %>%
    inner_join(clades, by = "component")
  
  ggraph::ggraph(g, "manual", x = igraph::V(g)$x, y = igraph::V(g)$y) +
    ggraph::geom_edge_fan(aes(label = dist), colour = "lightgrey") +
    geom_text(data = ann,
              aes(x = x, y = y, label = clade),
              vjust = -2.5) +
    scatterpie::geom_scatterpie(data = layout,
                                aes(x = x, y = y, r = 0.2),
                                cols = c("Germany", "Ireland", "UK"),
                                colour = "black",
                                show.legend = FALSE) +
    scatterpie::geom_scatterpie(data = layout,
                                aes(x = x, y = y, r = 0.2),
                                cols = c("Germany", "Ireland", "UK"),
                                colour = NA,
                                show.legend = FALSE) +
    coord_fixed() +
    scale_fill_manual(values = pal) +
    theme(panel.background = element_blank(),
          title = element_text(face = "bold"))
  
}

make_figure_3c <- function(tree, metadata, fastbaps, snp_dists, pal, threshold) {
  
  clades <- read_delim("data/international_clusters.tsv") %>%
    select(isolate, clade)
  
  ann_1 <- metadata %>%
    select(isolate) %>%
    left_join(clades, by = "isolate") %>%
    rename(Cluster = 2)
  
  ann_2 <- metadata %>%
    select(isolate, country) %>%
    rename(Country = 2) %>%
    column_to_rownames("isolate")
  
  ann_3 <- fastbaps %>%
    select(1, 2) %>%
    dplyr::rename(isolate = 1, Lineage = 2) %>%
    mutate(Lineage = str_pad(Lineage, 2, "left", "0")) %>%
    column_to_rownames("isolate")
  
  p0 <- ggtree::ggtree(tree, layout = "circular")
  
  p1 <- p0 %<+%
    ann_1 +
    ggtree::geom_tippoint(
      aes(colour = Cluster, subset = !is.na(Cluster)),
      size = 3
    ) +
    scale_colour_brewer(palette = "Dark2") +
    theme(legend.title = element_text(face = "bold"))
  
  p2 <- ggtree::gheatmap(p1, ann_2, offset = 1e-6, width = .1, colnames = F, color = NA) +
    scale_fill_manual(values = pal) +
    theme(legend.title = element_text(face = "bold")) +
    labs(fill = "Country") +
    ggnewscale::new_scale_fill()
  
  p3 <- ggtree::gheatmap(p2, ann_3, offset = 1.75 * 1e-3, width = .1, colnames = F, color = NA) +
    scale_fill_viridis_d() +
    theme(legend.title = element_text(face = "bold")) +
    labs(fill = "fastbaps lineage")  +
    ggnewscale::new_scale_fill() +
    guides(fill = guide_legend(ncol = 3, order = 2))
  
  p3 +
    theme(
      legend.box = "horizontal",
      legend.position.inside = c(1.05, 0.8125)
    )
  
}

make_figure_3 <- function() {
  
  par(mar = c(1, 1, 1, 1))
  
  fastbaps  <- read_delim("data/fastbaps.tsv")
  metadata  <- read_delim("data/metadata.tsv")
  tree      <- ggtree::read.tree("data/iqtree.treefile")
  snp_dists <- read_delim("data/snp_dists.cluster_wide.tsv")
  threshold <- as.numeric(readLines("data/threshold.txt"))
  
  pal <- c(
    "Germany" = "#FFCC00",
    "Ireland" = "#FF883E",
    "UK" = "#012169"
  )
  
  p_a <- make_figure_3a(fastbaps, metadata, pal)
  p_b <- make_figure_3b(metadata, snp_dists, pal, threshold)
  p_c <- make_figure_3c(tree, metadata, fastbaps, snp_dists, pal, threshold)
  
  part_1 <- cowplot::plot_grid(plot.new(), p_a, plot.new(),
                               nrow = 1,
                               rel_widths = c(1, 2, 1),
                               labels = c("", "A", ""))
  
  part_2 <- cowplot::plot_grid(p_b, p_c, nrow = 1, labels = c("B", "C"), rel_widths = c(1, 1.5))
  
  figure_2 <- cowplot::plot_grid(part_1, part_2, nrow = 2, rel_heights = c(1, 1.75), scale = 0.95)
  
  ggsave("plots/figure_3.png", figure_2, width = 16.25, height = 10, bg = "white")
  
}
