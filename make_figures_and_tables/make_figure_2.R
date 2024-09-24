library(tidyverse)

find_clades_with_shared_strains <- function(tree, snp_dists, metadata) {
  
  strains <- snp_dists %>%
    dplyr::rename(genome_1 = 1) %>%
    pivot_longer(!genome_1, names_to = "genome_2", values_to = "dist") %>%
    filter(genome_1 != genome_2) %>%
    filter(dist <= 20)
  
  clusters <- strains %>%
    tidygraph::as_tbl_graph(directed = FALSE) %>%
    igraph::cluster_leiden()
  
  clusters <- tibble(isolate = clusters$names, cluster = clusters$membership) 
  
  clusters_of_interest <- inner_join(clusters, metadata, by = "isolate") %>%
    select(cluster, study) %>%
    distinct() %>%
    group_by(cluster) %>%
    tally() %>%
    filter(n > 1) %>%
    pull(cluster)
  
  tips <- tibble(isolate = tree$tip.label, node = 1:length(tree$tip.label))
  
  max_dis <- clusters %>%
    filter(cluster %in% clusters_of_interest) %>%
    inner_join(strains[, c(1, 3)], by = c("isolate" = "genome_1")) %>%
    group_by(cluster) %>%
    filter(dist == max(dist)) %>%
    ungroup() %>%
    select(cluster, dist) %>%
    distinct() %>%
    mutate(dist = sprintf("%s%s SNPs", "\U2264", dist)) %>%
    arrange(cluster)
  
  clades <- clusters %>%
    filter(cluster %in% clusters_of_interest) %>%
    inner_join(tips, by = "isolate") %>%
    group_by(cluster) %>%
    summarise(nodes = paste0(node, collapse = ",")) %>%
    pull(nodes) %>%
    purrr::map(function(x) x %>% str_split(",") %>% unlist() %>% as.numeric())
  
  mrcas <- purrr::map(1:length(clades), function(x) ape::getMRCA(tree, clades[[x]])) %>%
    unlist()
  
  max_dis %>%
    mutate(mrca = mrcas)
  
}

plot_tree <- function(TREE, METADATA, FASTBAPS, SNP_DISTS, save = FALSE) {
  
  tree <- ggtree::read.tree(TREE) %>%
    ape::drop.tip("Reference")
  
  metadata <- read_delim(METADATA)
  
  snp_dists <- read_delim(SNP_DISTS)
  
  clades <- find_clades_with_shared_strains(tree, snp_dists, metadata)
  
  p0 <- ggtree::ggtree(tree, layout = "circular") +
    ggtree::geom_highlight(node = clades$mrca, extendto = 1.5 * max(tree$edge.length))
  
  ann_1 <- metadata %>%
    dplyr::rename(Study = study) %>%
    select(isolate, Study) %>%
    column_to_rownames("isolate")
  
  p1 <- ggtree::gheatmap(p0, ann_1, offset = 1, width = .1, colnames = F, color = NA) +
    scale_fill_manual(values = c("#ed4037", "#262161", "#faaf40")) +
    theme(legend.title = element_text(face = "bold")) +
    labs(fill = "Study") +
    ggnewscale::new_scale_fill()
  
  ann_2 <- read_delim(FASTBAPS) %>%
    select(1, 2) %>%
    dplyr::rename(isolate = 1, Lineage = 2) %>%
    mutate(Lineage = str_pad(Lineage, 2, "left", "0")) %>%
    column_to_rownames("isolate")
  
  ggtree::gheatmap(p1, ann_2, offset = .375 * 1e4, width = .1, colnames = F, color = NA) +
    scale_fill_viridis_d() +
    theme(legend.title = element_text(face = "bold")) +
    labs(fill = "fastbaps lineage")  +
    ggnewscale::new_scale_fill()
  
}

viz_putative_transmission_clusters <- function(METADATA, SNP_DISTS) {
  
  metadata <- read_delim(METADATA)
  
  snp_mat <- read_delim(SNP_DISTS) %>%
    rename(genome_1 = 1) %>%
    filter(genome_1 != "Reference") %>%
    select(-any_of("Reference")) %>%
    pivot_longer(!genome_1, names_to = "genome_2", values_to = "dist")
  
  identical_genomes <- snp_mat %>%
    filter(dist == 0 & genome_1 != genome_2)
  
  genome_groups <- identical_genomes %>%
    tidygraph::as_tbl_graph(directed = FALSE) %>%
    igraph::cluster_leiden()
  
  non_duplicates <- metadata %>%
    filter(!isolate %in% genome_groups$names) %>%
    pull(isolate)
  
  reps <- tibble(isolate = genome_groups$names, group = genome_groups$membership) %>%
    group_by(group) %>%
    sample_n(1) %>%
    pull(isolate)
  
  genomes <- c(non_duplicates, reps)
  
  g <- snp_mat %>%
    filter(dist <= 20) %>%
    tidygraph::as_tbl_graph()
  
  components <- igraph::components(g)$membership %>%
    as.data.frame() %>%
    rename(component = 1) %>%
    rownames_to_column("isolate")
  
  keep <- components %>%
    inner_join(metadata, by = "isolate") %>%
    select(component, study) %>%
    distinct() %>%
    group_by(component) %>%
    tally() %>%
    filter(n > 1) %>%
    pull(component)
  
  isolates <- components %>%
    filter(component %in% keep) %>%
    pull(isolate)
  
  components %>%
    filter(component %in% keep) %>%
    inner_join(metadata, by = "isolate") %>%
    arrange(study) %>%
    arrange(component)
  
  g <- snp_mat %>%
    filter(genome_1 %in% genomes & genome_2 %in% genomes) %>%
    filter(dist <= 20 & genome_1 != genome_2) %>%
    tidygraph::as_tbl_graph() %>%
    igraph::mst() %>%
    tidygraph::as_tbl_graph() %>%
    filter(name %in% isolates) %>%
    inner_join(metadata, by = c("name" = "isolate"))
  
  ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_fan(aes(label = dist), alpha = 0.25) +
    ggraph::geom_node_point(aes(fill = study), pch = 21, size = 5) +
    theme(panel.background = element_blank(),
          title = element_text(face = "bold")) +
    scale_fill_manual(values = c("#ed4037", "#262161", "#faaf40")) +
    labs(fill = "Study")
  
}

figure_2a <- function(FASTBAPS, METADATA) {

  fastbaps <- read_delim(FASTBAPS)
  metadata <- read_delim(METADATA)
  dat <- inner_join(fastbaps, metadata, by = "isolate") %>%
    mutate(lineage = str_pad(level_1, 2, "left", "0"))
  
  p_inp <- dat %>%
    group_by(lineage, study) %>%
    tally() %>%
    ungroup()
  
  x_ord <- p_inp %>%
    group_by(lineage) %>%
    summarise(n = sum(n)) %>%
    arrange(desc(n)) %>%
    pull(lineage)
  
  p_inp$lineage <- factor(p_inp$lineage, levels = x_ord)
  
  pal <- c("Present" = "#ed4037", "Van Tonder" = "#262161", "Wetzstein" = "#faaf40")
  
  ggplot(p_inp, aes(x = lineage, y = n)) +
    geom_bar(aes(fill = study), stat = "identity", colour = "black", show.legend = FALSE) +
    scale_fill_manual(values =  pal) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold")) +
    labs(x = "fastbaps lineage", y = "Number of genomes", fill = "Study")

}

figure_2b <- function(TREE, METADATA, FASTBAPS, SNP_DISTS) {
  
  plot_tree(TREE, METADATA, FASTBAPS, SNP_DISTS)
  
}

figure_2c <- function(METADATA, SNP_DISTS) {
  
  viz_putative_transmission_clusters(METADATA, SNP_DISTS)
  
}

make_figure_2 <- function() {
  
  fastbaps  <- "data/fastbaps.tsv"
  metadata  <- "data/metadata.tsv"
  tree      <- "data/mavium.final_tree.tre"
  snp_dists <- "data/snp_dists.tsv"å
  
  p_a <- figure_2a(fastbaps, metadata)
  p_b <- figure_2b(tree, metadata, fastbaps, snp_dists)
  p_c <- figure_2c(metadata, snp_dists)
  
  part_1 <- cowplot::plot_grid(plot.new(), p_a, plot.new(),
                               nrow = 1,
                               rel_widths = c(1, 2, 1),
                               labels = c("", "A", ""))
  
  part_2 <- cowplot::plot_grid(p_b, p_c, nrow = 1, labels = c("B", "C"))
  
  figure_2 <- cowplot::plot_grid(part_1, part_2, nrow = 2, rel_heights = c(1, 1.5), scale = 0.95)
  
  ggsave("plots/figure_2.png", figure_2, width = 15, height = 10, bg = "white")
  
}
