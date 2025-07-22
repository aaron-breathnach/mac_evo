library(tidyverse)

get_component_membership <- function(x) {
  
  igraph::components(x)$membership %>%
    data.frame() %>%
    rownames_to_column("name") %>%
    dplyr::rename(membership = 2)
  
}

make_figure_5a <- function(mut_dis, p_adj) {
  
  p_inp <- mut_dis %>%
    select(gene_id_prokka, pos, p, size, all_of(p_adj)) %>%
    dplyr::rename(p_adj = 5) %>%
    mutate(sig = ifelse(p_adj < 0.05, "y", "n"))
  
  ann <- p_inp %>%
    filter(p_adj <= 0.05) %>%
    select(gene_id_prokka, pos, p, p_adj)
  
  y_intercept <- ann %>%
    filter(p_adj <= 0.05) %>%
    filter(p == max(p)) %>%
    pull(p)
  
  ggplot(p_inp, aes(x = pos, y = -log10(p))) +
    geom_hline(yintercept = -log10(y_intercept),
               colour = "grey",
               linetype = "dashed") +
    geom_point(aes(colour = sig, fill = sig, size = size),
               pch = 21,
               # show.legend = FALSE,
               alpha = 0.75) +
    scale_colour_manual(values = c("grey", "black")) +
    scale_fill_manual(values = c("grey", "steelblue")) +
    ggrepel::geom_text_repel(data = ann,
                             aes(label = gene_id_prokka),
                             min.segment.length = 0,
                             max.overlaps = Inf,
                             fontface = "italic") +
    theme_classic(base_size = 12.5) +
    theme(axis.title.x = element_text(face = "bold"),
          axis.title.y = ggtext::element_markdown(),
          legend.title = element_text(face = "bold")) +
    labs(x = "Gene position (Mb)",
         y = "**&minus;log<sub>10</sub>(*p*)**",
         size = "Num. patients") +
    guides(colour = "none", fill = "none")
  
}

make_figure_5b <- function(sig_gen, network) {
  
  df <- sig_gen %>%
    select(gene_id_prokka, size) %>%
    setNames(c("name", "node_size"))
  
  gois <- unique(df$name)
  
  net <- network %>%
    select(3, 4, 6) %>%
    setNames(c("gene_1", "gene_2", "edge_size"))
  
  node_sizes <- tibble(name = unique(c(net$gene_1, net$gene_2))) %>%
    left_join(df, by = "name") %>%
    mutate(node_colour = ifelse(name %in% gois, "1", "0")) %>%
    mutate(node_size = replace_na(node_size, 1)) %>%
    select(name, node_colour, node_size)
  
  node_numb <- node_sizes %>%
    select(name) %>%
    mutate(node_numb = row_number()) %>%
    as_tibble()
  
  g <- tidygraph::as_tbl_graph(net, directed = FALSE) %>%
    inner_join(node_sizes, by = "name") %>%
    tidygraph::activate(edges) %>%
    left_join(node_numb, by = c("from" = "node_numb")) %>%
    left_join(node_numb, by = c("to" = "node_numb")) %>%
    select(from, to, edge_size) %>%
    tidygraph::activate(nodes)
  
  clusters <- get_component_membership(g) %>%
    mutate(goi = ifelse(name %in% gois, 1, 0)) %>%
    group_by(membership) %>%
    summarise(n = sum(goi)) %>%
    filter(n > 1) %>%
    pull(membership)
  
  genes <- get_component_membership(g) %>%
    filter(membership %in% clusters) %>%
    pull(name)
  
  g <- g %>%
    filter(name %in% genes)
  
  ggraph::ggraph(g, layout = "kk") + 
    ggraph::geom_edge_fan(aes(width = edge_size),
                          colour = "grey",
                          show.legend = FALSE) +
    ggraph::scale_edge_width(range = c(0.25, 1)) +
    scale_size(range = c(2, 6)) +
    ggraph::geom_node_point(aes(fill = node_colour, size = node_size),
                            pch = 21,
                            show.legend = FALSE) +
    ggraph::geom_node_text(aes(label = name), repel = TRUE, fontface = "italic") +
    theme(plot.background = element_blank(),
          panel.background = element_blank()) +
    scale_fill_manual(values = c("grey", "steelblue"))
  
}

make_figure_5c <- function(mut_dis, amr_gen) {

  p_inp <- mut_dis %>%
    mutate(dn_ds = dn / ds) %>%
    arrange(desc(dn_ds)) %>%
    filter(dn_ds > 1 & !is.infinite(dn_ds)) %>%
    filter(p < 0.05) %>%
    filter(gene_id_prokka %in% amr_gen$gene) %>%
    group_by(gene_id_prokka) %>%
    filter(p == min(p)) %>%
    ungroup() %>%
    mutate(group = case_when(
      bonferroni < 0.05 ~ "Bonferroni adjusted *p*-value < 0.05",
      fdr < 0.05 & bonferroni > 0.05 ~ "FDR adjusted *p*-value < 0.05",
      p < 0.05 & fdr > 0.05 ~ "Unadjusted *p*-value < 0.05"
    ))
    
  ggplot(p_inp, aes(x = -log10(p), y = reorder(gene_id_prokka, -log10(p)))) +
    geom_point(aes(size = size, fill = dn_ds),
               pch = 21,
               colour = "black",
               show.legend = TRUE) +
    theme_classic(base_size = 12.5) +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = element_text(face = "bold"),
          axis.text.y = element_text(face = "italic"),
          legend.title = ggtext::element_markdown()) +
    labs(x = "**-log<sub>10</sub>(*p*-value)**",
         y = "Gene",
         size = "**Num. patients**",
         fill = "***d*<sub>N</sub>/*d*<sub>S</sub> ratio**") +
    scale_fill_distiller(palette = "Blues")
  
}

make_figure_5 <- function() {
  
  par(mar = c(1, 1, 1, 1))
  
  mut_dis <- read_delim("data/mutational_distribution.tsv")
  amr_gen <- read_delim("data/card_mycobacterial_args.tsv")
  network <- read_delim("data/network.limit_30.score_400.tsv")
  
  p1 <- make_figure_5a(mut_dis, "bonferroni")
  
  top <- cowplot::plot_grid(plot.new(), p1, plot.new(),
                            nrow = 1,
                            labels = c("", "A", ""),
                            rel_widths = c(1, 8, 1))
  
  sig_gen <- mut_dis %>%
    filter(bonferroni < 0.05)
  
  p2 <- make_figure_5b(sig_gen, network)
  
  p3 <- make_figure_5c(mut_dis, amr_gen)
  p3 <- cowplot::plot_grid(plot.new(), p3, plot.new(), nrow = 3, rel_heights = c(1, 8, 1))
  
  bottom <- cowplot::plot_grid(
    plotlist = list(p2, p3),
    nrow = 1,
    scale = 0.95,
    labels = c("B", "C"),
    rel_widths = c(1, 1.25)
  )
  
  p <- cowplot::plot_grid(
    top,
    bottom,
    nrow = 2,
    scale = 0.95,
    rel_heights = c(1, 1.25)
  )
  
  ggsave("plots/figure_5.png", p, width = 12.5, height = 10, bg = "white")
  
}
