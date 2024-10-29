library(tidyverse)

make_manhattan_plot <- function(gen_pos, p_adj) {
  
  p_inp <- gen_pos %>%
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
               show.legend = FALSE,
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
          axis.title.y = ggtext::element_markdown()) +
    labs(x = "Gene position (Mb)", y = "**&minus;log<sub>10</sub>(*p*)**")
  
}

viz_dn_ds <- function(gen_pos, amr_gen) {

  p_inp <- gen_pos %>%
    mutate(dn_ds = dn / ds) %>%
    arrange(desc(dn_ds)) %>%
    filter(dn_ds > 1 & !is.infinite(dn_ds)) %>%
    filter(p < 0.05) %>%
    filter(gene_id_prokka %in% amr_gen$gene) %>%
    mutate(group = case_when(
      bonferroni < 0.05 ~ "Bonferroni adjusted *p*-value < 0.05",
      fdr < 0.05 & bonferroni > 0.05 ~ "FDR adjusted *p*-value < 0.05",
      p < 0.05 & fdr > 0.05 ~ "Unadjusted *p*-value < 0.05"
    ))
  
  x_intercept_1 <- p_inp %>%
    filter(fdr < 0.05) %>%
    filter(p == max(p)) %>%
    pull(p)
  
  x_intercept_2 <- p_inp %>%
    filter(bonferroni < 0.05) %>%
    filter(p == max(p)) %>%
    pull(p)
  
  x_intercepts <- c(x_intercept_1, x_intercept_2) %>%
    purrr::map(function(x) -log10(x)) %>%
    unlist()
    
  ggplot(p_inp, aes(x = -log10(p), y = reorder(gene_id_prokka, -log10(p)))) +
    geom_rect(
      aes(
        xmin = min(-log10(p)),
        xmax = -log10(x_intercept_1),
        ymin = 0,
        ymax = Inf
      ),
      fill = colorspace::lighten("lightgrey", 0.25)
    ) +
    geom_rect(
      aes(
        xmin = -log10(x_intercept_1),
        xmax = -log10(x_intercept_2),
        ymin = 0,
        ymax = Inf
      ),
      fill = colorspace::lighten("lightgrey", 0.75)
    ) +
    geom_vline(xintercept = x_intercepts, linetype = "dashed", colour = "black") +
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
         size = "**Number of patients**",
         fill = "***d*<sub>N</sub>/*d*<sub>S</sub> ratio**") +
    scale_fill_distiller(palette = "Blues")
  
}

get_component_membership <- function(x) {
  
  igraph::components(x)$membership %>%
    data.frame() %>%
    rownames_to_column("name") %>%
    dplyr::rename(membership = 2)
  
}

make_network_plot <- function(sig_gen, network) {
  
  df <- sig_gen %>%
    select(gene_id_prokka, size) %>%
    setNames(c("name", "node_size"))
  
  gois <- unique(df$name)
  
  net <- network %>%
    select(3, 4, 6) %>%
    dplyr::rename(edge_size = 3)
  
  node_sizes <- net %>%
    select(1) %>%
    dplyr::rename(name = 1) %>%
    distinct() %>%
    left_join(df, by = "name") %>%
    mutate(node_colour = ifelse(name %in% gois, "1", "0")) %>%
    mutate(node_size = replace_na(node_size, 1)) %>%
    select(name, node_colour, node_size)
  
  node_numb <- node_sizes %>%
    select(name) %>%
    mutate(node_numb = row_number()) %>%
    as_tibble()
  
  g <- tidygraph::as_tbl_graph(net) %>%
    inner_join(node_sizes, by = "name") %>%
    tidygraph::activate(edges) %>%
    left_join(node_numb, by = c("from" = "node_numb")) %>%
    left_join(node_numb, by = c("to" = "node_numb")) %>%
    mutate(edge_colour = ifelse(name.x %in% gois & name.y %in% gois, "1", "0")) %>%
    select(from, to, edge_size, edge_colour) %>%
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
    ggraph::geom_edge_fan(aes(width = edge_size, colour = edge_colour),
                          alpha = 0.5,
                          show.legend = FALSE) +
    ggraph::scale_edge_color_manual(values = c("0" = "grey", "1" = "steelblue")) +
    ggraph::scale_edge_width(range = c(0.1, 1)) +
    scale_size(range = c(2, 6)) +
    ggraph::geom_node_point(aes(fill = node_colour, size = node_size),
                            pch = 21,
                            show.legend = FALSE) +
    ggraph::geom_node_text(aes(label = name), repel = TRUE, fontface = "italic") +
    theme(plot.background = element_blank(),
          panel.background = element_blank()) +
    scale_fill_manual(values = c("grey", "steelblue"))
  
}

make_figure_4a <- function(gen_pos, p_adj) {
  
  make_manhattan_plot(gen_pos, p_adj)
  
}

make_figure_4b <- function(sig_gen, network) {
  
  make_network_plot(sig_gen, network)

}

make_figure_4c <- function(gen_pos, amr_gen) {
  
  viz_dn_ds(gen_pos, amr_gen)
  
}

make_figure_4 <- function() {
  
  gen_pos <- read_delim("data/gen_pos.tsv")
  amr_gen <- read_delim("data/card_mycobacterial_args.tsv")
  network <- read_delim("data/network.limit_30.score_400.tsv")
  
  p1 <- make_figure_4a(gen_pos, "bonferroni")
  
  top <- cowplot::plot_grid(plot.new(), p1, plot.new(),
                            nrow = 1,
                            labels = c("", "A", ""),
                            rel_widths = c(1, 8, 1))
  
  sig_gen <- gen_pos %>%
    filter(bonferroni < 0.05)
  
  p2 <- make_figure_4b(sig_gen, network)
  
  p3 <- make_figure_4c(gen_pos, amr_gen)
  p3 <- cowplot::plot_grid(plot.new(), p3, plot.new(), nrow = 3, rel_heights = c(1, 8, 1))
  
  bottom <- cowplot::plot_grid(
    plotlist = list(p2, p3),
    nrow = 1,
    scale = 0.95,
    labels = c("B", "C"),
    rel_widths = c(1, 1.75)
  )
  
  p <- cowplot::plot_grid(
    top,
    bottom,
    nrow = 2,
    scale = 0.95,
    rel_heights = c(1, 1.25)
  )
  
  ggsave("plots/figure_4.png", p, width = 10, height = 10, bg = "white")
  
}
