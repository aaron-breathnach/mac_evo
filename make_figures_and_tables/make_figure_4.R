library(tidyverse)

make_manhattan_plot <- function(gen_pos) {
  
  ann <- gen_pos %>%
    filter(fdr <= 0.05) %>%
    select(gene_id_prokka, pos, p, fdr)
  
  y_intercept <- ann %>%
    filter(fdr <= 0.05) %>%
    filter(p == max(p)) %>%
    pull(p)
  
  ggplot(gen_pos, aes(x = pos, y = -log10(p))) +
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

viz_dn_ds <- function(gen_pos, amr_gen, make_bar = FALSE) {
  
  p_inp <- gen_pos %>%
    mutate(dn_ds = dn / ds) %>%
    arrange(desc(dn_ds)) %>%
    filter(dn_ds > 1 & !is.infinite(dn_ds)) %>%
    filter(SIZE >= 1)
  
  ggplot(p_inp, aes(x = dn_ds, y = reorder(gene_id_prokka, dn_ds))) +
    geom_point(aes(size = SIZE, fill = sig),
               pch = 21,
               colour = "black",
               show.legend = FALSE) +
    theme_classic(base_size = 12.5) +
    theme(axis.title.x = ggtext::element_markdown(),
          axis.title.y = element_text(face = "bold"),
          axis.text.y = element_text(face = "italic"),
          legend.title = element_text(face = "bold")) +
    labs(x = "***d*<sub>N</sub>/*d*<sub>S</sub> ratio**",
         y = "Gene") +
    scale_fill_manual(values = c("grey", "steelblue"))
  
}

get_component_membership <- function(x) {
  
  igraph::components(x)$membership %>%
    data.frame() %>%
    rownames_to_column("name") %>%
    dplyr::rename(membership = 2)
  
}

make_network_plot <- function(sig_gen, network) {
  
  df <- sig_gen[,1:2] %>%
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

make_figure_4a <- function(gen_pos) {
  
  make_manhattan_plot(gen_pos)
  
}

make_figure_4b <- function(sig_gen, network) {
  
  make_network_plot(sig_gen, network)

}

make_figure_4c <- function(gen_pos, amr_gen) {
  
  viz_dn_ds(gen_pos, amr_gen)
  
}

make_figure_4 <- function() {
  
  gen_pos <- read_delim("data/gen_pos.tsv")
  sig_gen <- read_delim("data/sig_gen.tsv")
  amr_gen <- read_delim("data/card_mycobacterial_args.tsv")
  network <- read_delim("data/network.limit_30.score_400.tsv")
  
  p1 <- make_figure_4a(gen_pos)
  
  top <- cowplot::plot_grid(plot.new(), p1, plot.new(),
                            nrow = 1,
                            labels = c("", "A", ""),
                            rel_widths = c(1, 8, 1))
  
  p2 <- make_figure_4b(sig_gen, network)
  
  p3 <- make_figure_4c(gen_pos, amr_gen)
  p3 <- cowplot::plot_grid(plot.new(), p3, plot.new(), nrow = 3, rel_heights = c(1, 8, 1))
  
  bottom <- cowplot::plot_grid(p2, p3, nrow = 1, scale = 0.95, labels = c("B", "C"))
  
  p <- cowplot::plot_grid(top, bottom,
                          nrow = 2,
                          scale = 0.95,
                          rel_heights = c(1, 1.25))
  
  ggsave("plots/figure_4.png", p, width = 10, height = 8.75, bg = "white")
  
}
