library(tidyverse)

make_figure_s1 <- function() {
  
  mash_dists <- read_delim("data/mash.dist.txt") %>%
    rename(query = 1) %>%
    setNames(str_replace(basename(names(.)), ".fna.*", "")) %>%
    mutate(query = str_replace(basename(query), ".fna.*", ""))
  
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
  
  p1 <- ggplot(p_inp_1, aes(x = abbreviation, y = n)) +
    geom_bar(aes(fill = subspecies),
             stat = "identity",
             colour = "black",
             show.legend = FALSE) +
    geom_text(aes(label = n), vjust = -0.25) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold")) +
    labs(x = "Subspecies", y = "Number of isolates") +
    scale_fill_manual(values = c("#869F77", "#E7B5AC")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  ##########
  ## umap ##
  ##########
  
  mash_dists <- mash_dists %>%
    column_to_rownames("query")
  
  mash_dists <- mash_dists[
    rownames(mash_dists) %in% subspecies$query,
    colnames(mash_dists) %in% subspecies$query
  ]
  
  set.seed(321)
  snps_umap <- umap::umap(mash_dists)
  
  points <- snps_umap$layout %>%
    as.data.frame() %>%
    rownames_to_column("isolate") %>%
    rename(UMAP1 = 2, UMAP2 = 3) %>%
    filter(isolate != "Reference")
  
  p_inp_2 <- inner_join(points, subspecies, by = c("isolate" = "query")) %>%
    mutate(subspecies = recode(
      subspecies,
      "Mycobacterium avium subsp. avium" = "*M. avium* subsp. *avium*",
      "Mycobacterium avium subsp. hominissuis" = "*M. avium* subsp. *hominissuis*"
    ))
  
  p2 <- ggplot(p_inp_2, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(fill = subspecies), pch = 21, size = 2) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          legend.text = ggtext::element_markdown()) +
    labs(fill = "Subspecies") +
    scale_fill_manual(values = c("#869F77", "#E7B5AC"))
  p2
  
  #######################
  ## combine the plots ##
  #######################
  
  p <- cowplot::plot_grid(p1, p2,
                          labels = "AUTO",
                          rel_widths = c(1, 3))
  
  ggsave("plots/figure_s1.png", p, width = 10, height = 4)
  
}
