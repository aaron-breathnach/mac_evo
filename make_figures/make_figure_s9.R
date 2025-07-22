library(tidyverse)

make_figure_s9 <- function() {
  
  metadata <- read_delim("data/metadata.tsv")
  
  rpoB <- read_delim("data/panaroo/gene_data.csv") %>%
    filter(grepl("rpoB", gene_name))
  
  variants <- read_delim("data/filtered_variants.tsv")
  
  variants$annotation_id <- variants$INFO %>%
    purrr::map(function(x) {
      x %>%
        str_split("\\|") %>%
        unlist() %>%
        nth(5)
    }) %>%
    unlist()
  
  variants <- variants %>%
    inner_join(rpoB, by = "annotation_id") %>%
    select(GENOME, INFO, annotation_id) %>%
    inner_join(metadata, by = c("GENOME" = "isolate")) %>%
    select(patient, INFO, annotation_id) %>%
    unique() %>%
    filter(grepl("missense", INFO))
  
  variants$pos <- variants$INFO %>%
    purrr::map(function(x) {
      x %>%
        str_split("\\|") %>%
        unlist() %>%
        nth(12) %>%
        str_replace("/.*", "") %>%
        as.numeric()
    }) %>%
    unlist()
  
  gen_len <- 3435
  
  tmp <- tibble(pos = 1:gen_len)
  
  p_inp <- variants %>%
    select(patient, pos) %>%
    as.data.frame() %>%
    arrange(pos) %>%
    group_by(pos) %>%
    tally() %>%
    ungroup()
  
  p_inp <- tmp %>%
    left_join(p_inp, by = "pos") %>%
    mutate(n = ifelse(is.na(n), 0, n))
  
  ann <- p_inp %>%
    filter(n > 1)
  
  p <- ggplot(p_inp, aes(x = pos, y = n)) +
    geom_rect(
      aes(
        xmin = 1521,
        xmax = 1599,
        ymin = 0,
        ymax = Inf
      ),
      fill = colorspace::lighten("grey", 0.25)
    ) +
    geom_line(colour = "steelblue") +
    geom_text(data = ann,
              aes(x = pos, y = 1.05 * n, label = scales::comma(pos))) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(breaks = c(1:max(p_inp$n))) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold")) +
    labs(x = "Nucleotide position", y = "Num. patients")
  
  ggsave("plots/figure_s9.png", p, width = 10, height = 3.75)
  
}
