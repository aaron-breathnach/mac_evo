library(tidyverse)

make_upset_plot <- function(dat) {
  
  pal <- c(
    "#ffa500", #1
    "#800080", #2
    "#a52929", #3
    "#0000ff", #4
    "#008000"  #5
  )
  
  queries <- list(
    ComplexUpset::upset_query(
      intersect = c("Ethambutol", "Rifampicin", "Clarithromycin"),
      color = pal[1],
      fill = pal[1],
      only_components = c("intersections_matrix", "Number of patients")
    ),
    ComplexUpset::upset_query(
      intersect = c("Ethambutol", "Clarithromycin"),
      color = pal[2],
      fill = pal[2],
      only_components = c("intersections_matrix", "Number of patients")
    ),
    ComplexUpset::upset_query(
      intersect = c("Ethambutol", "Rifampicin", "Clarithromycin", "Amikacin", "Levofloxacin"),
      color = pal[3],
      fill = pal[3],
      only_components = c("intersections_matrix", "Number of patients")
    ),
    ComplexUpset::upset_query(
      intersect = c("Ethambutol", "Azithromycin"),
      color = pal[4],
      fill = pal[4],
      only_components = c("intersections_matrix", "Number of patients")
    ),
    ComplexUpset::upset_query(
      intersect = c("Ethambutol", "Rifampicin"),
      color = pal[5],
      fill = pal[5],
      only_components = c("intersections_matrix", "Number of patients")
    )
  )
  
  base_annotations <- list(
    "Number of patients" = ComplexUpset::intersection_size(counts = FALSE)
  )
  
  ComplexUpset::upset(
    dat,
    intersect = colnames(dat),
    base_annotations = base_annotations,
    queries = queries,
    set_sizes = FALSE
  )
  
}

make_figure_1a <- function(patient_metadata) {
  
  p_inp <- patient_metadata %>%
    select(1:7) %>%
    pivot_longer(!c(patient, Gender)) %>%
    group_by(Gender, name, value) %>%
    tally() %>%
    ungroup() %>%
    arrange(desc(n)) %>%
    mutate(name = factor(name, levels = c("Nationality",
                                          "Bronchiectasis",
                                          "HIV",
                                          "Smoker",
                                          "BMI"))) %>%
    group_by(name) %>%
    mutate(perc = 100 * n / sum(n)) %>%
    ungroup() %>%
    mutate(label = paste0(round(perc, 2), "%"))
  
  ggplot(p_inp, aes(x = n, y = value)) +
    facet_grid(name ~ ., scales = "free", space = "free", switch = "y") +
    geom_bar(aes(fill = Gender), stat = "identity", colour = "black") +
    theme_bw(base_size = 12.5) +
    theme(panel.grid = element_blank(),
          strip.text.y.left = element_text(face = "bold", angle = 0),
          axis.title.x = element_text(face = "bold"),
          axis.title.y = element_blank(),
          legend.title = element_text(face = "bold")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    scale_y_discrete(position = "right") +
    labs(x = "Number", fill = "Sex") +
    scale_fill_manual(values = c("pink", "blue"))
  
}

make_figure_1b <- function(patient_metadata) {
  
  dat <- patient_metadata %>%
    select(patient, Treated, Treatment) %>%
    setNames(str_to_lower(names(.))) %>%
    filter(treated == "Yes") %>%
    separate_longer_delim(treatment, ", ") %>%
    mutate(treated = recode(treated, "Yes" = 1)) %>%
    pivot_wider(names_from = treatment, values_from = treated, values_fill = 0) %>%
    column_to_rownames("patient")
  
  make_upset_plot(dat)
  
}

make_figure_1 <- function() {
  
  patient_metadata <- read_delim("data/patient_info.tsv")
  
  p_a <- make_figure_1a(patient_metadata)
  p_b <- make_figure_1b(patient_metadata)
  p_ab <- cowplot::plot_grid(p_a, p_b, labels = c("A", "B"), scale = 0.95)
  
  img <- png::readPNG("data/fig_1_c.png")

  p_c <- ggplot2::ggplot() +
    ggplot2::annotation_custom(
      grid::rasterGrob(
        img,
        width = ggplot2::unit(1, "npc"),
        height = ggplot2::unit(1, "npc")
      ),
      -Inf, Inf, -Inf, Inf
    )
  
  p <- cowplot::plot_grid(p_ab, p_c,
                          nrow = 2,
                          rel_heights = c(1, 1.5),
                          labels = c("", "C"))
  
  ggsave("plots/figure_1.png", p, width = 11.25, height = 12.5, bg = "white")
  
}
