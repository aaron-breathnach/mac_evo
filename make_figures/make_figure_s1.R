library(tidyverse)

make_boxplot <- function(dat, variable, y_lab) {
  
  p_inp <- dat %>%
    select(study, patient, all_of(variable)) %>%
    rename(n = 3)
  
  ggplot(p_inp, aes(x = study, y = n)) +
    geom_jitter(
      aes(colour = study, fill = study),
      alpha = 0.5,
      show.legend = FALSE
    ) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold")) +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = pal) +
    labs(x = "Study", y = y_lab)
  
}

make_figure_s1 <- function() {

  pal <- c("#FF883E", "#012169", "#FFCC00")
  
  #############
  ## panel a ##
  #############
  
  dat_a <- tibble(
    study = c("Present", "Van Tonder", "Wetzstein"),
    n = c(50, 135, 102)
  )
  
  p_a <- ggplot(dat_a, aes(x = study, y = n)) +
    geom_bar(
      aes(fill = study),
      stat = "identity",
      colour = "black",
      width = 0.75,
      show.legend = FALSE
    ) +
    geom_text(aes(label = paste0("n=", n)), vjust = -0.25) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold")) +
    scale_fill_manual(values = pal) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Study", y = "Number of isolates per study")
  
  #############
  ## panel b ##
  #############
  
  metadata <- read_delim("data/metadata.tsv")
  
  dat_b <- metadata %>%
    group_by(patient) %>%
    filter(time_from_diagnosis == max(time_from_diagnosis)) %>%
    ungroup()
  
  p_b <- make_boxplot(
    dat = dat_b,
    variable = "time_from_diagnosis",
    y_lab = "Max time from diagnosis (years) per patient"
  )
  
  #############
  ## panel c ##
  #############
  
  dat_c <- metadata %>%
    group_by(study, patient) %>%
    summarise(number_of_isolates = n()) %>%
    ungroup()
  
  p_c <- make_boxplot(
    dat = dat_c,
    variable = "number_of_isolates",
    y_lab = "Number of isolates per patient"
  )
  
  ########################
  ## combine the panels ##
  ########################
  
  plot_list <- list(p_a, p_b, p_c)
  
  p <- patchwork::wrap_plots(plot_list) +
    patchwork::plot_annotation(tag_levels = "A")
  
  ggsave("plots/figure_s1.png", p, width = 12.5, height = 5)
  
}
