library(tidyverse)

make_figure_s2 <- function() {
  
  metadata <- read_delim("data/metadata.tsv")
  
  ## isolates per study
  df1 <- metadata %>%
    filter(bracken_perc > 70) %>%
    group_by(study) %>%
    tally()
  
  ## max time from diagnosis
  df2 <- metadata %>%
    group_by(patient, study) %>%
    filter(time_from_diagnosis == max(time_from_diagnosis)) %>%
    ungroup() %>%
    select(patient, study, time_from_diagnosis) %>%
    rename(n = 3)
  
  ## num samples per patient
  df3 <- metadata %>%
    group_by(patient, study) %>%
    tally() %>%
    ungroup()
  
  ## plots
  
  pal <- c(
    "Present" = "#FF883E",
    "Van Tonder" = "#012169",
    "Wetzstein" = "#FFCC00"
  )
  
  
  p1 <- ggplot(df1, aes(x = study, y = n)) +
    geom_bar(stat = "identity",
             aes(fill = study),
             colour = "black",
             width = 0.75,
             show.legend = FALSE) +
    geom_text(aes(label = paste0("n=", n), vjust = -0.25)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold")) +
    labs(x = "Study", y = "Number of isolates") +
    scale_fill_manual(values = pal)
  
  p2 <- ggplot(df2, aes(x = study, y = n)) +
    geom_jitter(aes(colour = study, fill = study), alpha = 0.5, show.legend = F) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold")) +
    labs(x = "Study", y = "Max time from diagnosis (years) per patient") +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = pal)
  
  p3 <- ggplot(df3, aes(x = study, y = n)) +
    geom_jitter(aes(colour = study, fill = study), alpha = 0.5, show.legend = F) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold")) +
    labs(x = "Study", y = "Number of isolates per patient") +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = pal)
  
  plot_list <- list(p1, p2, p3)
  
  p <- patchwork::wrap_plots(plot_list) +
    patchwork::plot_annotation(tag_levels = "A")
  
  ggsave("plots/figure_s2.png", p, height = 5, width = 10)
  
}
