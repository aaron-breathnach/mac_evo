library(tidyverse)

source("src/compare_first_vs_last.R")

make_figure_3a <- function(relab, metadata) {
  
  p_inp <- inner_join(relab, metadata, by = "isolate") %>%
    group_by(patient, time_from_diagnosis, haplotype) %>%
    summarise(rel_contribution = mean(rel_contribution)) %>%
    mutate(patient = patient %>%
             str_replace("P", "P0") %>%
             str_replace("VT_", "V") %>%
             str_replace("Wetzstein_", "W"))
  
  p <- ggplot(p_inp, aes(x = time_from_diagnosis, y = rel_contribution)) +
    facet_grid(~ patient, scales = "free_x", space = "free") +
    geom_area(aes(fill = haplotype),
              stat = "identity",
              position = "fill") +
    labs(x = "Time from diagnosis (years)",
         y = "Relative haplotype frequency",
         colour = "Haplotype",
         fill = "Haplotype") +
    scale_x_continuous(expand = expansion(mult = c(0, 0)),
                       breaks = 1:12) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    scale_fill_brewer(palette = "Set3") +
    theme_bw(base_size = 12.5) +
    theme(panel.grid = element_blank(),
          strip.text.x = element_text(face = "bold", angle = 90),
          axis.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.title = element_text(face = "bold"))
  
  return(p)
  
}

make_figure_3b <- function(relab, metadata) {
  
  df <- run_check_haplo_status(relab, metadata, threshold = 1) %>%
    group_by(patient) %>%
    filter(time_from_diagnosis == max(time_from_diagnosis) | time_from_diagnosis == min(time_from_diagnosis)) %>%
    ungroup() %>%
    mutate(timepoint = ifelse(time_from_diagnosis == 0, "baseline", "endpoint")) %>%
    select(patient, timepoint, haplotype) %>%
    distinct() %>%
    pivot_wider(names_from = "timepoint", values_from = "haplotype")
  
  first_vs_last <- compare_first_vs_last(relab, metadata, long = TRUE)
  
  p_inp <- first_vs_last %>%
    mutate(status = factor(status, levels = c("Persisted", "Cleared", "Acquired"))) %>%
    mutate(patient = patient %>%
             str_replace("P", "P0") %>%
             str_replace("VT_", "V") %>%
             str_replace("Wetzstein_", "W"))
  
  ggplot(p_inp, aes(x = patient, y = n, label = n, group = status)) +
    geom_bar(aes(fill = status),
             stat = "identity",
             position = position_dodge2(width = 0.9, preserve = "single", padding = 0),
             colour = "black") +
    geom_text(position = position_dodge2(width = 0.9), vjust = -0.25) +
    labs(x = "Patient", y = "Number of haplotypes", fill = "Status") +
    theme_classic(base_size = 12.5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold")) +
    palettetown::scale_fill_poke(pokemon = 6, spread = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
}

make_figure_3 <- function() {
  
  relab <- readRDS("data/haplotype_deconstructor.RDS")[[2]]
  
  metadata <- read_delim("data/metadata.tsv") %>%
    filter(bracken_pass & multiple_carriage)
  
  p_a <- make_figure_3a(relab, metadata)
  p_b <- make_figure_3b(relab, metadata)
  
  p <- cowplot::plot_grid(p_a, p_b, nrow = 2, labels = "AUTO", scale = 0.95)
  
  ggsave("plots/figure_3.png", p, width = 12.5, height = 8.75, bg = "white")
  
}
