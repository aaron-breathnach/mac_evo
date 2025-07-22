library(tidyverse)

make_figure_s7 <- function() {
  
  p_inp_1 <- list.files(
    "muttui/output",
    pattern = "mutational_spectrum.csv",
    full.names = TRUE
  ) %>%
    purrr::map(function(x) {
      
      isolate <- x %>%
        basename() %>%
        str_replace("_mutational_spectrum.csv", "")
      
      read_delim(x) %>%
        mutate(isolate = isolate) %>%
        mutate(group = str_sub(Substitution, 3, 5)) %>%
        select(3, 4, 1, 2)
      
    }) %>%
    bind_rows() %>%
    group_by(isolate, group) %>%
    summarise(num_mut = sum(Number_of_mutations)) %>%
    ungroup() %>%
    group_by(isolate) %>%
    mutate(per_mut = 100 * num_mut / sum(num_mut)) %>%
    ungroup()
  
  p1 <- ggplot(p_inp_1, aes(x = group, y = per_mut)) +
    geom_jitter(aes(colour = group, fill = group), show.legend = FALSE, alpha = 0.5) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold")) +
    labs(x = "Mutation", y = "Percentage of mutations per sample") +
    scale_colour_viridis_d(option = "turbo") +
    scale_fill_viridis_d(option = "turbo")
  
  p_inp_2 <- read_delim("muttui/combined_mutational_spectra.csv") %>%
    mutate(group = str_sub(Substitution, 3, 5))
  
  p2 <- ggplot(p_inp_2, aes(x = Substitution, y = Number_of_mutations)) +
    facet_grid(~group, scales = "free_x", space = "free_x") +
    geom_bar(stat = "identity", aes(fill = group), show.legend = FALSE) +
    theme_bw(base_size = 12.5) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(face = "bold"),
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold")) +
    scale_x_discrete(expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                       labels = scales::comma) +
    labs(x = "Mutation", y = "Combined number of mutations") +
    scale_fill_viridis_d(option = "turbo")
  
  p <- cowplot::plot_grid(p1, p2, nrow = 2, labels = "AUTO", align = "v", axis = "lr")
  
  ggsave("plots/figure_s7.png", p, width = 10, height = 7.5)
  
}