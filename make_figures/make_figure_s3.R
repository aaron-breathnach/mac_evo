library(tidyverse)

make_figure_s3 <- function() {
  
  dat <- read_delim("data/muttui_combined.csv") %>%
    mutate(group = Substitution %>%
             str_replace(".*\\[", "") %>%
             str_replace("\\].*", "")) %>%
    arrange(group) %>%
    mutate(Substitution = factor(Substitution, levels = .$Substitution)) %>%
    mutate(prop = Number_of_mutations / sum(Number_of_mutations)) %>%
    group_by(group) %>%
    summarise(Number_of_mutations = sum(Number_of_mutations)) %>%
    mutate(prop = Number_of_mutations / sum(Number_of_mutations))
  
  p <- ggplot(dat, aes(x = group, y = prop)) +
    geom_bar(aes(fill = group),
             stat = "identity",
             colour = "black",
             show.legend = FALSE) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_viridis_d(option = "turbo") +
    labs(x = "Subtitution", y = "Proportion of mutations")
  
  ggsave("plots/figure_s3.png", p, width = 6, height = 4)
  
}
