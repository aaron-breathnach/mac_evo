library(tidyverse)

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

list_meds <- function(dat, i) {
  
  cols <- colnames(dat)
  
  values <- t(dat[i,]) %>% as.data.frame() %>% pull(1)
  
  n <- 0
  out <- c()
  for (i in values) {
    n <- n + 1
    if (i == 1) {
      out <- c(out, cols[n])
    }
  }
  
  return(out)
  
}

make_upset_plot <- function(dat) {
  
  dat$meds <- map(1:nrow(dat), function(x) list_meds(dat, x))
  
  dat %>%
    ggplot(aes(x = meds)) +
    geom_bar() +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -1) +
    ggupset::scale_x_upset(n_intersections = 20) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold")) +
    labs(x = "Antibiotics", y = "Intersection size")
  
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
  
  p <- cowplot::plot_grid(p_a, p_b, labels = c("A", "B"), scale = 0.95)
  
  ggsave("plots/figure_1.png", p, width = 11.25, height = 5, bg = "white")
  
}
