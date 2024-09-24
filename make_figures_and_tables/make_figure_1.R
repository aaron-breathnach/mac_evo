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

make_boxplot <- function(dat, value, title) {
  
  p_inp <- dat %>%
    select(study, all_of(value)) %>%
    rename(value = 2)
  
  ggplot(p_inp, aes(x = study, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(colour = study, fill = study),
                alpha = 0.50,
                width = 0.25,
                show.legend = FALSE) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold")) +
    labs(x = "Study", y = title) +
    scale_colour_manual(values = c("#ed4037", "#262161", "#faaf40")) +
    scale_fill_manual(values = c("#ed4037", "#262161", "#faaf40"))
  
}

make_figure_1 <- function() {
  
  patient_metadata <- read_delim("data/patient_info.tsv")
  isolate_metadata <- read_delim("data/metadata.tsv")
  
  p_a <- make_figure_1a(patient_metadata)
  
  p_b <- make_figure_1b(patient_metadata)
  
  ## Figure 1C
  
  dat_c <- tibble(study = c("Present", "Van Tonder", "Wetzstein"),
                  n = c(48, 135, 102))
  
  p_c <- ggplot(dat_c, aes(x = study, y = n)) +
    geom_bar(aes(fill = study),
             stat = "identity",
             show.legend = FALSE,
             colour = "black") +
    geom_text(aes(label = n), vjust = -0.25, size = 5) +
    scale_fill_manual(values = c("#ed4037", "#262161", "#faaf40")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold")) +
    labs(x = "Study", y = "Number of genomes")
  
  ## Figure 1D
  
  dat_d <- isolate_metadata %>%
    filter(time_from_diagnosis >= 1) %>%
    group_by(patient) %>%
    filter(time_from_diagnosis == max(time_from_diagnosis)) %>%
    ungroup()
  
  p_d <- make_boxplot(dat_d,
                      "time_from_diagnosis",
                      "Max time from diagnosis (years)")
  
  ## Figure 1E
  
  dat_e <- isolate_metadata %>%
    group_by(patient, study) %>%
    tally() %>%
    ungroup()
  
  p_e <- make_boxplot(dat_e,
                      "n",
                      "Number of isolates per patient")
  
  ## combine the plots
  
  top <- cowplot::plot_grid(p_a, p_b, labels = c("A", "B"), scale = 0.95)
  
  bottom <- cowplot::plot_grid(p_c, p_d, p_e,
                               nrow = 1,
                               scale = 0.95,
                               labels = c("C", "D", "E"))
  
  p <- cowplot::plot_grid(top, bottom, nrow = 2, rel_heights = c(1, 0.75))
  ggsave("plots/figure_1.png", p, width = 11.25, height = 8.75, bg = "white")
  
}
