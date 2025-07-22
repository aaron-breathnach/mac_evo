library(tidyverse)

## part a
make_figure_s3a <- function(snp_dists) {
  
  within <- snp_dists %>%
    filter(patient.x == patient.y & genome_1 != genome_2) %>%
    mutate(comparison = "Within patient") %>%
    select(comparison, dist)
  
  between <- snp_dists %>%
    filter(patient.x != patient.y) %>%
    select(patient.x, dist) %>%
    mutate(comparison = "Between patient") %>%
    select(comparison, dist)
  
  p_inp <- rbind(within, between) %>%
    mutate(comparison = factor(comparison, levels = c("Within patient", "Between patient")))
  
  ggplot(p_inp, aes(x = log10(dist + 1))) +
    geom_density(aes(colour = comparison, fill = comparison), alpha = 0.5) +
    labs(x = "log(Distance)", y = "Density", colour = "", fill = "") +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold"),
          legend.position.inside = c(2/3, 0.9)) +
    scale_colour_manual(values = c("midnightblue", "orangered")) +
    scale_fill_manual(values = c("midnightblue", "orangered"))
  
}

make_figure_s3b <- function(snp_dists, metadata, mut_rat) {
  
  tmp <- snp_dists %>%
    filter(patient.x == patient.y) %>%
    select(patient.x, dist) %>%
    distinct() %>%
    group_by(patient.x) %>%
    filter(dist == max(dist)) %>%
    ungroup() %>%
    mutate(
      IQR = IQR(dist),
      O_upper = quantile(dist, probs = 0.75, na.rm = FALSE) + 3 * IQR,
      O_lower = quantile(dist, probs = 0.25, na.rm = FALSE) - 3 * IQR
    ) %>%
    filter(O_lower <= dist & dist <= O_upper)
  
  max_dis <- as.numeric(quantile(tmp$dist, probs = 0.95))
  pids <- unique(tmp$patient.x)
  
  p_inp <- snp_dists %>%
    select(genome_1, genome_2, dist) %>%
    inner_join(metadata, by = c("genome_1" = "isolate")) %>%
    inner_join(metadata, by = c("genome_2" = "isolate")) %>%
    filter(genome_1 != genome_2 & patient.x == patient.y) %>%
    filter(patient.x %in% pids) %>%
    filter(time_from_diagnosis.x == 0) %>%
    select(patient.x, time_from_diagnosis.y, dist) %>%
    setNames(gsub("\\..*", "", names(.))) %>%
    arrange(patient, time_from_diagnosis)
  
  label <- sprintf("%s SNPs/year", round(mut_rat, 2))
  y_pos <- 1.1 * max(p_inp$dist)
  
  ggplot(p_inp, aes(x = time_from_diagnosis, y = dist)) +
    geom_point(colour = "lightgrey") +
    geom_smooth(method = "lm", colour = "steelblue") +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold")) +
    labs(x = "Time from diagnosis (years)", y = "Number of SNPs") +
    annotate("text", x = 0, y = y_pos, label = label, hjust = 0)
  
}

make_figure_s3 <- function() {
  
  metadata <- read_delim("data/metadata.tsv") %>%
    filter(!multiple_carriage) 
  
  snp_dists <- read_delim("data/snp_dists.cluster_wide.tsv")
  
  mut_rat <- readLines("data/mutation_rate.txt") %>%
    as.numeric()
  
  p1 <- make_figure_s3a(snp_dists)
  p2 <- make_figure_s3b(snp_dists, metadata, mut_rat)
  p  <- cowplot::plot_grid(p1, p2,
                          nrow = 1,
                          labels = "AUTO",
                          rel_widths = c(1.5, 1))
  
  ggsave("plots/figure_s3.png", p, width = 10, height = 3.75)
  
}
