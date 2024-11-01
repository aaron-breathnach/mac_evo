library(tidyverse)

source("src/compare_first_vs_last.R")

make_figure_s2 <- function() {
  
  relab <- readRDS("data/haplotype_deconstructor.RDS")[[2]]
  metadata <- read_delim("data/metadata.tsv")
  
  tmp <-  compare_first_vs_last(relab, metadata, threshold = 0, long = FALSE)
  
  dat <- tmp %>%
    select(patient, persist) %>%
    filter(nchar(persist) > 0) %>%
    separate_longer_delim(persist, ",") %>%
    distinct() %>%
    rename(haplotype = 2)
  
  patients_a <- tmp %>%
    filter(n_persist == 0) %>%
    pull(patient)
  
  patients_b <- inner_join(relab, metadata, by = "isolate") %>%
    filter(patient != "VT_198") %>%
    group_by(patient) %>%
    filter(time_from_diagnosis == max(time_from_diagnosis)) %>%
    ungroup() %>%
    inner_join(dat, by = c("patient", "haplotype")) %>%
    group_by(patient) %>%
    summarise(rel_contribution = sum(rel_contribution)) %>%
    filter(rel_contribution < 50) %>%
    pull(patient)
  
  patients <- c(patients_a, patients_b)
  
  meta <- metadata %>%
    rename(time = time_from_diagnosis) %>%
    select(patient, isolate, time) %>%
    group_by(patient) %>%
    filter(time == min(time) | time == max(time)) %>%
    ungroup() %>%
    mutate(timepoint = ifelse(time == min(time), "First", "Last")) %>%
    filter(patient %in% patients) %>%
    select(patient, isolate, timepoint)
  
  df <- inner_join(meta, relab, by = "isolate") %>%
    filter(patient != "VT_198") %>%
    mutate(patient = patient %>%
             str_replace("P", "P0") %>%
             str_replace("VT_", "V") %>%
             str_replace("Wetzstein_", "W"))
  
  p <- ggplot(df, aes(x = timepoint, y = rel_contribution)) +
    facet_wrap(~ patient) +
    geom_bar(aes(fill = haplotype),
             stat = "identity",
             position = "fill",
             colour = "black") +
    scale_fill_brewer(palette = "Set3") +
    theme_bw(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold")) +
    labs(x = "Time-point",
         y = "Relative haplotype frequency",
         fill = "Haplotype")
  
  ggsave("plots/figure_s2.png", p, width = 7.5, height = 5)
  
}
