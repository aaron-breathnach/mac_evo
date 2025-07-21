library(tidyverse)

get_mut_rat_and_transm_thresh <- function() {
  
  metadata <- read_delim("data/metadata.tsv")
  
  snp_dists <- read_delim("data/snp_dists.cluster_wide.tsv")
  
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
  
  dat <- snp_dists %>%
    select(genome_1, genome_2, dist) %>%
    inner_join(metadata, by = c("genome_1" = "isolate")) %>%
    inner_join(metadata, by = c("genome_2" = "isolate")) %>%
    filter(genome_1 != genome_2 & patient.x == patient.y) %>%
    filter(patient.x %in% pids) %>%
    filter(time_from_diagnosis.x == 0) %>%
    select(patient.x, time_from_diagnosis.y, dist) %>%
    setNames(gsub("\\..*", "", names(.))) %>%
    arrange(patient, time_from_diagnosis)
  
  formula <- dist ~ time_from_diagnosis + (1 | patient)
  
  mod <- lme4::lmer(formula, data = dat)
  
  mut_rat <- broom.mixed::tidy(
    mod,
    conf.int = TRUE,
    exponentiate = TRUE,
    effects = "fixed"
  ) %>%
    pull(3) %>%
    nth(2)
  
  transm_thresh <- floor(max_dis + mut_rat)
  
  writeLines(as.character(mut_rat), "data/mutation_rate.txt")
  
  writeLines(as.character(transm_thresh), "data/threshold.txt")
  
}

if (sys.nframe() == 0) {
  get_mut_rat_and_transm_thresh()
}