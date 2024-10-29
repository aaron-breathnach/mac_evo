compare_timepoints <- function(previous, current) {
  
  cleared <- c()
  persist <- c()
  
  for (p in previous) {
    if (p %in% current) {
      persist <- c(persist, p)
    } else {
      cleared <- c(cleared, p)
    }
  }
  
  new_inf <- c()
  
  for (c in current) {
    if (!c %in% previous) {
      new_inf <- c(new_inf, c)
    }
  }
  
  previous <- paste0(previous, collapse = ",")
  current  <- paste0(current, collapse = ",")
  cleared  <- paste0(cleared, collapse = ",")
  persist  <- paste0(persist, collapse = ",")
  new_inf  <- paste0(new_inf, collapse = ",")
  
  tibble(haplotype = current,
         cleared = cleared,
         persist = persist,
         new_inf = new_inf)
  
}

.check_haplo_status <- function(df, i, threshold) {
  
  time_from_diagnosis <- df[[i, "time_from_diagnosis"]]
  
  current <- df[i, 4:ncol(df)] %>%
    pivot_longer(cols = everything()) %>%
    filter(value > threshold) %>%
    pull(name)
  
  previous <- df[i - 1, 4:ncol(df)] %>%
    pivot_longer(cols = everything()) %>%
    filter(value > threshold) %>%
    pull(name)
  
  out <- compare_timepoints(previous, current) %>%
    mutate(time_from_diagnosis = time_from_diagnosis) %>%
    select(ncol(.), 2:ncol(.) - 1)
  
  if (i == 2) {
    
    row_one <- tibble(time_from_diagnosis = 0,
                      haplotype = paste0(previous, collapse = ","),
                      cleared = "",
                      persist = "",
                      new_inf = "")
    
    out <- bind_rows(row_one, out)
    
  }
  
  return(out)
  
}

check_haplo_status <- function(df, pid, threshold) {
  
  df <- df %>%
    select(isolate, patient, haplotype, time_from_diagnosis, rel_contribution) %>%
    pivot_wider(names_from = "haplotype",
                values_from = "rel_contribution",
                values_fill = 0) %>%
    arrange(time_from_diagnosis) %>%
    filter(patient == pid)
  
  n <- nrow(df)
  
  purrr::map(2:n, function(x) .check_haplo_status(df, x, threshold)) %>%
    bind_rows() %>%
    mutate(patient = pid) %>%
    select(ncol(.), 1:ncol(.) - 1)
  
}

run_check_haplo_status <- function(haplotype_abundance, metadata, threshold = 1) {
  
  df <- inner_join(haplotype_abundance, metadata, by = "isolate")
  
  pids <- unique(df$patient)
  
  purrr::map(pids, function(x) check_haplo_status(df, x, threshold)) %>%
    bind_rows()
  
}

.compare_first_vs_last <- function(x, i) {

  patient <- x[[i, "patient"]]

  compare_timepoints(str_split(x[[i, 2]], ",") %>% unlist(),
                     str_split(x[[i, 3]], ",") %>% unlist) %>%
    mutate(patient = patient) %>%
    select(ncol(.), 1:ncol(.) - 1)

}

compare_first_vs_last <- function(haplotype_abundance, metadata, threshold = 1, long = FALSE) {
  
  df <- run_check_haplo_status(haplotype_abundance, metadata, threshold) %>%
    group_by(patient) %>%
    filter(time_from_diagnosis == max(time_from_diagnosis) | time_from_diagnosis == min(time_from_diagnosis)) %>%
    ungroup() %>%
    mutate(timepoint = ifelse(time_from_diagnosis == 0, "baseline", "endpoint")) %>%
    select(patient, timepoint, haplotype) %>%
    distinct() %>%
    pivot_wider(names_from = "timepoint", values_from = "haplotype")
  
  first_vs_last <- purrr::map(1:nrow(df), function(x) .compare_first_vs_last(df, x)) %>%
    bind_rows() %>%
    mutate(n_cleared = str_count(cleared, "H")) %>%
    mutate(n_persist = str_count(persist, "H")) %>%
    mutate(n_new_inf = str_count(new_inf, "H"))
  
  if (long) {
  
    first_vs_last <- first_vs_last %>%
      select(patient, n_cleared, n_persist, n_new_inf) %>%
      pivot_longer(!patient, names_to = "status", values_to = "n") %>%
      mutate(status = recode(status,
                             "n_cleared" = "Cleared",
                             "n_persist" = "Persisted",
                             "n_new_inf" = "Acquired")) %>%
      mutate(status = factor(status, levels = c("Persisted", "Cleared", "Acquired")))
    
  }
  
  return(first_vs_last)
  
}

first_vs_last_bar <- function(pid, first_vs_last) {
  
  p_inp <- first_vs_last %>%
    filter(patient %in% pid)
  
  ggplot(p_inp, aes(x = status, y = n, label = n)) +
    geom_bar(stat = "identity") +
    geom_text(vjust = -0.25) +
    labs(x = "Status", y = "Number of haplotypes") +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold")) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.1)),
      breaks = function(x) unique(floor(pretty(seq(min(x), (max(x) + 1) * 1.1))))
    )
  
}
