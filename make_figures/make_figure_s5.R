library(tidyverse)

calc_perc_pers <- function(pid, metadata, dat) {
  
  meta <- metadata %>%
    filter(patient == pid & time_from_diagnosis > 0) %>%
    arrange(time_from_diagnosis)
  
  n <- nrow(meta) - 1
  
  out <- c()
  
  for (i in 1:n) {
    
    isolate_1 <- meta[[i, "isolate"]]
    isolate_2 <- meta[[i + 1, "isolate"]]
    
    tmp <- dat %>%
      filter(GENOME %in% c(isolate_1, isolate_2)) %>%
      mutate(snp = sprintf("%s___%s", CHROM, POS)) %>%
      select(snp, GENOME) %>%
      arrange(snp) %>%
      mutate(timepoint = ifelse(GENOME == isolate_1, "T1", "T2")) %>%
      select(-GENOME) %>%
      mutate(value = 1) %>%
      distinct() %>%
      pivot_wider(names_from = "timepoint", values_from = "value", values_fill = 0)
    
    if (nrow(tmp) > 0) {
      
      tmp <- tmp %>%
        mutate(value = ifelse(T1 == T2, 1, 0)) %>%
        summarise(perc = 100 * sum(value) / n()) %>%
        mutate(timepoint = sprintf("T%s to T%s", i, i + 1)) %>%
        mutate(patient = pid) %>%
        select(3, 2, 1)
      
    } else {
      
      tmp <- c()
      
    }
    
    out <- rbind(out, tmp)
    
  }
  
  return(out)
  
}

make_figure_s3 <- function() {
  
  metadata <- read_delim("data/metadata.tsv") %>%
    filter(bracken_pass & !multiple_carriage)
  
  dat <- read_delim("data/filtered_variants.tsv") %>%
    filter(GENOME %in% metadata$isolate & grepl("missense", INFO))
  
  res <- map(patients, function(x) calc_perc_pers(x, metadata, dat)) %>%
    bind_rows()
  
  MEAN <- mean(res$perc)
  
  label <- paste0("Mean: ", round(MEAN, 2), "%")
  
  dens <- density(res$perc)
  
  p <- ggplot(res, aes(x = perc)) +
    geom_density(fill = "steelblue") +
    geom_vline(xintercept = MEAN,
               linetype = "dashed",
               colour = "grey",
               linewidth = 0.75) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    annotate(geom = "text",
             x = 1.25 * MEAN,
             y = max(dens$y),
             label = label,
             hjust = 0) +
    labs(x = "Percentage of persistent nonsynonymous variants", y = "Density")
  
  ggsave("plots/figure_s5.png", p, width = 6, height = 4)
  
}
