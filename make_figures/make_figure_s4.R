library(tidyverse)

shift_axis <- function(p, xmin, xmax, y = 0){
  
  g <- ggplotGrob(p)
  
  ax <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  
  p + 
    annotation_custom(
      grid::grobTree(ax, vp = grid::viewport(y = 1, height = sum(ax$height))), 
      ymax = y, ymin = y
    ) +
    annotate("segment", y = 0, yend = 0, x = xmin, xend = xmax, 
             arrow = arrow(length = unit(0.1, "inches"))) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank())
  
}

plot_timelines <- function(pid, timelines) {
  
  p_inp <- timelines %>%
    filter(patient == pid) %>%
    group_by(patient, date) %>%
    summarise(event = paste0(rev(event), collapse = "<br>")) %>%
    ungroup() %>%
    mutate(type = ifelse(grepl("treatment", event), "treatment", "isolate")) %>%
    mutate(shape = ifelse(type == "isolate", "1f9ec", "1f48a"))
  
  n <- nrow(p_inp)
  
  p_inp$y <- rep(c(1, -1), ceiling(n / 2))[1:n]
  
  vjust <- p_inp %>%
    mutate(vjust = case_when(
      y == -1 & !grepl("<br>", event) ~ +2.75,
      y == 1 &  !grepl("<br>", event) ~ -1.75,
      y == -1 & grepl("<br>", event)  ~ +1.75,
      y == 1 &  grepl("<br>", event)  ~ -0.75
    )) %>%
    pull(vjust)
  
  # vjust <- map(p_inp$y, function(x) if (x < 0) {2.5} else {-1.5})
  
  min_year <- lubridate::year(min(p_inp$date) - 365)
  max_year <- lubridate::year(max(p_inp$date) + 365)
  
  years <- max_year - min_year
  
  x_breaks <- c()
  for (i in 0:years) {
    x_breaks <- c(x_breaks, as.Date(min_year + i))
  }
  
  x_min <- lubridate::as_date(sprintf("%s-01-01", min_year))
  x_max <- lubridate::as_date(sprintf("%s-06-01", max_year))
  
  for (i in 1:nrow(p_inp)) {
    if (grepl("First", p_inp[[i, 3]])) {
      date_1 <- p_inp[[i, 2]]
    }
    if (grepl("Mutation", p_inp[[i, 3]])) {
      date_2 <- p_inp[[i, 2]]
    }
  }
  
  y_min <- min(p_inp$y)
  y_max <- max(p_inp$y)
  
  p <- p_inp %>%
    ggplot(aes(date, y)) +
    annotate("rect",
             xmin = date_1, xmax = date_2, 
             ymin = y_min, ymax = y_max,
             alpha = .1,
             fill = "firebrick") +
    ggalt::geom_lollipop(point.size = 0, colour = "grey") +
    emoGG::geom_emoji(emoji = p_inp$shape) +
    ggtext::geom_richtext(aes(x = date, y = y, label = event),
                          hjust = 0.5,
                          vjust = vjust, 
                          size = 4,
                          fill = NA, 
                          label.color = NA,
                          label.padding = grid::unit(rep(0, 4), "pt")) +
    theme(aspect.ratio = 1,
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_text(size = 12.5),
          plot.title = element_text(size = 15, face = "bold")) +
    scale_x_date(date_breaks = "1 year",
                 date_labels = "%Y") +
    scale_y_continuous(limits = c(1.5 * min(p_inp$y), 1.5 * max(p_inp$y))) +
    expand_limits(x = c(x_min, x_max), y = 1.2)
  
  title <- str_replace(pid, "P", "Patient ")
  
  shift_axis(p, x_min, x_max) +
    theme(panel.grid = element_blank(),
          panel.background = element_blank()) +
    ggtitle(title)
  
}

make_figure_s4 <- function() {
  
  patient_metadata <- read_delim("data/patient_info.tsv") %>%
    filter(Treated == "Yes") %>%
    mutate(patient = gsub("P", "P0", patient))
  
  isolate_metadata <- read_delim("data/metadata.tsv") %>%
    mutate(patient = gsub("P", "P0", patient)) %>%
    filter(study == "Present" & !multiple_carriage & patient %in% patient_metadata$patient)
  
  args <- read_delim("data/card_mycobacterial_args.tsv")
  
  dat <- read_delim("data/filtered_variants.tsv") %>%
    filter(GENOME %in% isolate_metadata$isolate)
  
  dat$allele_frequency <- dat$INFO %>%
    str_replace(".SB=.*", "") %>%
    str_replace(".*=", "") %>%
    as.numeric()
  
  dat$gene <- map(dat$INFO, function(x) {
    str_split(x, "\\|") %>%
      unlist() %>%
      nth(4)
  }) %>%
    unlist()
  
  antibiotics <- c("Azithromycin", "Clarithromycin", "Ethambutol", "Rifampicin")
  
  df <- dat %>%
    rename(isolate = 1) %>%
    filter(grepl("missense_variant", INFO)) %>%
    filter(nchar(gene) < 14) %>%
    select(isolate, gene, allele_frequency) %>%
    mutate(gene = str_replace(gene, "_.*", "")) %>%
    inner_join(args, by = "gene") %>%
    inner_join(isolate_metadata, by = "isolate") %>%
    inner_join(patient_metadata, by = "patient") %>%
    separate_longer_delim(antibiotic, ";") %>%
    filter(antibiotic %in% antibiotics)
  
  DF <- df %>%
    select(patient, antibiotic) %>%
    distinct() %>%
    separate_longer_delim(antibiotic, ";") %>%
    separate_longer_delim(antibiotic, " & ") %>%
    distinct() %>%
    filter(antibiotic %in% antibiotics) %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = "antibiotic", values_from = "value", values_fill = 0) %>%
    column_to_rownames("patient")
  
  acquisitions <- df %>%
    select(patient, gene, antibiotic, date_of_collection) %>%
    group_by(patient, gene) %>%
    filter(date_of_collection == min(date_of_collection)) %>%
    ungroup() %>%
    mutate(gene = paste0("*", gene, "*")) %>%
    group_by(patient, date_of_collection) %>%
    summarise(gene = paste0(unique(gene), collapse = " + ")) %>%
    ungroup() %>%
    mutate(event = sprintf(
      "Mutation%s in %s",
      ifelse(str_count(gene, "\\+") > 0, "s", ""),
      gene)) %>%
    select(patient, event, date_of_collection) %>%
    rename(date = 3)
  
  collections <- isolate_metadata %>%
    filter(patient %in% unique(df$patient)) %>%
    select(patient, date_of_collection) %>%
    arrange(date_of_collection) %>%
    group_by(patient) %>%
    mutate(event = sprintf("%s isolate", str_to_sentence(english::ordinal(row_number())))) %>%
    ungroup() %>%
    select(1, 3, 2) %>%
    rename(date = 3)
  
  treatments <- patient_metadata %>%
    select(patient, `Start of treatment`, `End of treatment`) %>%
    filter(patient %in% rownames(DF)) %>%
    pivot_longer(!patient, names_to = "event", values_to = "date") %>%
    mutate(date = lubridate::dmy(date))
  
  timelines <- rbind(acquisitions, treatments, collections) %>%
    arrange(date) %>%
    distinct()
  
  pids <- unique(timelines$patient)
  
  plot_list <- map(pids, function(x) plot_timelines(x, timelines))
  
  p <- cowplot::plot_grid(plotlist = plot_list,
                          nrow = 1,
                          labels = "AUTO", scale = 0.9)
  
  ggsave("plots/figure_s4.png", p, height = 8.75, width = 17.5, bg = "white")
  
}
