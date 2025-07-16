library(tidyverse)

add_ann_ids <- function(variants) {
  
  variants$annotation_id <- variants$INFO %>%
    purrr::map(function(x) x %>%
                 str_split("\\|") %>%
                 unlist() %>%
                 nth(5)) %>%
    unlist()
  
  return(variants)
  
}

add_maf <- function(variants) {
  
  variants$maf <- variants$INFO %>%
    purrr::map(function(x) x %>%
                 str_replace(".*AF=", "") %>%
                 str_replace(";.*", "") %>%
                 as.numeric()) %>%
    unlist()
  
  return(variants)
  
}

viz_maf <- function(pid_x_snp, metadata, dat) {
  
  x <- str_split(pid_x_snp, "___") %>%
    purrr::map(function(x) {
      string <- unlist(x)
      c(string[1], paste0(string[2], "___", string[3]))
    }) %>%
    unlist()
  
  pid <- x[1]
  mut <- x[2]
  
  isolates <- metadata %>%
    filter(patient == pid) %>%
    pull(isolate)
  
  gene <- dat %>%
    filter(snp == mut) %>%
    pull(gene_name) %>%
    unique()
  
  output <- c()
  
  for (iso in isolates) {
    
    time_from_diagnosis <- metadata %>%
      filter(isolate == iso) %>%
      pull(time_from_diagnosis)
    
    if (iso %in% dat$isolate) {
      
      maf <- dat %>%
        filter(isolate == iso & snp == mut) %>%
        pull(maf)
      
    } else {
      
      maf <- 0
      
    }
    
    tmp <- tibble(patient = pid,
                  isolate = iso,
                  time_from_diagnosis = time_from_diagnosis,
                  snp = mut,
                  gene = gene,
                  maf = maf)
    
    output <- rbind(output, tmp)
    
  }
  
  p_inp <- output %>%
    mutate(colour = ifelse(maf > 0, "Detected", "Not detected"))
  
  tidy_gene <- paste0("*", gsub("_.*", "", gene), "*")
  
  tidy_pid <- pid %>%
    str_replace("P", "P0") %>%
    str_replace("VT_", "V") %>%
    str_replace("Wetzstein_", "W")
  
  title <- sprintf("**%s in %s**", tidy_gene, tidy_pid)
  
  ggplot(p_inp, aes(x = time_from_diagnosis, y = maf)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
    geom_line() +
    geom_point(aes(fill = colour), pch = 21, show.legend = FALSE) +
    theme_classic(base_size = 12.5) +
    theme(axis.title = element_text(face = "bold"),
          plot.title = ggtext::element_markdown()) +
    scale_fill_manual(values = c("steelblue", "grey")) +
    labs(title = title,
         x = "Time from diagnosis (years)",
         y = "Minor allele frequency",
         colour = "")
  
}

## output from `run_mutational_distribution_analysis.R`
res <- read_delim("data/gen_pos.tsv") %>%
  filter(bonferroni <= 0.05) %>%
  select(gene_id, gene_id_prokka, annotation_id, size) %>%
  rename(gene_name = 1)

annotation_ids <- res$annotation_id
## metadata
tmp <- read_delim("data/metadata.tsv") %>%
  filter(!multiple_carriage & time_from_diagnosis > 0)
## patients with >3 isolates
patients <- tmp %>%
  group_by(patient) %>%
  tally() %>%
  ungroup() %>%
  filter(n > 1) %>%
  pull(patient)

metadata <- tmp %>%
  filter(patient %in% patients) %>%
  arrange(time_from_diagnosis) %>%
  group_by(patient) %>%
  mutate(t = row_number()) %>%
  ungroup() %>%
  mutate(timepoint = sprintf("T%s-T%s",
                             str_pad(t - 1, 2, "left", "0"),
                             str_pad(t, 2, "left", "0"))) %>%
  select(-t)

## panaroo gene presence/absence matrix
presence_absence <- read_delim("data/panaroo/gene_presence_absence.csv")

non_unique_gene_names <- presence_absence %>%
  pivot_longer(
    cols = !c(Gene, `Non-unique Gene name`, Annotation),
    names_to = "isolate", 
    values_to = "annotation_id"
  ) %>%
  filter(annotation_id %in% annotation_ids) %>%
  pull(`Non-unique Gene name`)

gene2anno <- presence_absence %>%
  filter(`Non-unique Gene name` %in% non_unique_gene_names) %>%
  pivot_longer(
    cols = !c(Gene, `Non-unique Gene name`, Annotation),
    names_to = "isolate", 
    values_to = "annotation_id"
  )

## variants
dat <- read_delim("data/filtered_variants.tsv") %>%
  filter(grepl("missense_variant", INFO) & GENOME %in% metadata$isolate) %>%
  add_ann_ids() %>%
  add_maf() %>%
  inner_join(gene2anno, by = "annotation_id") %>%
  select(1:3, Gene, maf) %>%
  rename(gene_name = 4) %>%
  dplyr::rename(isolate = 1) %>%
  inner_join(metadata, by = "isolate") %>%
  inner_join(res, by = "gene_name") %>%
  arrange(timepoint) %>%
  mutate(snp = sprintf("%s___%s", CHROM, POS)) %>%
  mutate(pid_x_snp = sprintf("%s___%s", patient, snp)) %>%
  select(-CHROM, -POS)

write_tsv(dat, "data/sig_gen_maf.tsv")

pids_x_snps <- dat %>%
  group_by(gene_id_prokka, pid_x_snp) %>%
  tally() %>%
  ungroup() %>%
  filter(n > 1)

## number of variants that persisted
numb_pers <- pids_x_snps %>%
  group_by(gene_id_prokka) %>%
  summarise(numb_pers = n())

missing <- res$gene_id_prokka[!res$gene_id_prokka %in% numb_pers$gene_id_prokka] %>%
  purrr::map(function(x) tibble(gene_id_prokka = x, numb_pers = 0)) %>%
  bind_rows()

numb_pers <- rbind(numb_pers, missing)

## visualise the percentage of variants that persisted
p_inp <- dat %>%
  select(gene_id_prokka, patient, size) %>%
  distinct() %>%
  group_by(gene_id_prokka, size) %>%
  summarise(total = n()) %>%
  ungroup() %>%
  inner_join(numb_pers, by = "gene_id_prokka") %>%
  mutate(perc_pers = 100 * numb_pers / total) %>%
  arrange(perc_pers) %>%
  mutate(gene_id_prokka = factor(gene_id_prokka, levels = .$gene_id_prokka))

p1 <- ggplot(p_inp, aes(x = perc_pers, y = gene_id_prokka)) +
  geom_bar(stat = "identity",
           fill = "steelblue",
           colour = "black",
           width = 1) +
  geom_text(aes(label = paste0(round(perc_pers), "%")),
            hjust = -0.125) +
  theme_classic(base_size = 12.5) +
  theme(axis.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "italic")) +
  labs(x = "Percent persistence", y = "Gene") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0)))

ggsave("output/persistence.png", p1, width = 5, height = 3.75)

## visualise the allele frequencies of persistent variants
pids <- gsub("___.*", "", pids_x_snps)

maf <- dat %>%
  filter(pid_x_snp %in% pids_x_snps$pid_x_snp) %>%
  select(patient, isolate, time_from_diagnosis, snp, maf)

plot_list <- pids_x_snps$pid_x_snp %>%
  purrr::map(function(x) viz_maf(x, metadata, dat))

p2 <- cowplot::plot_grid(plotlist = plot_list, nrow = 2)

ggsave("output/persistence_maf.png", p2, width = 17.5, height = 7.5)
