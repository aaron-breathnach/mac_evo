library(tidyverse)

parse_info <- function(x, n) {
  
  x %>%
    str_split("\\|") %>%
    unlist() %>%
    nth(n)
  
}

make_table_s2 <- function() {
  
  vir_gen <- read_delim("data/abricate_vfdb.tsv") %>%
    pull(GENE)
  
  amr_gen <- read_delim("data/card_mycobacterial_args.tsv")  %>%
    pull(gene)
  
  genes <- c(amr_gen, vir_gen)
  
  metadata <- read_delim("data/metadata.tsv") %>%
    filter(study == "Present")
  
  dat <- read_delim("data/minor_vafs/filtered_vcfs.txt") %>%
    filter(GENOME %in% metadata$isolate)
  
  dat$allele_frequency <- map(dat$INFO, function(x) {
    x %>%
      str_replace(".*;AF=", "") %>%
      str_replace(";.*", "") %>%
      as.numeric()
  }) %>%
    unlist()
  
  dat$gene <- map(dat$INFO, function(x) parse_info(x, 4)) %>%
    unlist()
  
  dat$variant <- map(dat$INFO, function(x) parse_info(x, 2)) %>%
    unlist()
  
  dat$position <- map(dat$INFO, function(x) parse_info(x, 12)) %>%
    unlist()
  
  dat$substitution_nucl <- map(dat$INFO, function(x) parse_info(x, 10)) %>%
    unlist()
  
  dat$substitution_prot <- map(dat$INFO, function(x) parse_info(x, 11)) %>%
    unlist()
  
  dat <- dat %>%
    mutate(category = case_when(
      str_replace(gene, "_.*", "") %in% amr_gen ~ "AMR",
      str_replace(gene, "_.*", "") %in% vir_gen ~ "Virulence"
    ))
  
  meta <- metadata %>%
    select(1:4)
  
  gene_products <- metadata %>%
    group_by(patient) %>%
    filter(time_from_diagnosis == min(time_from_diagnosis)) %>%
    pull(isolate) %>%
    map(function(x) {
      gff <- sprintf("%s/%s.gff", "data/gffs", x)
      rtracklayer::readGFF(gff) %>%
        as.data.frame() %>%
        select(gene, product) %>%
        drop_na() %>%
        as_tibble()
    }) %>%
    bind_rows() %>%
    group_by(gene) %>%
    sample_n(1) %>%
    ungroup() %>%
    distinct()
  
  dat <- dat %>%
    filter(nchar(gene) < 14 & variant == "missense_variant") %>%
    left_join(gene_products, by = "gene", relationship = "many-to-many") %>%
    select(1, gene, category, product, position, substitution_nucl, substitution_prot, allele_frequency) %>%
    rename(isolate = 1) %>%
    inner_join(meta, ., by = "isolate")
  
  wb <- openxlsx::createWorkbook()
  
  for (pid in unique(metadata$patient)) {
    
    x <- dat %>%
      filter(patient == pid) %>%
      arrange(desc(allele_frequency))
    
    openxlsx::addWorksheet(wb, sheetName = pid)
    openxlsx::writeDataTable(wb, sheet = pid, x)
    
  }
  
  openxlsx::saveWorkbook(wb, "tables/table_s2.xlsx", overwrite = TRUE)
  
}
  