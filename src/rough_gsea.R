library(tidyverse)

get_go_terms <- function(x) {
  
  x %>%
    str_split(", ") %>%
    unlist() %>%
    purrr::map(function(x) if (grepl("GO:", x)) x) %>%
    unlist() %>%
    paste0(collapse = ",")
  
}

get_cog <- function(x) {
  
  x %>%
    str_split(", ") %>%
    unlist() %>%
    purrr::map(function(x) if (grepl("COG:", x) & nchar(x) == 11) str_replace(x, "COG.", "")) %>%
    unlist() %>%
    paste0(collapse = ",")
  
}

ec2go <- readLines("https://current.geneontology.org/ontology/external2go/ec2go")[-c(1:2)] %>%
  tibble() %>%
  dplyr::rename(x = 1) %>%
  mutate(ec = x %>%
           str_replace(" >.*", "") %>%
           str_replace("EC.", "")) %>%
  mutate(go_accession = x %>%
           str_replace(".* ", "")) %>%
  mutate(go_annotation = x %>%
           str_replace(".* > ", "") %>%
           str_replace(" \\; .*", "")) %>%
  select(2:4)

ref_gff <- "~/Desktop/pilot_study/mac_evo/data/gffs/reference.gff"

prokka <- rtracklayer::readGFF(ref_gff) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(gene = gsub("_.*", "", gene)) %>%
  filter(type == "CDS" & !is.na(eC_number)) %>%
  select(gene, eC_number)

gene2go <- inner_join(prokka, ec2go, by = c("eC_number" = "ec"), relationship = "many-to-many")

parent_nodes <- GOfuncR::get_parent_nodes(gene2go$go_accession) %>%
  as_tibble()

cols <- c("gene", "p", "dn_ds")

raw <- read_delim("data/gen_pos.tsv")

uniq_gene <- raw %>%
  group_by(gene_id_prokka) %>%
  tally() %>%
  ungroup() %>%
  filter(n == 1) %>%
  pull(gene_id_prokka)

gene_stat <- raw %>%
  filter(gene_id_prokka %in% uniq_gene) %>%
  select(gene_id_prokka, p) %>%
  arrange(p) %>%
  rename(gene = 1)

vec <- gene_stat$p
names(vec) <- gene_stat$gene
vec <- sort(vec, decreasing = TRUE)

res <- clusterProfiler::GSEA(vec, TERM2GENE = go_map, pvalueCutoff = 1)
res@result


go_map <- gene2go %>%
  drop_na() %>%
  select(3, 1) %>%
  dplyr::rename(GO = 1) %>%
  distinct() %>%
  as.data.frame() %>%
  clusterProfiler::buildGOmap()

clusterProfiler::enricher(genes, TERM2GENE = go_map, pvalueCutoff = 1, qvalueCutoff = 1)
