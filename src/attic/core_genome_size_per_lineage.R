files <- list.files("ska_x_fastbaps",
                    pattern = "gubbins.core.aln",
                    recursive = TRUE,
                    full.names = TRUE)

get_len <- function(x) {
  
  fasta <- ape::read.FASTA(x)
  
  cluster <- x %>%
    str_split("/") %>%
    unlist() %>%
    nth(2) %>%
    str_replace(".*_", "")
  
  tibble(
    cluster = cluster,
    core_size = length(fasta[[1]])
  )
  
}

fastbaps <- read_delim("data/fastbaps.tsv") %>%
  select(1, 2) %>%
  rename(cluster = 2) %>%
  mutate(cluster = str_pad(cluster, 2, "left", "0"))

isolates <- read_delim("data/international_clusters.tsv") %>%
  pull(isolate)

p_inp <- purrr::map(files, function(x) get_len(x)) %>%
  bind_rows() %>%
  inner_join(fastbaps, by = "cluster") %>%
  filter(isolate %in% isolates) %>%
  select(cluster, core_size) %>%
  distinct() %>%
  arrange(core_size) %>%
  mutate(cluster = factor(cluster, levels = .$cluster))

p <- ggplot(p_inp, aes(x = core_size, y = cluster)) +
  geom_bar(stat = "identity", colour = "black", fill = "steelblue") +
  geom_text(aes(label = sprintf("n=%s", scales::comma(core_size))),
            hjust = -0.05) +
  scale_x_log10(expand = expansion(mult = c(0, 1/3)),
                labels = scales::comma) +
  labs(x = "Num. core SNPs (post-Gubbins)", y = "Lineage") +
  theme_classic(base_size = 12.5) +
  theme(axis.title = element_text(face = "bold"))

ggsave("output/core_size_per_cluster.png", p, height = 3.75, width = 5)
