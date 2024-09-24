metadata <- read_delim("data/metadata.tsv")

dat <- readRDS("data/haplotype_deconstructor.RDS")[[2]]

inp <- dat %>%
  select(1, 2, 4) %>%
  pivot_wider(names_from = "haplotype", values_from = "rel_contribution") %>%
  column_to_rownames("isolate")

vegdist <- vegan::vegdist(inp) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("genome_1") %>%
  pivot_longer(!genome_1, names_to = "genome_2", values_to = "dist") %>%
  inner_join(metadata, by = c("genome_1" = "isolate")) %>%
  inner_join(metadata, by = c("genome_2" = "isolate")) %>%
  filter(genome_1 != genome_2 & patient.x == patient.y)

p_inp <- vegdist %>%
  select(genome_1, genome_2, patient.x, multiple_carriage.x, dist) %>%
  rename(patient = 3, multiple_carriage = 4) %>%
  group_by(patient) %>%
  filter(dist == max(dist)) %>%
  ungroup() %>%
  select(patient, multiple_carriage, dist) %>%
  distinct() %>%
  mutate(group = ifelse(multiple_carriage, "Yes", "No") %>%
           factor(levels = c("Yes", "No")))

p <- ggplot(p_inp, aes(x = group, y = dist)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    aes(colour = multiple_carriage, fill = multiple_carriage),
    width = 0.25,
    show.legend = FALSE,
    alpha = 0.5
  ) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold")) +
  labs(x = "Multiple carriage", y = "Maximum Bray-Curtis distance per patient")

ggsave("plots/max_dist.png", p, width = 5, height = 5)
