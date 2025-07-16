library(tidyverse)

dat <- read_delim("data/filtered_variants_summary.tsv") %>%
  pivot_longer(!genome) %>%
  group_by(genome) %>%
  mutate(perc = 100 * value / sum(value)) %>%
  ungroup() %>%
  mutate(name = recode(
    name,
    "num_high_qual" = "Passed",
    "num_prob_area" = "Failed: too dense",
    "num_prob_dist" = "Failed: too close"
  ) %>%
    factor(levels = c("Passed", "Failed: too dense", "Failed: too close")
    )
  )

p <- ggplot(dat, aes(x = name, y = perc)) +
  geom_violin(aes(fill = name), show.legend = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(x = "Outcome", y = "Percentage of variants (%)") +
  theme_classic(base_size = 12.5) +
  theme(axis.title = element_text(face = "bold")) +
  scale_fill_manual(values = c("#fdf289", "#46edc8", "#374d7c"))

filename <- "output/snp_qc_figure.png"

ggsave(filename, p, width = 5, height = 4)
