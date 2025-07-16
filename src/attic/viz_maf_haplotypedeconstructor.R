library(tidyverse)

haplo <- readRDS("data/haplotype_deconstructor.RDS")

mat <- haplo$decomposed$fitted

mafs <- tibble(
  snp = rownames(mat),
  maf = rowSums(ifelse(mat > 0.05, 1, 0)) / ncol(mat)
)

p <- ggplot(mafs, aes(x = maf)) +
  geom_histogram(fill = "steelblue") +
  theme_classic(base_size = 12.5) +
  theme(axis.title = element_text(face = "bold")) +
  labs(x = "Minor allele frequency", y = "Count")

ggsave("output/minor_allele_frequencies.png",
       p,
       width = 5,
       height = 3.75)
