library(tidyverse)

dat <- read_delim("data/filtered_variants.tsv")

metadata <- read_delim("data/metadata.tsv")

dat <- read_delim("data/filtered_variants.tsv") %>%
  mutate(AF = INFO %>%
           str_replace(".*AF=", "") %>%
           str_replace("\\;.*", "") %>%
           as.numeric()) %>%
  mutate(DP = INFO %>%
           str_replace(".*DP=", "") %>%
           str_replace("\\;.*", "") %>%
           as.numeric()) %>%
  mutate(DP4 = INFO %>%
           str_replace(".*DP4=", "") %>%
           str_replace("\\;.*", "")) %>%
  separate_wider_delim(DP4, ",", names = c("ref_f", "ref_r", "alt_f", "alt_r")) %>%
  mutate(alt_num = as.numeric(alt_f) + as.numeric(alt_r))

dat %>%
  filter(AF < 0.01)

###