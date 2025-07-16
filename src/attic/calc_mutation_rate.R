library(tidyverse)

get_cum_num_mut <- function(pid, dat) {
  
  print(pid)
  
  df <- dat %>%
    filter(patient == pid)
  
  if (nrow(df) > 0) {
    
    times <- sort(unique(df$time_from_diagnosis))
    
    SNPS <- c()
    
    CUM_NUM_MUT <- tibble(
      time_from_diagnosis = 0,
      num_inp = 0,
      num_out = 0
    )
    
    for (time in times) {
      
      snps <- df %>%
        filter(time_from_diagnosis == time) %>%
        pull(snp)
      
      n1 <- length(SNPS)
      SNPS <- unique(c(SNPS, snps))
      n2 <- length(SNPS)
      
      cum_num_mut <- tibble(
        time_from_diagnosis = time,
        num_inp = n1,
        num_out = n2
      )
      
      CUM_NUM_MUT <- rbind(CUM_NUM_MUT, cum_num_mut)
      
    }
    
    CUM_NUM_MUT %>%
      mutate(patient = pid) %>%
      select(4, 1:3)
    
  }
}

metadata <- read_delim("data/metadata.tsv") %>%
  filter(!multiple_carriage)

patients <- metadata %>%
  filter(time_from_diagnosis > 1) %>%
  pull(patient) %>%
  unique()

metadata <- metadata %>%
  filter(patient %in% patients)

variants <- read_delim("data/filtered_variants.tsv")

variants$af <- variants$INFO %>%
  purrr::map(function(x) x %>%
               str_replace(".*AF=", "") %>%
               str_replace(".SB.*", "") %>%
               as.numeric()) %>%
  unlist()

dat <- variants %>%
  filter(af > 0.5) %>%
  mutate(snp = paste0(CHROM, "___", POS)) %>%
  dplyr::rename(isolate = 1) %>%
  select(isolate, snp, af) %>%
  inner_join(metadata, by = "isolate")

inp <- metadata$patient %>%
  purrr::map(function(x) get_cum_num_mut(x, dat)) %>%
  bind_rows()

mod <- lme4::glmer(
  num_out ~ time_from_diagnosis + (1 | patient),
  data = inp,
  family = poisson
)

values <- mod %>%
  summary() %>%
  coef() %>%
  as.data.frame() %>%
  pull(1)

mutation_rate <- exp(values[1] + values[2])
writeLines(as.character(mutation_rate), "data/mutation_rate.txt")

slope <- round(mutation_rate, 2)

label <- sprintf("%s SNPs per year", slope)

p <- ggplot(inp, aes(x = time_from_diagnosis, y = num_out)) +
  geom_point(colour = "lightgrey") +
  geom_smooth(method = "glm", colour = "steelblue") +
  theme_classic(base_size = 12.5) +
  theme(axis.title = element_text(face = "bold")) +
  labs(x = "Time from diagnosis (years)", y = "Number of SNPs") +
  annotate("text", x = 0, y = 87.5, label = label, hjust = 0)

ggsave("output/mutation_rate.png", p, height = 5, width = 5)
