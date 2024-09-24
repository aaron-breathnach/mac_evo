get_vcf <- function(vcf) {
  
  vcfR::read.vcfR(vcf, verbose = FALSE) %>%
    slot("fix") %>%
    as_tibble() %>%
    mutate(QUAL = as.numeric(QUAL)) %>%
    filter(!grepl("INDEL", INFO) & QUAL >= 30) %>%
    dplyr::rename(`#CHROM` = 1) %>%
    mutate(POS = as.numeric(POS))
  
}