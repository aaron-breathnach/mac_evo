This README outlines the steps involved in the analysis described by Walsh *et al.*. For queries, please contact [Aaron Walsh](mailto:awalsh12@tcd.ie?subject=within-host%20evolution%20of%20MAC).

## 01
`python src/download_databases.py`

## 02 (lists mycobacterial genes in CARD)
`Rscript src/list_mycobacterial_args.R`

## 03: run QC
`python src/run_qc.py`

## 04
`python src/run_assembly_and_annotation.py`

## 05
`python src/run_panaroo.py`

## 06 (extracts a list of single-copy core genes from the Panaroo outputs)
`python src/parse_panaroo_output.py`

## 07
`python src/run_lofreq_within_host.py`

## 08
`python src/run_snpeff.py`

## 09 (get genome lengths)
`Rscript src/get_genome_lengths.R`

## 10 (filter LoFreq output)
`Rscript src/filer_variants.R`

## 11 (makes a dNDdScv formatted RefCDS file)
`Rscript src/make_dndscv_fmt_ref_tab.R`

## 12
`Rscript src/run_mutational_distribution_analysis.R`

## 13 (determines if variants of significant genes persisted over time)
`Rscript src/persistence_of_sig_gen.R`

## 14 (runs muttui to determine mutational spectrum)
`python src/run_muttui.py`

## 15 queries significant genes against STRING
`Rscript src/run_string_analysis.R`

## 16 (generate a species-wide genome alignment using ska) 
`python src/run_ska.py`

## 17
`python src/run_iqtree.py`

## 18
`Rscript src/run_fastbaps.R`

## 19 (runs ska within each fastbaps cluster)
`python src/run_ska_x_fastbaps.py`

## 20 (combines cluster-wide snp matrices)
`Rscript src/parse_dists.R`

## 21 (calculates mutation rate and transmission threshold)

`Rscript src/get_mut_rat_and_transm_thresh.R`
## 22: 

`python src/run_lofreq_common_reference.py`
## 23
`Rscript src/run_haplotype_deconstructor.R`
