This README outlines the steps involved in the analysis described by Walsh *et al.*. For queries, please contact [Aaron Walsh](mailto:awalsh12@tcd.ie?subject=within-host%20evolution%20of%20MAC).

## Step 01: download the reference genome, Kraken2 database, and CARD database

```
python src/download_databases.py
```

## Step 02: annotate the reference genome using Prokka

```
docker run -v $PWD:/data staphb/prokka:latest \
prokka \
databases/GCF_009741445.1/GCF_009741445.1_ASM974144v1_genomic.fna \
--outdir prokka/reference \
--prefix reference \
--genus Mycobacterium \
--usegenus \
--force \
--centre X \
--compliant
```

## Step 03: run QC

```
python src/run_qc.py
sh run_qc.sh
```

## Step 04: run *de novo* assembly using SPAdes and annotate the assemblies using Prokka

```
python src/run_assembly_and_annotation.py
sh run_spades.sh
sh run_prokka.sh
```

## Step 05: extract genome lengths from Prokka output

```
Rscript src/get_genome_lengths.R
```

## Step 06: run Snippy

```
python src/run_snippy.py
sh run_snippy.py
```

## Step 07: run whole genome alignment using snippy-core

```
python src/run_snippy_core.py
sh run_snippy_core.sh
```

## Step 08: calculate pairwise SNP distance using snp-dists

```
docker run -v $PWD:/data staphb/snp-dists snp-dists snippy/core.aln > data/snp_dists.tsv
```

## Step 09: build a phylogentic tree using Gubbins

```
python src/run_run_gubbins.py
cp gubbins/mavium.final_tree.tre data/
```

## Step 10: run fastbaps

```
Rscript src/run_fastbaps.R --fasta snippy/gubbins.aln --out_dir data
```

## Step 11: list patients containing multiple strains

```
Rscript src/list_reinfected_patients.R \
--snp_dists data/snp_dists.tsv \
--metadata data/metadata.tsv \
--threshold 20
```

## Step 12: run within-host LoFreq

```
python src/run_lofreq_common_reference.py
```

## Step 13: run SnpEff

```
python src/run_snpeff.py data/metadata.tsv .lofreq.vcf.gz
cp within_host_evolution/*/*.lofreq.snpeff.vcf data/lofreq/within_host/
```

## Step 14: filter vcfs using the sliding window method described by [Tonkin-Hill et al. (2022)](https://doi.org/10.1038/s41564-022-01238-1)

```
Rscript src/filter_vcf.R \
--patients data/patients.txt \
--metadata data/metadata.tsv \
--gen_len data/genome_lengths.tsv \
--out_dir data
```

## Step 15: make a dndscv formatted RefCDS object using the reference GFF

```
Rscript src/make_dndscv_fmt_ref_tab.R --gff prokka/reference/reference.gff --out_dir data
```

## Step 16: list mycobacterial antibiotic resistance genes in the CARD database

```
Rscript src/list_card_mycobacterial_args.R --card data/card --out_dir data
```

## Step 17: run Panaroo on the genomes of the first isolate from each patient

```
mkdir -p panaroo/input/
cp prokka/*/*.gff panaroo/input/
python src/run_panaroo.py
sh run_panaroo.sh
mkdir data/panaroo
cp panaroo/gene_presence_absence.csv panaroo/gene_data.csv data/panaroo/
```

## Step 18: run the mutational distribution analysis to identify genes under positive selection

```
run_mutational_distribution_analysis.R \
--panaroo data/panaroo \
--metadata data/metadata.tsv \
--vcf data/filtered_variants.tsv \
--ref_cds data/ref_cds.tsv \
--myco_args data/card_mycobacterial_args.tsv \
--out_dir data
```

## Step 19: parse the Panaroo output to make a list of single-copy core genes

```
python src/parse_panaroo_output.py
```

## Step 20: run within-host LoFreq

```
python src/run_lofreq_common_reference.py
```

##  Step 21: run HaplotypeDeconstructor

```
Rscript src/run_haplotype_deconstructor.R
```

## Step 22: run Mash to identify *M. avium* subspecies

```
python src/run_mash.py
```

## Step 23: run MutTui

```
python src/run_muttui.py
sh run_muttui.sh
```
