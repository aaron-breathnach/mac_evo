## Step 1: QC

```
python src/run_qc.py
```

## Step 2: assembly and annotation

```
python src/run_assembly_and_annotation.py
```

## Step 3: run Snippy

```
python src/run_snippy.py
```

## Step 4: run snippy-core

```
python src/run_snippy_core.py
```

## Step 5: run Gubbins

```
python src/build_tree.py
```

## Step 6: run fastbaps

```
Rscript src/run_fastbaps.R --fasta snippy/gubbins.aln --out_dir data
```

## Step 7: list reinfected patients
```
Rscript src/list_reinfected_patients.R --snp_dists data/snp_dists.tsv --metadata data/metadata.tsv --threshold 20
```

## Step 8:
`docker run -v $PWD:/data staphb/snp-dists snp-dists -j 16 snippy/core.aln > snp_dists.tsv`