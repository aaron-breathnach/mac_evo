import configparser

def get_dictionary():
    #####################
    ## quality control ##
    #####################
    HOSTILE = '''
    hostile 
    clean 
    --fastq1 reads/raw/{prefix}_R1_001.fastq.gz 
    --fastq2 reads/raw/{prefix}_R2_001.fastq.gz 
    --out-dir reads/hostile 
    --threads {threads}
    '''.replace('\n', '').replace('    ', '')
    FASTP = '''
    docker run -v $PWD:/data staphb/fastp 
    fastp 
    -i reads/hostile/{prefix}_R1_001.clean_1.fastq.gz 
    -I reads/hostile/{prefix}_R2_001.clean_2.fastq.gz 
    -o reads/processed/{prefix}_R1_001.qcd.fastq.gz 
    -O reads/processed/{prefix}_R2_001.qcd.fastq.gz 
    -j reads/fastp/{prefix}.json 
    -h reads/fastp/{prefix}.html
    '''.replace('\n', '').replace('    ', '')
    KRAKEN2 = '''
    docker run -v $PWD:/data staphb/kraken2 
    kraken2 
    --db databases/kraken/k2_standard_08gb/ 
    --threads {threads} 
    --output kraken/{prefix}.output 
    --report kraken/{prefix}.report 
    --paired reads/processed/{prefix}_R1_001.qcd.fastq.gz reads/processed/{prefix}_R2_001.qcd.fastq.gz 
    --gzip-compressed
    '''.replace('\n', '').replace('    ', '')
    BRACKEN = '''
    docker run -v $PWD:/data staphb/bracken 
    bracken 
    -d databases/kraken/k2_standard_08gb/  
    -i kraken/{prefix}.report 
    -o kraken/{prefix}.bracken 
    -l G
    '''.replace('\n', '').replace('    ', '')
    qc = {
        'hostile': HOSTILE,
        'fastp': FASTP,
        'kraken': KRAKEN2,
        'bracken': BRACKEN
    }
    #############################
    ## assembly and annotation ##
    #############################
    SPADES = '''
    docker run -v $PWD:/data staphb/spades 
    spades.py 
    -1 reads/processed/{prefix}_R1_001.qcd.fastq.gz 
    -2 reads/processed/{prefix}_R2_001.qcd.fastq.gz 
    -o spades/{prefix} 
    --isolate 
    -t {threads} 
    -m {ram}
    '''.replace('\n', '').replace('    ', '')
    PROKKA = '''
    docker run -v $PWD:/data staphb/prokka:latest 
    prokka 
    spades/{prefix}/contigs.fasta 
    --outdir prokka/{prefix} 
    --prefix {prefix} 
    --genus Mycobacterium 
    --usegenus 
    --cpus {threads} 
    --force 
    --centre X 
    --compliant
    '''.replace('\n', '').replace('    ', '')
    assembly_and_annotation = {
        'spades': SPADES,
        'prokka': PROKKA
    }
    #########################
    ## population analysis ##
    #########################
    SKA_BUILD = '''
    ska build 
    -f {out_dir}/file_list.txt 
    -k 31 
    -o {out_dir}/index 
    --threads {threads} 
    -v
    '''.replace('\n', ' ').replace('    ', '')
    SKA_ALIGN = '''
    ska align 
    --min-freq 1 
    --filter no-filter 
    {out_dir}/index.skf 
    -o {out_dir}/alignment.aln 
    --ambig-mask 
    --threads {threads} 
    -v
    '''.replace('\n', ' ').replace('    ', '')
    IQTREE = 'iqtree -s {out_dir}/alignment.aln -m MFP --prefix iqtree/iqtree'
    SKA_MAP = '''
    ska map 
    -o {out_dir}/alignment.aln 
    --ambig-mask 
    --threads {threads} 
    {ref} 
    {out_dir}/index.skf
    '''.replace('\n', ' ').replace('    ', '')
    FASTBAPS = 'Rscript src/run_fastbaps.R && Rscript src/list_fastbaps_clusters.R'
    RUN_GUBBINS = '''
    run_gubbins.py 
    --prefix {out_dir}/gubbins 
    {out_dir}/alignment.aln 
    --threads 
    {threads}
    '''.replace('\n', ' ').replace('    ', '')
    SNP_DISTS = '''
    docker run -v $PWD:/data staphb/snp-dists 
    snp-dists 
    {out_dir}/gubbins.filtered_polymorphic_sites.fasta > {out_dir}/snp_dists.txt
    '''.replace('\n', ' ').replace('    ', '')
    population_analysis = {
        'ska_build': SKA_BUILD,
        'ska_align': SKA_ALIGN,
        'iqtree': IQTREE,
        'ska_map': SKA_MAP,
        'fastbaps': FASTBAPS,
        'run_gubbins': RUN_GUBBINS,
        'snp_dists': SNP_DISTS
    }
    ############
    ## MutTui ##
    ############
    MUTTUI = '''
    MutTui korimuto 
    -v muttui/vcfs/{query}.lofreq.tidy.vcf 
    -r muttui/fnas/{reference}.tidy.fna 
    -o muttui/output/{query} 
    --multi_contig
    '''.replace('\n', '').replace('    ', '')
    muttui = {'muttui': MUTTUI}
    #####################
    ## variant calling ##
    #####################
    BWA_INDEX = 'bwa index {index}'
    BWA_MEM = '''
    bwa mem -x intractg -t {threads} {index} reads/processed/{prefix}_R1_001.qcd.fastq.gz reads/processed/{prefix}_R2_001.qcd.fastq.gz | 
    samtools fixmate - - | 
    samtools view -bq 1 | 
    samtools sort -m 4g -o {out_dir}/{prefix}.bwa.bam -T {out_dir}/{prefix}.bwa.bam.tmp -
    '''.replace('\n', '').replace('    ', '')
    SAMTOOLS_INDEX_1 = 'samtools index {out_dir}/{prefix}.bwa.bam'
    BCFTOOLS_MPILEUP = '''
    bcftools mpileup -f {index} -Ou {out_dir}/{prefix}.bwa.bam | 
    bcftools call -Ob --ploidy 1 -v -m | 
    bcftools view - | 
    vcfutils.pl varFilter - > {out_dir}/{prefix}.mpileup.vcf
    '''.replace('\n', '').replace('    ', '')
    LOFREQ_PART_1 = '''
    lofreq viterbi -f {index} {out_dir}/{prefix}.bwa.bam | 
    lofreq indelqual --dindel -f {index} - | 
    samtools sort -m 4g -o {out_dir}/{prefix}.lofreq.bam -T {out_dir}/{prefix}.lofreq.bam.tmp -
    '''.replace('\n', '')
    SAMTOOLS_INDEX_2 = 'samtools index {out_dir}/{prefix}.lofreq.bam'
    LOFREQ_PART_2 = '''
    lofreq call-parallel 
    --sig 1E-4 
    --min-cov 25 
    --min-bq 25 
    --min-alt-bq 25 
    --min-mq 60 
    -d 10000 
    {bed}
    --pp-threads {threads} 
    -f {index} 
    -o {out_dir}/{prefix}.lofreq.vcf 
    {out_dir}/{prefix}.lofreq.bam
    '''.replace('\n', '').replace('    ', '')
    BGZIP = 'bgzip {out_dir}/{prefix}.lofreq.vcf'
    TABIX = 'tabix {out_dir}/{prefix}.lofreq.vcf.gz'
    variant_calling = {
        'bwa_index': BWA_INDEX,
        'bwa_mem': BWA_MEM,
        'samtools_index_1': SAMTOOLS_INDEX_1,
        'bcftools_mpileup': BCFTOOLS_MPILEUP,
        'lofreq_part_1': LOFREQ_PART_1,
        'samtools_index_2': SAMTOOLS_INDEX_2,
        'lofreq_part_2': LOFREQ_PART_2,
        'bgzip': BGZIP,
        'tabix': TABIX
    }
    ##############################
    ## create output dictionary ##
    ##############################
    dictionary = {
        'qc': qc,
        'assembly_and_annotation': assembly_and_annotation,
        'population_analysis': population_analysis,
        'muttui': muttui,
        'variant_calling': variant_calling
    }
    return(dictionary)

if __name__ == '__main__':
    dictionary = get_dictionary()
    with open('src/config.ini', 'w') as f:
        config = configparser.ConfigParser()
        for i in dictionary:
            config.add_section(i)
            for j in dictionary[i]:
                config.set(i, j, dictionary[i][j])
        config.write(f)
