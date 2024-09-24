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
    QUAST = '''
    docker run -v $PWD:/data staphb/quast 
    quast.py 
    prokka/{prefix}/{prefix}.fna 
    -r databases/GCF_009741445.1/GCF_009741445.1_ASM974144v1_genomic.fna 
    -g databases/GCF_009741445.1/genomic.gff 
    -1 reads/processed/{prefix}_R1_001.qcd.fastq.gz 
    -2 reads/processed/{prefix}_R1_001.qcd.fastq.gz 
    -o quast/{prefix} 
    -t {threads}
    '''.replace('\n', '').replace('    ', '')
    CHECKM = '''
    docker run -v $PWD:/data staphb/checkm:latest 
    checkm 
    lineage_wf 
    -t {threads} 
    checkm/input.txt 
    checkm 
    -f checkm/checkm.tsv 
    --tab
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
        'quast': QUAST,
        'checkm': CHECKM,
        'prokka': PROKKA
    }
    ############
    ## snippy ##
    ############
    SNIPPY = '''
    docker run -v $PWD:/data staphb/snippy 
    snippy 
    --cpus {threads} 
    --outdir snippy/{prefix} 
    --ref databases/GCF_009741445.1/GCF_009741445.1_ASM974144v1_genomic.fna 
    --R1 reads/processed/{prefix}_R1_001.qcd.fastq.gz 
    --R2 reads/processed/{prefix}_R2_001.qcd.fastq.gz
    '''.replace('\n', '').replace('    ', '')
    SNIPPY_CORE = '''
    docker run -v $PWD:/data staphb/snippy 
    snippy-core 
    --ref databases/GCF_009741445.1/GCF_009741445.1_ASM974144v1_genomic.fna 
    --prefix snippy/core 
    {}
    '''.replace('\n', '').replace('    ', '')
    SNIPPY_CLEAN_FULL_ALN = 'docker run -v $PWD:/data staphb/snippy snippy-clean_full_aln snippy/core.full.aln > snippy/gubbins.aln'
    snippy = {
        'snippy': SNIPPY,
        'snippy_core': SNIPPY_CORE,
        'snippy_clean_full_aln': SNIPPY_CLEAN_FULL_ALN
    }
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
        'snippy': snippy,
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
