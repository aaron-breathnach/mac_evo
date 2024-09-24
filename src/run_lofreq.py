import configparser
import os

def run_lofreq(configuration_file, prefix, out_dir, index, threads, use_bed, rerun):
    config = configparser.ConfigParser()
    config.read(configuration_file)
    targets = ['{index}.bwt',
               '{out_dir}/{prefix}.bwa.bam',
               '{out_dir}/{prefix}.bwa.bam.bai',
               '{out_dir}/{prefix}.lofreq.bam',
               '{out_dir}/{prefix}.lofreq.bam.bai',
               '{out_dir}/{prefix}.lofreq.vcf',
               '{out_dir}/{prefix}.lofreq.vcf.gz',
               '{out_dir}/{prefix}.lofreq.vcf.gz.tbi']
    targets = [target.format(prefix=prefix, out_dir=out_dir, index=index) for target in targets]
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    cmd = []
    if not os.path.exists(targets[0]) or rerun:
        cmd.append(config['variant_calling']['bwa_index'])
    if not os.path.exists(targets[1]) or rerun:
        cmd.append(config['variant_calling']['bwa_mem'])
    if not os.path.exists(targets[2]) or rerun:
        cmd.append(config['variant_calling']['samtools_index_1'])
    if not os.path.exists(targets[3]) or rerun:
        cmd.append(config['variant_calling']['lofreq_part_1'])
    if not os.path.exists(targets[4]) or rerun:
        cmd.append(config['variant_calling']['samtools_index_2'])
    if not os.path.exists(targets[5]) and not os.path.exists(targets[3]) or rerun:
        if use_bed:
            bed = ' -l panaroo/output/single_copy_core_genes.bed '
        else:
            bed = ' '
        cmd.append(config['variant_calling']['lofreq_part_2'].format(bed=bed, prefix=prefix, out_dir=out_dir, threads=threads, index=index))
    if not os.path.exists(targets[6]) or rerun:
        cmd.append(config['variant_calling']['bgzip'])
    if not os.path.exists(targets[7]) or rerun:
        cmd.append(config['variant_calling']['tabix'])
    cmd.append(config['variant_calling']['rm'])
    cmd = ' && '.join(cmd).format(prefix=prefix, out_dir=out_dir, threads=threads, index=index)
    return(cmd)