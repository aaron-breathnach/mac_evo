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
    cmds = []
    if not os.path.exists(targets[0]) or rerun:
        bwa_index = config['variant_calling']['bwa_index'].format(index=index)
        cmds.append(bwa_index)
    if not os.path.exists(targets[1]) or rerun:
        bwa_mem = config['variant_calling']['bwa_mem'].format(threads=threads, index=index, prefix=prefix, out_dir=out_dir)
        cmds.append(bwa_mem)
    if not os.path.exists(targets[2]) or rerun:
        samtools_index_1 = config['variant_calling']['samtools_index_1'].format(out_dir=out_dir, prefix=prefix)
        cmds.append(samtools_index_1)
    if not os.path.exists(targets[3]) or rerun:
        lofreq_part_1 = config['variant_calling']['lofreq_part_1'].format(index=index, out_dir=out_dir, prefix=prefix)
        cmds.append(lofreq_part_1)
    if not os.path.exists(targets[4]) or rerun:
        samtools_index_2 = config['variant_calling']['samtools_index_2'].format(out_dir=out_dir, prefix=prefix)
        cmds.append(samtools_index_2)
    if not os.path.exists(targets[5]) and not os.path.exists(targets[3]) or rerun:
        if use_bed:
            bed = ' -l panaroo/output/single_copy_core_genes.bed '
        else:
            bed = ' '
        lofreq_part_2 = config['variant_calling']['lofreq_part_2'].format(bed=bed, threads=threads, index=index, out_dir=out_dir, prefix=prefix)
        cmds.append(lofreq_part_2)
    if not os.path.exists(targets[6]) or rerun:
        bgzip = config['variant_calling']['bgzip'].format(out_dir=out_dir, prefix=prefix)
        cmds.append(bgzip)
    if not os.path.exists(targets[7]) or rerun:
        tabix = config['variant_calling']['tabix'].format(out_dir=out_dir, prefix=prefix)
        cmds.append(tabix)
    rm = config['variant_calling']['rm'].format(out_dir=out_dir)
    cmds.append(rm)
    cmd = ' '.join(' && '.join(cmds).split())
    return(cmd)
