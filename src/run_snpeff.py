import configparser
import os
import pandas as pd
from pathlib import Path
from utils import get_ref

info = '''
# Mycobacterium avium strain {reference}
{reference}.genome : Mycobacterium avium strain {reference}
    {reference}.chromosomes : {chromosomes}
    {reference}.{chromosome_0}.codonTable : Bacterial_and_Plant_Plastid
'''.strip('\n').split('\n')

def get_chromosomes(fna):
    with open(fna, 'r') as f:
        chromosomes = ', '.join([header.split('|')[-1:][0].strip('\n') for header in f.readlines() if '>' in header])
    return(chromosomes)

def edit_snpeff_config(reference, snpeff_config):
    target = '# Mycobacterium avium strain {}\n'.format(reference)
    f = open(snpeff_config, 'r')
    lines = f.readlines()
    if target in lines:
        print(f'Reference {reference} alreay in {snpeff_config}...')
    else:
        print(f'Adding reference {reference} to {snpeff_config}...')
        n = 0
        dictionary = {'part_1': [], 'part_3': []}
        for l in lines:
            n += 1
            if n < 235:
                dictionary['part_1'].append(l)
            else:
                dictionary['part_3'].append(l)
        fna = 'prokka/{reference}/{reference}.fna'.format(reference=reference)
        chromosomes = get_chromosomes(fna)
        chromosome_0 = chromosomes.split(', ')[0]
        part_2 = '\n'.join(info).format(reference=reference, chromosomes=chromosomes, chromosome_0=chromosome_0) + '\n'
        lines = dictionary['part_1'] + [part_2] + dictionary['part_3']
        with open(snpeff_config, 'w') as f:
            [f.write(line) for line in lines]

def snpEff_build(reference, snpeff_build):
    out_dir = '{}/snpEff/data/{}'.format(Path.home().as_posix(), reference)
    cmds = []
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(f'{out_dir}/cds.fa'):
        cds = 'cp prokka/{ref}/{ref}.ffn {out_dir}/cds.fa'.format(ref=reference, out_dir=out_dir)
        cmds.append(cds)
    if not os.path.exists(f'{out_dir}/protein.faa'):
        faa = 'cp prokka/{ref}/{ref}.faa {out_dir}/protein.faa'.format(ref=reference, out_dir=out_dir)
        cmds.append(faa)
    if not os.path.exists(f'{out_dir}/genes.gbk'):
        gbk = 'cp prokka/{ref}/{ref}.gbk {out_dir}/genes.gbk'.format(ref=reference, out_dir=out_dir)
        cmds.append(gbk)
    if not os.path.exists(f'{out_dir}/snpEffectPredictor.bin'):
        snpEff_build = snpeff_build.format(Path.home().as_posix(), reference)
        cmds.append(snpEff_build)
    cmd = ' && '.join(cmds)
    return(cmd)

def snpEff_main(inp_vcf, out_dir, reference, snpeff_main):
    if not inp_vcf.endswith('.gz'):
        cmd_1 = 'gsed -i \'s/gnl|X|//g\' {}'.format(inp_vcf)
    else:
        cmd_1 = 'zcat < {} | gsed \'s/gnl|X|//g\' > {}'.format(inp_vcf, inp_vcf.rstrip('.gz'))
        inp_vcf = inp_vcf.rstrip('.gz')
    out_vcf = '{}/{}'.format(out_dir, os.path.basename(inp_vcf).replace('.vcf', '.snpeff.vcf'))
    cmd_2 = snpeff_main.format(Path.home().as_posix(), reference, inp_vcf, out_vcf)
    cmd = f'{cmd_1} && {cmd_2}'
    return(cmd)

def run_snpEff(dictionary):
    config = configparser.ConfigParser()
    config.read('src/config.ini')
    snpeff_build = config['snpeff']['snpeff_build']
    snpeff_main  = config['snpeff']['snpeff_main']
    cmds = []
    reference = dictionary['reference']
    edit_snpeff_config(reference, '{}/snpEff/snpEff.config'.format(Path.home().as_posix()))
    cmds.append(snpEff_build(reference, snpeff_build))
    patient = dictionary['patient']
    queries = dictionary['query']
    for query in queries:
        out_dir = 'within_host_evolution/{}'.format(patient)
        inp_vcf = '{}/{}.lofreq.vcf'.format(out_dir, query)
        out_vcf = inp_vcf.strip('.vcf') + '.snpeff.vcf'
        if not os.path.exists(out_vcf):
            cmds.append(snpEff_main(inp_vcf, out_dir, reference, snpeff_main))
    cmd = ' && '.join([cmd for cmd in cmds if len(cmd) > 0])
    return(cmd)

def run_snpeff():
    dat = pd.read_csv('data/metadata.tsv', sep='\t')
    patients = list(set(dat['patient'].to_list()))
    cmds = []
    for patient in patients:
        dictionary = get_ref(patient, dat)
        cmd = run_snpEff(dictionary)
        cmds.append(cmd)
    with open('run_within_host_snpeff.sh', 'w') as f:
        for cmd in cmds:
            if len(cmd.replace(' ', '')) > 0:
                f.write(cmd + '\n')

if __name__ == '__main__':
    run_snpeff()
