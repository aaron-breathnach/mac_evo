from Bio import SeqIO
import configparser
import os
import pandas as pd
from pathlib import Path
import subprocess

def tidy_fasta(fasta):
    if not os.path.exists('muttui/fnas'):
        os.makedirs('muttui/fnas')
    tidy_fasta = 'muttui/fnas/{}'.format(os.path.basename(fasta.replace('.fna', '.tidy.fna')))
    with open(tidy_fasta, 'w') as f:
        for record in SeqIO.parse(fasta, "fasta"):
            header = record.id.split('|')[2]
            sequence = str(record.seq)
            f.write(f'>{header}\n{sequence}\n')

def tidy_vcf(vcf):
    if not os.path.exists('muttui/vcfs'):
        os.makedirs('muttui/vcfs')
    out = 'muttui/vcfs/{}'.format(vcf.split("/")[-1].replace(".vcf.gz", ".tidy.vcf"))
    gunzip = f'gunzip -c {vcf} > {out}'
    gsed = f'gsed -i "s/gnl|X|//" {out}'
    cmd = f'{gunzip} && {gsed}'
    subprocess.run(cmd, shell=True)

def run_muttui():
    config = configparser.ConfigParser()
    config.read('src/config.ini')
    muttui = config['muttui']['muttui']
    fastas = [x.as_posix() for x in Path('within_host_evolution').rglob('*.fna') if 'tidy' not in x.name]
    for fasta in fastas:
        tidy_fasta(fasta)
    vcfs = [x.as_posix() for x in Path('within_host_evolution').rglob('*.vcf.gz')] 
    for vcf in vcfs:
        tidy_vcf(vcf)
    if not os.path.exists('muttui/output'):
        os.makedirs('muttui/output')
    dat = pd.read_csv('data/metadata.tsv', sep='\t')
    dat = dat[dat['bracken_pass']]
    dat = dat[dat['multiple_carriage'] == False]
    patients = list(set(dat['patient'].to_list()))
    with open('run_muttui.sh', 'w') as f:
        for patient in patients:
            reference = dat[(dat['time_from_diagnosis'] == 0) & (dat['patient'] == patient)].iloc[0, 1]
            queries = dat[(dat['time_from_diagnosis'] > 0) & (dat['patient'] == patient)]['isolate'].to_list()
            for query in queries:
                cmd = muttui.format(query=query, reference=reference)
                f.write(cmd + '\n')

if __name__ == '__main__':
    run_muttui()
