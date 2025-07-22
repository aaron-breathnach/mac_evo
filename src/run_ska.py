import configparser
import glob
import os
from utils import get_prefixes

def make_file_list():
    with open('ska/file_list.txt', 'w') as f:
        isolates = get_prefixes()
        for isolate in isolates:
            reads = '\t'.join(sorted(glob.glob(f'reads/processed/{isolate}*')))
            row = f'{isolate}\t{reads}\n'
            f.write(row)

def run_ska():
    targets = ['ska/index.skf', 'ska/alignment.aln', 'data/fastbaps.tsv', 'data/snp_dists.species_wide.txt']
    if not os.path.exists('ska'):
        os.makedirs('ska')
    config = configparser.ConfigParser()
    config.read('src/config.ini')
    cmds = []
    if not os.path.exists(targets[0]):
        make_file_list()
        ska_build = config['population_analysis']['ska_build'].format(out_dir='ska', threads=os.cpu_count())
        cmds.append(ska_build)
    if not os.path.exists(targets[1]):
        ska_align = config['population_analysis']['ska_align'].format(out_dir='ska', threads=os.cpu_count())
        cmds.append(ska_align)
    if not os.path.exists(targets[2]):
        fastbaps = config['population_analysis']['fastbaps']
        cmds.append(fastbaps)
    if not os.path.exists(targets[3]):
        snp_dists = 'docker run -v $PWD:/data staphb/snp-dists snp-dists ska/alignment.aln > data/snp_dists.species_wide.txt'
        cmds.append(snp_dists)
    with open('run_ska.sh', 'w') as f:
        for cmd in cmds:
            f.write(cmd + '\n')

if __name__ == '__main__':
    run_ska()
