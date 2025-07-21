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
    if not os.path.exists('ska'):
        os.makedirs('ska')
    config = configparser.ConfigParser()
    config.read('src/config.ini')
    ska_build = config['population_analysis']['ska_build'].format(out_dir='ska', threads=os.cpu_count())
    ska_align = config['population_analysis']['ska_align'].format(out_dir='ska', threads=os.cpu_count())
    fastbaps  = config['population_analysis']['fastbaps']
    snp_dists = 'docker run -v $PWD:/data staphb/snp-dists snp-dists ska/alignment.aln > data/snp_dists.species_wide.txt'
    make_file_list()
    with open('run_ska.sh', 'w') as f:
        f.write(ska_build + '\n')
        f.write(ska_align + '\n')
        f.write(fastbaps + '\n')
        f.write(snp_dists + '\n')

if __name__ == '__main__':
    run_ska()
