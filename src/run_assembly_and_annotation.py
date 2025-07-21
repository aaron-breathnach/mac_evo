import argparse
import configparser
from math import floor
import os
import pandas as pd
from psutil import virtual_memory

cpu_count = os.cpu_count()

parser = argparse.ArgumentParser(prog='run_assembly_and_annotation.py', description='Does what it says on the tin')
parser.add_argument('--threads', help='Number of threads to use for alignment', default=cpu_count)
parser.add_argument('--cfg', help='Config file', default='src/config.ini')
parser.add_argument('--rerun', help='Rerun the analysis', action='store_true')
args = parser.parse_args()

def list_isolates():
    df = pd.read_csv('data/metadata.tsv', sep='\t')
    isolates = df[df['time_from_diagnosis'] == 0]['isolate'].to_list()
    return(isolates)

def run_assembly_and_annotation(threads, cfg, rerun):
    config = configparser.ConfigParser()
    config.read(cfg)
    ram = floor(0.95 * virtual_memory().total / 1e9)
    cmds = []
    for prefix in list_isolates():
        targets = ['spades/{}/contigs.fasta'.format(prefix), 'prokka/{}/{}.fna'.format(prefix, prefix)]
        if not os.path.exists(targets[0]) or rerun:
            cmd = config['assembly_and_annotation']['spades'].format(prefix=prefix, threads=threads, ram=ram)
            cmds.append(cmd)
        if not os.path.exists(targets[1]) or rerun:
            cmd = config['assembly_and_annotation']['prokka'].format(prefix=prefix, threads=threads)
            cmds.append(cmd)
    with open('run_assembly_and_annotation.sh', 'w') as f:
        for cmd in cmds:
            f.write(cmd + '\n')

if __name__ == '__main__':
    run_assembly_and_annotation(args.threads, args.cfg, args.rerun)
