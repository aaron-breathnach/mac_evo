import argparse
import configparser
from math import floor
import os
import pandas as pd
from psutil import virtual_memory
from utils import stop_instance

cpu_count = os.cpu_count()

parser = argparse.ArgumentParser(prog='run_assembly_and_annotation.py', description='Does what it says on the tin')
parser.add_argument('--threads', help='Number of threads to use for alignment', default=cpu_count)
parser.add_argument('--cfg', help='Config file', default='src/config.ini')
parser.add_argument('--rerun', help='Rerun the analysis', action='store_true')
args = parser.parse_args()

def list_isolates():
    df = pd.read_csv('data/metadata.tsv', sep='\t')
    isolates = df[(df['time_from_diagnosis'] == 0) | (df['study'] == 'Present')]['isolate'].to_list()
    return(isolates)

def run_spades(prefix, threads, config):
    ram = floor(0.95 * virtual_memory().total / 1e9)
    cmd = config['assembly_and_annotation']['spades'].format(prefix=prefix, threads=threads, ram=ram)
    return(cmd)

def write_spades_shell_script(threads, config, rerun):
    cmds = []
    for prefix in list_isolates():
        target = 'spades/{}/contigs.fasta'.format(prefix)
        if not os.path.exists(target) or rerun:
            cmd = run_spades(prefix, threads, config)
            cmds.append(cmd)
    with open('run_spades.sh', 'w') as f:
        for cmd in cmds:
            f.write(cmd + '\n')

def run_prokka(prefix, threads, config):
    cmd = config['assembly_and_annotation']['prokka'].format(prefix=prefix, threads=threads)
    return(cmd)

def write_prokka_shell_script(threads, config, rerun):
    cmds = []
    for prefix in list_isolates():
        target = 'prokka/{}/{}.fna'.format(prefix, prefix)
        if not os.path.exists(target) or rerun:
            cmd = run_prokka(prefix, threads, config)
            cmds.append(cmd)
    with open('run_prokka.sh', 'w') as f:
        for cmd in cmds:
            f.write(cmd + '\n')

def run_quast(prefix, threads, config):
    cmd = config['assembly_and_annotation']['quast'].format(prefix=prefix, threads=threads)
    return(cmd)

def write_quast_shell_script(threads, config, rerun):
    cmds = []
    for prefix in list_isolates():
        target = 'quast/{}/report.tsv'.format(prefix)
        if not os.path.exists(target) or rerun:
            cmd = run_quast(prefix, threads, config)
            cmds.append(cmd)
    with open('run_quast.sh', 'w') as f:
        for cmd in cmds:
            f.write(cmd + '\n')

def write_checkm_shell_script(threads, config, rerun):
    target = 'checkm/checkm.tsv'
    if not os.path.isdir('checkm'):
        os.makedirs('checkm')
    with open('checkm/input.tsv', 'w') as f:
        for prefix in list_isolates():
            col1 = prefix
            col2 = 'spades/{}/contigs.fasta'.format(prefix)
            cols = '{}\t{}\n'.format(col1, col2)
            f.write(cols)
    with open('run_checkm.sh', 'w') as f:
        if not os.path.exists(target) or rerun:
            cmd = config['assembly_and_annotation']['checkm'].format(threads=threads)
            f.write(cmd + '\n')

def wrapper(cfg, threads, rerun):
    config = configparser.ConfigParser()
    config.read(cfg)
    write_spades_shell_script(threads, config, rerun)
    write_prokka_shell_script(threads, config, rerun)
    write_quast_shell_script(threads, config, rerun)
    write_checkm_shell_script(threads, config, rerun)
    stop = stop_instance()
    cmd = 'sh run_spades.sh && sh run_prokka.sh && sh run_quast.sh && sh run_checkm.sh && {}\n'.format(stop)
    with open('run_assembly_and_annotation.sh', 'w') as f:
        f.write(cmd)

if __name__ == '__main__':
    wrapper(cfg=args.cfg, threads=args.threads, rerun=args.rerun)
