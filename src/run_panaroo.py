import configparser
import os
import pandas as pd
import subprocess
from utils import get_ref

def copy_gffs(directory):
    out_dir = '{}/input'.format(directory)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    dat = pd.read_csv('data/metadata.tsv', sep='\t')
    patients = list(set(dat['patient'].to_list()))
    gffs = []
    for patient in patients:
        dictionary = get_ref(patient, dat)
        ref = dictionary['reference']
        target = '{}/input/{}.gff'.format(directory, ref)
        if not os.path.exists(target):
            gff = 'prokka/{}/{}.gff'.format(ref, ref)
            gffs.append(gff)
    if len(gffs) > 0:
        cp = 'cp {} {}/'.format(' '.join(gffs), out_dir)
        subprocess.run(cp, shell=True)

def run_panaroo():
    copy_gffs('panaroo')
    target = 'panaroo/output/gene_presence_absence.csv'
    if not os.path.exists(target):
        config = configparser.ConfigParser()
        config.read('src/config.ini')
        threads = os.cpu_count()
        cmd_1 = config['panaroo']['panaroo'].format(threads=threads)
        cmd_2 = config['panaroo']['panaroo_generate_gffs']
        cmd_3 = config['panaroo']['panaroo_filter_pa']
        cmd_4 = config['panaroo']['parse']
        cmds = [cmd_1, cmd_2, cmd_3, cmd_4]
        for cmd in cmds:
            with open('run_within_host_panaroo.sh', 'w') as f:
                f.write(cmd + '\n')

if __name__ == '__main__':
    run_panaroo()
