import argparse
import os
import pandas as pd
from run_lofreq import run_lofreq
import subprocess
from utils import get_ref

cpu_count = os.cpu_count()

parser = argparse.ArgumentParser(prog='run_lofreq_within_host.py', description='Run variant calling for each patient using LoFreq')
parser.add_argument('--metadata', help='Metadata file', default='metadata.tsv')
parser.add_argument('--config', help='Config file', default='src/config.ini')
parser.add_argument('--threads', help='Number of threads to use for alignment', default=cpu_count)
parser.add_argument('--rerun', help='Rerun the analysis', action='store_true')
args = parser.parse_args()

def copy(dictionary):
    reference = dictionary['reference']
    fna = 'prokka/{reference}/{reference}.fna'.format(reference=reference)
    out_dir = 'within_host_evolution/{}'.format(dictionary['patient'])
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists('{}/{}.fna'.format(out_dir, reference)):
        cp = 'cp {} {}'.format(fna, out_dir)
        subprocess.run(cp, shell=True)

def run_variant_calling(dictionary, config, threads, rerun):
    patient = dictionary['patient']
    out_dir = 'within_host_evolution/{}'.format(patient)
    reference = dictionary['reference']
    index = '{}/{}.fna'.format(out_dir, reference)
    queries = dictionary['query']
    cmds = []
    for query in queries:
        cmd = run_lofreq(config, query, out_dir, index, threads, False, rerun)
        cmds.append(cmd)
    return(cmds)

def run_lofreq_within_host(metadata, config, threads, rerun):
    dat = pd.read_csv(metadata, sep='\t')
    patients = list(set(dat['patient'].to_list()))
    cmds = []
    for patient in patients:
        dictionary = get_ref(patient, dat)
        copy(dictionary)
        cmd = run_variant_calling(dictionary, config, threads, rerun)
        cmds = cmds + cmd
    with open('run_lofreq_within_host.sh', 'w') as f:
        for cmd in cmds:
            if len(cmd.replace(' ', '')) > 0:
                f.write(cmd + '\n')

if __name__ == '__main__':
    run_lofreq_within_host(args.metadata, args.config, args.threads, args.rerun)
