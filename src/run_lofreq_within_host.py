import argparse
import os
import pandas as pd
from run_lofreq import run_lofreq
from utils import get_ref

cpu_count = os.cpu_count()

parser = argparse.ArgumentParser(prog='run_lofreq_within_host.py', description='Run variant calling for each patient using LoFreq')
parser.add_argument('--metadata', help='Metadata file', default='data/metadata.tsv')
parser.add_argument('--config', help='Config file', default='src/config.ini')
parser.add_argument('--threads', help='Number of threads to use for alignment', default=cpu_count)
parser.add_argument('--rerun', help='Rerun the analysis', action='store_true')
args = parser.parse_args()

def run_lofreq_within_host(metadata, config, threads, rerun):
    dat = pd.read_csv(metadata, sep='\t')
    patients = list(set(dat['patient'].to_list()))
    cmds = []
    for patient in patients:
        out_dir = 'within_host_evolution/{}'.format(patient)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        dictionary = get_ref(patient, dat)
        reference = dictionary['reference']
        ## copy reference genome to output directory
        old = 'prokka/{reference}/{reference}.fna'.format(reference=reference)
        new = '{}/{}.fna'.format(out_dir, reference)
        if not os.path.exists(new):
            cp = 'cp {} {}'.format(old, new)
            cmds.append(cp)
        ## run LoFreq for each query
        queries = dictionary['query']
        for query in queries:
            variant_calling = run_lofreq(config, query, out_dir, new, threads, False, rerun)
            cmds.append(variant_calling)
        with open('run_lofreq_within_host.sh', 'w') as f:
            for cmd in cmds:
                f.write(cmd + '\n')

if __name__ == '__main__':
    run_lofreq_within_host(args.metadata, args.config, args.threads, args.rerun)
