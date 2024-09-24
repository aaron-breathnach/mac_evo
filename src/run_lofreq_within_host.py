import argparse
import os
import pandas as pd
import run_lofreq
import utils

cpu_count = os.cpu_count()

parser = argparse.ArgumentParser(prog='run_lofreq_within_host.py', description='Run variant calling for each patient using LoFreq')
parser.add_argument('--metadata', help='Metadata file', default='metadata.tsv')
parser.add_argument('--config', help='Config file', default='src/config.ini')
parser.add_argument('--threads', help='Number of threads to use for alignment', default=cpu_count)
parser.add_argument('--rerun', help='Rerun the analysis', action='store_true')
args = parser.parse_args()

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

if __name__ == '__main__':
    ## read metadata
    dat = pd.read_csv(args.metadata, sep='\t')
    patients = list(set(dat['patient'].to_list()))
    ## write commands for each patient
    cmds = []
    for patient in patients:
        dictionary = utils.get_ref(patient, dat)
        utils.copy(dictionary)
        cmd = run_variant_calling(dictionary, args.config, args.threads, args.rerun)
        cmds = cmds + cmd
    with open('run_lofreq_within_host.sh'.format(args.method), 'w') as f:
        for cmd in cmds:
            if len(cmd.replace(' ', '')) > 0:
                f.write(cmd + '\n')
