import argparse
import os
from run_lofreq import run_lofreq
import utils

cpu_count = os.cpu_count()

parser = argparse.ArgumentParser(prog='run_lofreq_common_reference.py', description='Run LoFreq with the type strain as the reference')
parser.add_argument('--metadata', help='Metadata file', default='metadata.tsv')
parser.add_argument('--config', help='Config file', default='src/config.ini')
parser.add_argument('--index', help='Reference genome', default='src/config.ini')
parser.add_argument('--threads', help='Number of threads to use for alignment', default=cpu_count)
parser.add_argument('--rerun', help='Rerun the analysis', action='store_true')
args = parser.parse_args()

if __name__ == '__main__':
    index = 'prokka/reference/reference.fna'
    with open('run_lofreq_common_reference.sh', 'w') as f:
        prefixes = utils.get_prefixes()
        cmds = []
        for prefix in prefixes:
            cmd = run_lofreq(configuration_file = args.config,
                             prefix=prefix,
                             out_dir='variant_calling',
                             index=index,
                             threads=args.threads,
                             use_bed=True,
                             rerun=args.rerun)
            f.write(cmd + '\n')
