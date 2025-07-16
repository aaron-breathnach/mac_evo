import argparse
import configparser
import os
from utils import get_prefixes

cpu_count = os.cpu_count()

parser = argparse.ArgumentParser(prog='run_assembly_and_annotation.py', description='Assembles genomes with SPAdes and annotates them with Prokka')
parser.add_argument('--threads', help='Number of threads to use', default=cpu_count)
parser.add_argument('--cfg', help='Config file', default='src/config.ini')
parser.add_argument('--rerun', help='Rerun the analysis', action='store_true')
args = parser.parse_args()

if __name__ == '__main__':
    prefixes = get_prefixes()
    config = configparser.ConfigParser()
    config.read('src/config.ini')
    SNIPPY = config['snippy']['snippy']
    with open('run_snippy.sh', 'w') as f:
        for prefix in prefixes:
            target = 'snippy/{}/snps.tab'.format(prefix)
            if not os.path.exists(target):
                cmd = SNIPPY.format(prefix=prefix, threads=args.threads)
                f.write(cmd + '\n')
