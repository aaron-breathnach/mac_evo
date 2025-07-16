import argparse
import configparser
import os
from utils import get_prefixes

parser = argparse.ArgumentParser(prog='run_snippy_core.py', description='Does what it says on the tin')
parser.add_argument('--cfg', help='Config file', default='src/config.ini')
parser.add_argument('--rerun', help='Rerun the analysis', action='store_true')
args = parser.parse_args()

if __name__ == '__main__':
    config = configparser.ConfigParser()
    config.read(args.cfg)
    with open('run_snippy_core.sh', 'w') as f:
        cmds = []
        targets = ['snippy/core.full.aln', 'snippy/gubbins.aln']
        if not os.path.exists(targets[0]) or args.rerun:
            directories = ' '.join(['snippy/{}'.format(prefix) for prefix in get_prefixes()])
            SNIPPY_CORE = config['snippy']['snippy_core'].format(directories)
            cmds.append(SNIPPY_CORE)
        if not os.path.exists(targets[1]) or args.rerun:
            SNIPPY_CLEAN_FULL_ALN = config['snippy']['snippy_clean_full_aln']
            cmds.append(SNIPPY_CLEAN_FULL_ALN)
        cmd = ' && '.join(cmds)
        f.write(cmd)
