import argparse
import configparser
import os
from utils import get_prefixes

cpu_count = os.cpu_count()

parser = argparse.ArgumentParser(prog='run_qc.py', description='Does what it says on the tin')
parser.add_argument('--threads', help='Number of threads to use', default=cpu_count)
parser.add_argument('--cfg', help='Config file', default='src/config.ini')
parser.add_argument('--run_kraken_bracken', help='Run Kraken + Bracken', action='store_true')
parser.add_argument('--rerun', help='Rerun the analysis', action='store_true')
args = parser.parse_args()

def run_qc(prefix, threads, config, run_kraken_bracken, rerun):
    output = {
        'hostile': 'reads/hostile/{}_R1_001.clean_1.fastq.gz'.format(prefix),
        'fastp': 'reads/processed/{}_R1_001.qcd.fastq.gz'.format(prefix),
        'kraken': 'kraken/{}.output'.format(prefix),
        'bracken': 'kraken/{}.bracken'.format(prefix)
    }
    cmd = []
    if not os.path.exists(output['fastp']) or rerun:
        if not os.path.exists(output['hostile']) or rerun:
            HOSTILE = config['qc']['hostile']
            cmd_1 = HOSTILE.format(prefix=prefix, threads=threads)
            cmd.append(cmd_1)
        FASTP = config['qc']['fastp']
        cmd_2 = '{} && rm reads/raw/{}* reads/hostile/{}*'.format(FASTP.format(prefix=prefix), prefix, prefix)
        cmd.append(cmd_2)
    if run_kraken_bracken:
        if not os.path.exists('kraken'):
            os.makedirs('kraken')
        if not os.path.exists(output['kraken']) or rerun:
            KRAKEN2 = config['qc']['kraken']
            cmd_3 = KRAKEN2.format(prefix=prefix, threads=threads)
            cmd.append(cmd_3)
        if not os.path.exists(output['bracken']) or rerun:
            BRACKEN = config['qc']['bracken']
            cmd_4 = '{} && rm kraken/{}.output'.format(BRACKEN.format(prefix=prefix), prefix)
            cmd.append(cmd_4)
    cmd = ' && '.join(cmd)
    return(cmd)

def write_qc_shell_script(threads, cfg, run_kraken_bracken, rerun):
    config = configparser.ConfigParser()
    config.read(cfg)
    cmds = [run_qc(prefix, threads, config, run_kraken_bracken, rerun) for prefix in get_prefixes()]
    with open('run_qc.sh', 'w') as f:
        for cmd in cmds:
            if len(cmd) > 0:
                f.write(cmd + '\n')

if __name__ == '__main__':
    write_qc_shell_script(
        threads=args.threads,
        cfg=args.cfg,
        run_kraken_bracken=args.run_kraken_bracken,
        rerun=args.rerun
    )
