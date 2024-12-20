import argparse
import os
import subprocess

cpu_count = os.cpu_count()

parser = argparse.ArgumentParser(prog='run_run_gubbins.py', description='Build a phylogenetic tree with Gubbins')
parser.add_argument('--inp', help='Whole genome alignment generated by Snippy', default='snippy/gubbins.aln')
parser.add_argument('--prefix', help='Prefix for the output files', default='gubbins/mavium')
parser.add_argument('--threads', help='Number of threads to use', default=cpu_count)
parser.add_argument('--rerun', help='Rerun the analysis', action='store_true')
args = parser.parse_args()

def run_gubbins(inp, prefix, threads, rerun):
    target = '{}.final_tree.tre'.format(prefix)
    if not os.path.exists(target) or rerun:
        cmd = 'run_gubbins.py {} --prefix {} --threads {} --tree-builder iqtree-fast'.format(inp, prefix, threads)
        subprocess.run(cmd, shell=True)

if __name__ == '__main__':
    run_gubbins(args.inp, args.prefix, args.threads, args.rerun)
