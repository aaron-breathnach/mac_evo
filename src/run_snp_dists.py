import argparse
import os
import subprocess

cpu_count = os.cpu_count()

parser = argparse.ArgumentParser(prog='run_snp_dists.py', description='Run snp-dists on whole genome alignment')
parser.add_argument('--inp', help='Alignment generated by Snippy', default='snippy/core.aln')
parser.add_argument('--out_dir', help='Output directory', default='data')
parser.add_argument('--threads', help='Number of threads to use', default=cpu_count)
parser.add_argument('--rerun', help='Rerun the analysis', action='store_true')
args = parser.parse_args()

snp_dists = 'docker run -v $PWD:/data staphb/snp-dists snp-dists -j {} {} > {}/snp_dists.tsv'

def run_snp_dists(inp, out_dir, threads, rerun):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    target = '{}/snp_dists.tsv'.format(out_dir)
    if not os.path.exists(target) or rerun:
        cmd = snp_dists.format(threads, inp, out_dir)
        subprocess.run(cmd, shell=True)

if __name__ == '__main__':
    run_snp_dists(args.inp, args.out_dir, args.threads, args.rerun)