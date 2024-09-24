import argparse
import os
import pandas as pd
import subprocess

cpu_count = os.cpu_count()

parser = argparse.ArgumentParser()
parser.add_argument('--out_dir', help='The output directory', default='panaroo')
parser.add_argument('--threads', help='Number of threads to use for alignment', default=cpu_count)
parser.add_argument('--rerun', help='Rerun the analysis', action='store_true')
args = parser.parse_args()

def get_ref(patient, df):
    df = df[df['patient'] == patient]
    reference = df.nsmallest(1, 'time_from_diagnosis').iloc[0, 1]
    query = df[df['isolate'] != reference]['isolate'].to_list()
    output = {'patient': patient, 'reference': reference, 'query': query}
    return(output)

def copy_gffs(directory):
    out_dir = '{}/input'.format(directory)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    dat = pd.read_csv('metadata.tsv', sep='\t')
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

def run_panaroo(directory, threads, rerun):
    copy_gffs(directory)
    target = '{}/output/gene_presence_absence.csv'.format(directory)
    if not os.path.exists(target) or rerun:
        panaroo = 'panaroo -i {}/input/*.gff -o {}/output -t {} --clean-mode sensitive\n'.format(directory, directory, threads)
        with open('run_within_host_panaroo.sh', 'w') as f:
            f.write(panaroo)

if __name__ == '__main__':
    run_panaroo(args.out_dir, args.threads, args.rerun)
