import os
import pandas as pd
from pathlib import Path
import subprocess

def run_mash():
    df = pd.read_csv('data/metadata.tsv', sep='\t')
    isolates = df[df['study'] == 'Present']['isolate'].to_list()
    queries = ['prokka/{}/{}.fna'.format(x, x) for x in isolates]
    references = [x.as_posix() for x in Path('databases/references_genomes').rglob('*')]
    genomes = queries + references
    if not os.path.exists('mash'):
        os.makedirs('mash')
    with open('mash/mash.list.txt', 'w') as f:
        for genome in genomes:
            f.write(genome + '\n')
    step_1 = 'docker run -v $PWD:/data staphb/mash mash sketch -l mash/mash.list.txt -o mash/mash.sketch'
    step_2 = 'docker run -v $PWD:/data staphb/mash mash dist mash/mash.sketch.msh mash/mash.sketch.msh -t > mash/mash.dist.txt'
    step_3 = 'cp mash/mash.dist.txt data/'
    mash = '{} && {} && {}'.format(step_1, step_2, step_3)
    subprocess.run(mash, shell=True)

if __name__ == '__main__':
    run_mash()
