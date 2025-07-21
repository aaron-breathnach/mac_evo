import os
import subprocess

def run_iqtree():
    target = 'data/iqtree.treefile'
    if not os.path.exists(target):
        if not os.path.exists('tmp'):
            os.makedirs('tmp')
        iqtree = 'iqtree -s ska/alignment.aln -m MFP --prefix tmp/iqtree'
        cp = 'cp tmp/iqtree.treefile data/iqtree.treefile'
        rm = 'rm -r tmp'
        cmd = f'{iqtree} && {cp} && {rm}'
        subprocess.run(cmd, shell=True)

if __name__ == '__main__':
    run_iqtree()
