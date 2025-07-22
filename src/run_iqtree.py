import os

def run_iqtree():
    target = 'data/iqtree.treefile'
    if not os.path.exists(target):
        if not os.path.exists('tmp'):
            os.makedirs('tmp')
        iqtree = 'iqtree -s ska/alignment.aln -m MFP --prefix tmp/iqtree'
        cp = 'cp tmp/iqtree.treefile data/iqtree.treefile'
        rm = 'rm -r tmp'
        cmd = f'{iqtree} && {cp} && {rm}'
        with open('run_iqtree.sh', 'w') as f:
            f.write(cmd + '\n')

if __name__ == '__main__':
    run_iqtree()
