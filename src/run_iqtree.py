import os
import configparser

def run_iqtree():
    if not os.path.exists('iqtree'):
        os.makedirs('iqtree')
    config = configparser.ConfigParser()
    config.read('src/config.ini')
    iqtree = config['population_analysis']['iqtree'].format(out_dir='ska')
    with open('run_iqtree.sh', 'w') as f:
        f.write(iqtree + '\n')

if __name__ == '__main__':
    run_iqtree()
