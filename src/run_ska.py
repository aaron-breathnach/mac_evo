import glob
import os
import configparser

def make_file_list():
    with open('ska/file_list.txt', 'w') as f:
        assemblies = [x for x in glob.glob('assemblies/*') if not 'Reference' in x]
        for assembly in assemblies:
            isolate = os.path.basename(assembly.split('.')[0])
            reads = '\t'.join(sorted(glob.glob(f'reads/processed/{isolate}*')))
            row = f'{isolate}\t{reads}\n'
            f.write(row)

def run_ska():
    if not os.path.exists('ska'):
        os.makedirs('ska')
    config = configparser.ConfigParser()
    config.read('src/config.ini')
    ska_build = config['population_analysis']['ska_build'].format(out_dir='ska', threads=os.cpu_count())
    ska_align = config['population_analysis']['ska_align'].format(out_dir='ska', threads=os.cpu_count())
    fastbaps = config['population_analysis']['fastbaps']
    make_file_list('ska')
    with open('run_ska.sh', 'w') as f:
        f.write(ska_build + '\n')
        f.write(ska_align + '\n')
        f.write(fastbaps + '\n')

if __name__ == '__main__':
    run_ska()
