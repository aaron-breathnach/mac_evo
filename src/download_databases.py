import os
import pandas as pd
import subprocess

def download_type_strain():
    url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/741/445/GCF_009741445.1_ASM974144v1/GCF_009741445.1_ASM974144v1_genomic.fna.gz'
    out_dir = 'databases/GCF_009741445.1'
    filename = os.path.basename(url)
    target = '{}/{}'.format(out_dir, filename)
    if not os.path.exists(target):
        if not os.path.exists(target.rstrip('.gz')):
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            cmd = 'wget {} -O {} && gunzip {}'.format(url, target, target)
            subprocess.run(cmd, shell=True)

def download_kraken_database():
    url = 'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240605.tar.gz'
    out_dir = 'databases/kraken/k2_standard_08gb'
    filename = os.path.basename(url)
    target = '{}/{}'.format(out_dir, filename)
    if not os.path.exists(target):
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        cmd = 'wget {} -O {} && tar -xzvf {} -C {}'.format(url, target, target, out_dir)
        print(cmd)
        subprocess.run(cmd, shell=True)

def download_card_data():
    out_dir = 'data/card'
    filenames = ['shortname_pathogens.tsv', 'shortname_antibiotics.tsv', 'aro_index.tsv']
    targets = ['{}/{}'.format(out_dir, f) for f in filenames]
    if not all([os.path.exists(target) for target in targets]):
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        url = 'https://card.mcmaster.ca/download/0/broadstreet-v3.3.0.tar.bz2'
        tar_bz2 = 'data/card/{}'.format(os.path.basename(url))
        cmd_1 = 'wget {} -O {}'.format(url, tar_bz2)
        cmd_2 = 'tar -xvjf {} -C {}'.format(tar_bz2, out_dir)
        cmd = '{} && {}'.format(cmd_1, cmd_2)
        subprocess.run(cmd, shell=True)

def download_all_ref():
    if not os.path.exists('databases/reference_genomes'):
        os.makedirs('databases/reference_genomes')
    if not os.path.exists('data/assembly_summary.txt'):
        subprocess.run('wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -O data/assembly_summary.txt', shell=True)
    df = pd.read_csv('data/assembly_summary.txt', sep='\t', skiprows=1)
    subspecies = ['Mycobacterium avium subsp. avium', 'Mycobacterium avium subsp. hominissuis']
    ftp_paths = df[df['organism_name'].isin(subspecies)]['ftp_path'].to_list()
    for ftp_path in ftp_paths:
        basename = os.path.basename(ftp_path)
        url = '{}/{}_genomic.fna.gz'.format(ftp_path, basename)
        filename = 'databases/reference_genomes/{}'.format(os.path.basename(url))
        if not any([os.path.exists(x) for x in [filename, filename.replace('.gz', '')]]):
            subprocess.run('wget {} -O {} && gunzip {}'.format(url, filename, filename), shell=True)

def download_databases():
    download_type_strain()
    download_kraken_database()
    download_card_data()
    download_all_ref()

if __name__ == '__main__':
    download_databases()
