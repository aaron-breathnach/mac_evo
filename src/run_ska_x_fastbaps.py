import configparser
import os
import pandas as pd

## list fastbaps clusters with isolates from more than one patient
def list_fastbaps_clusters(metadata, fastbaps):
    df = pd.merge(metadata, fastbaps, on='isolate')
    df = df[['patient', 'level_1']].drop_duplicates()
    df = df.groupby('level_1').size().reset_index(name='n')
    clusters = df[df['n'] > 1]['level_1'].tolist()
    return(clusters)

## calculates within-cluster pairwise distances
def calc_within_cluster_dist(out_dir='ska_x_fastbaps', threads=os.cpu_count()):
    config = configparser.ConfigParser()
    config.read('src/config.ini')
    ska_build   = config['population_analysis']['ska_build']
    ska_map     = config['population_analysis']['ska_map']
    run_gubbins = config['population_analysis']['run_gubbins']
    snp_dists   = config['population_analysis']['snp_dists']
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    metadata = pd.read_csv('data/metadata.tsv', sep='\t')
    fastbaps = pd.read_csv('data/fastbaps.tsv', sep='\t')
    clusters = list_fastbaps_clusters(metadata, fastbaps)
    cmds = []
    for cluster in clusters:
        name = 'fastbaps_cluster_' + str(cluster).zfill(2)
        outdir = f'{out_dir}/{name}'
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        isolates = sorted(fastbaps[fastbaps['level_1'] == int(cluster)]['isolate'].to_list())
        targets = [f'{outdir}/index.skf',
                   f'{outdir}/alignment.aln',
                   f'{outdir}/gubbins.filtered_polymorphic_sites.fasta',
                   f'{outdir}/snp_dists.txt']
        if not os.path.exists(targets[0]):
            with open(f'{outdir}/file_list.txt', 'w') as file_list:
                for isolate in isolates:
                    reads = '\t'.join([f'reads/processed/{isolate}_R{i}_001.qcd.fastq.gz' for i in range(1, 3)])
                    row = f'{isolate}\t{reads}\n'
                    file_list.write(row)
            skabuild = ska_build.format(out_dir=outdir, threads=threads)
            cmds.append(skabuild)
        if not os.path.exists(targets[1]):
            ref = 'assemblies/{}.contigs.fasta'.format(isolates[0])
            skamap = ska_map.format(out_dir=outdir, threads=threads, ref=ref)
            cmds.append(skamap)
        if not os.path.exists(targets[2]):
            rungubbins = run_gubbins.format(out_dir=outdir, threads=threads)
            cmds.append(rungubbins)
        if not os.path.exists(targets[3]):
            snpdists = snp_dists.format(out_dir=outdir)
            cmds.append(snpdists)
    with open('run_ska_x_fastbaps.sh', 'w') as f:
        [f.write(cmd + '\n') for cmd in cmds]

if __name__ == '__main__':
    calc_within_cluster_dist()
