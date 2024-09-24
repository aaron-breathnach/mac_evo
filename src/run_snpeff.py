import os
import pandas as pd
import sys

info = '''
# Mycobacterium avium strain {reference}
{reference}.genome : Mycobacterium avium strain {reference}
    {reference}.chromosomes : {chromosomes}
    {reference}.{chromosome_0}.codonTable : Bacterial_and_Plant_Plastid
'''.strip('\n').split('\n')

def get_chromosomes(fna):
    with open(fna, 'r') as f:
        chromosomes = ', '.join([header.split('|')[-1:][0].strip('\n') for header in f.readlines() if '>' in header])
    return(chromosomes)

def edit_config(reference):
    config = '/home/ec2-user/snpEff/snpEff.config'
    target = '# Mycobacterium avium strain {}'.format(reference)
    f = open(config, 'r')
    lines = f.readlines()
    if not target in lines:
        n = 0
        dictionary = {'part_1': [], 'part_3': []}
        for l in lines:
            n += 1
            if n < 235:
                dictionary['part_1'].append(l)
            else:
                dictionary['part_3'].append(l)
        fna = 'prokka/{reference}/{reference}.fna'.format(reference=reference)
        chromosomes = get_chromosomes(fna)
        chromosome_0 = chromosomes.split(', ')[0]
        part_2 = '\n'.join(info).format(reference=reference, chromosomes=chromosomes, chromosome_0=chromosome_0) + '\n'
        lines = dictionary['part_1'] + [part_2] + dictionary['part_3']
        with open(config, 'w') as f:
            [f.write(line) for line in lines]

def snpEff_build(reference):
    out_dir = '/home/ec2-user/snpEff/data/{}'.format(reference)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    cds = 'cp prokka/{ref}/{ref}.ffn {out_dir}/cds.fa'.format(ref=reference, out_dir=out_dir)
    faa = 'cp prokka/{ref}/{ref}.faa {out_dir}/protein.faa'.format(ref=reference, out_dir=out_dir)
    gbk = 'cp prokka/{ref}/{ref}.gbk {out_dir}/genes.gbk'.format(ref=reference, out_dir=out_dir)
    snpEff_build = 'java -jar /home/ec2-user/snpEff/snpEff.jar build -genbank -v {}'.format(reference)
    cmd = '{} && {} && {} && {}'.format(cds, faa, gbk, snpEff_build)
    return(cmd)

def snpEff_main(inp_vcf, out_dir, reference):
    if not inp_vcf.endswith('.gz'):
        cmd_1 = 'sed -i \'s/gnl|X|//g\' {}'.format(inp_vcf)
    else:
        cmd_1 = 'zcat {} | sed \'s/gnl|X|//g\' > {}'.format(inp_vcf, inp_vcf.rstrip('.gz'))
        inp_vcf = inp_vcf.rstrip('.gz')
    out_vcf = '{}/{}'.format(out_dir, os.path.basename(inp_vcf).replace('.vcf', '.snpeff.vcf'))
    cmd_2 = 'java -jar /home/ec2-user/snpEff/snpEff.jar {} {} > {}'.format(reference, inp_vcf, out_vcf)
    cmd = '{} && {}'.format(cmd_1, cmd_2)
    return(cmd)

def run_snpEff(dictionary, extension='vcf'):
    cmds = []
    patient = dictionary['patient']
    reference = dictionary['reference']
    ## make the database
    edit_config(reference)
    cmds.append(snpEff_build(reference))
    ## run snpEff
    queries = dictionary['query']
    for query in queries:
        out_dir = 'within_host_evolution/{}'.format(patient)
        inp_vcf = '{}/{}.{}'.format(out_dir, query, extension)
        cmds.append(snpEff_main(inp_vcf, out_dir, reference))
    cmd = ' && '.join(cmds)
    return(cmd)

def get_ref(patient, df):
    df = df[df['patient'] == patient]
    reference = df.nsmallest(1, 'time_from_diagnosis').iloc[0, 1]
    query = df[df['isolate'] != reference]['isolate'].to_list()
    output = {'patient': patient, 'reference': reference, 'query': query}
    return(output)

if __name__ == '__main__':
    dat = pd.read_csv(sys.argv[1], sep='\t')
    extension = sys.argv[2]
    patients = list(set(dat['patient'].to_list()))
    cmds = []
    for patient in patients:
        dictionary = get_ref(patient, dat)
        cmd = run_snpEff(dictionary, extension)
        cmds.append(cmd)
    with open('run_within_host_snpeff.sh', 'w') as f:
        for cmd in cmds:
            if len(cmd.replace(' ', '')) > 0:
                f.write(cmd + '\n')
