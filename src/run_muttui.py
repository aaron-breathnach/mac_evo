import pandas as pd

muttui = '''
MutTui korimuto 
-v within_host_evolution/{patient}/{query}.lofreq.vcf 
-r within_host_evolution/{patient}/{reference}.fna 
-o within_host_evolution/{patient}/{query} 
--multi_contig
'''.replace('\n', '')

def run_muttui(metadata='data/metadata.tsv'):
    dat = pd.read_csv(metadata, sep='\t')
    patients = list(set(dat['patient'].to_list()))
    with open('run_muttui.sh', 'w') as f:
        for patient in patients:
            reference = dat[(dat['time_from_diagnosis'] == 0) & (dat['patient'] == patient)].iloc[0, 1]
            queries = dat[(dat['time_from_diagnosis'] > 0) & (dat['patient'] == patient)]['isolate'].to_list()
            for query in queries:
                cmd = muttui.format(patient=patient, reference=reference, query=query)
                f.write(cmd + '\n')

if __name__ == '__main__':
    run_muttui()
