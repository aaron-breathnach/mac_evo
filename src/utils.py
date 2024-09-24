from ec2_metadata import ec2_metadata
import os
import pandas as pd
from pathlib import Path
import subprocess
import itertools

def get_prefixes():
    prefixes = pd.read_csv('metadata.tsv', sep='\t')['isolate'].tolist()
    return(prefixes)

def get_ref(patient, df):
    df = df[df['patient'] == patient]
    reference = df.nsmallest(1, 'time_from_diagnosis').iloc[0, 1]
    query = df[df['isolate'] != reference]['isolate'].to_list()
    output = {'patient': patient, 'reference': reference, 'query': query}
    return(output)

def get_reads(prefix):
    r1 = 'reads/processed/{}_R1_001.qcd.fastq.gz'.format(prefix)
    r2 = 'reads/processed/{}_R2_001.qcd.fastq.gz'.format(prefix)
    reads = [r1, r2]
    return(reads)

def copy(dictionary):
    reference = 'prokka/{}/{}.fna'.format(dictionary['reference'], dictionary['reference'])
    reads = [get_reads(query) for query in dictionary['query']]
    files = ' '.join([reference] + list(itertools.chain.from_iterable(reads)))
    out_dir = 'within_host_evolution/{}/'.format(dictionary['patient'])
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    cp = 'cp {} {}'.format(files, out_dir)
    subprocess.run(cp, shell=True)

def stop_instance():
    instance_id = ec2_metadata.instance_id
    cmd = 'aws ec2 stop-instances --instance-ids {}'.format(instance_id)
    return(cmd)