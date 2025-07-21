import pandas as pd

def get_prefixes():
    prefixes = pd.read_csv('data/metadata.tsv', sep='\t')['isolate'].tolist()
    return(prefixes)

def get_ref(patient, df):
    df = df[df['patient'] == patient]
    reference = df.nsmallest(1, 'time_from_diagnosis').iloc[0, 1]
    query = df[df['isolate'] != reference]['isolate'].to_list()
    output = {'patient': patient, 'reference': reference, 'query': query}
    return(output)
