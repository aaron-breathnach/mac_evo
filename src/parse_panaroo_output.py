import networkx

def get_core_gets(gene_presense_absence):
    f = open(gene_presense_absence, 'r')
    core_genes = []
    n = 0
    for l in f.readlines():
        n += 1
        if n > 1:
            cols = l.split('\t')
            gene = cols[0]
            vals = [int(x.replace('\n', '')) for x in cols[1:len(cols)]]
            if sum(vals) == len(cols) - 1:
                core_genes.append(gene)
    return(core_genes)

def get_single_copy_core_genes(core_genes, final_graph):
    g = networkx.read_gml(final_graph)
    paralogs = networkx.get_node_attributes(g, 'paralog')
    names = networkx.get_node_attributes(g, 'name')
    single_copy_core_genes = []
    for p in paralogs:
        if paralogs[p] == 0:
            if names[p] in core_genes:
                single_copy_core_genes.append(p)
    return(single_copy_core_genes)

def get_pangenome_id(l):
    pangenome_id = l.split('\t')[8].split(';')[4].replace('pangenome_id=', '')
    return(pangenome_id)

def make_bed(l):
    cols = l.split('\t')
    chrom = cols[0]
    chromStart = cols[3]
    chromEnd = cols[4]
    name = 'pangenome_id={}'.format(get_pangenome_id(l))
    bed = '\t'.join([chrom, chromStart, chromEnd, name])
    return(bed)

def parse_panaroo_output(gene_presence_absence, final_graph, gff, bed):
    core_genes = get_core_gets(gene_presence_absence)
    single_copy_core_genes = get_single_copy_core_genes(core_genes, final_graph)
    with open(gff, 'r') as f:
        lines = [make_bed(l) for l in f.readlines() if 'pangenome_id' in l and get_pangenome_id(l) in single_copy_core_genes]
    with open(bed, 'w') as f:
        [f.write(l + '\n') for l in lines]

if __name__ == 'main':
    gene_presence_absence = 'panaroo/results/gene_presence_absence_filt_pseudo_length.Rtab'
    final_graph = 'panaroo/results/final_graph.gml'
    gff = 'panaroo/results/postpanaroo_gffs/reference_panaroo.gff'
    bed = 'panaroo/results/single_copy_core_genes.bed'
    parse_panaroo_output(gene_presence_absence, final_graph, gff, bed)
