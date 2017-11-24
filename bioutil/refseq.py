import pandas as pd
from Bio.GenBank import parse
from functools import wraps

def _subst(f):
    @wraps(f)
    def wrapper(*kwargs, **kwds):
        t, d = f(*kwargs, **kwds)
        if isinstance(t , str):
            return d[t]
        elif isinstance(t, pd.Index) or isinstance(t, pd.Series):
            return t.map(lambda x: d[x])
        else:
            return [d[x] for x in t]
    return wrapper

@_subst
def target_accession_to_version(index, rna_gbff):
    """Convert target accessions into versions"""
    v = {}
    with open(rna_gbff) as f:
        for record in parse(f):
            v[record.accession[0]] = record.version
    return index, v

@_subst
def target_accession_to_gene_name(index, rna_gbff):
    """Convert target accessions to their corresponding gene names"""
    v = {}
    with open(rna_gbff) as f:
        for record in parse(f):
            v[record.accession[0]] = next(q.value[1:-1]
                    for q in next(f.qualifiers
                        for f in record.features if f.key =='gene')
                            if q.key == '/gene=')
    return index, v

@_subst
def target_version_to_gene_name(index, rna_gbff):
    """Convert target versions (i.e. accession + '.' + version) to their corresponding gene names"""
    v = {}
    with open(rna_gbff) as f:
        for record in parse(f):
            v[record.version] = next(q.value[1:-1]
                    for q in next(f.qualifiers
                        for f in record.features if f.key =='gene')
                            if q.key == '/gene=')
    return index, v
