from Bio.SeqIO import parse
from io import BytesIO, TextIOWrapper
import pandas as pd
from urllib.request import urlopen
import gzip
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

def _retrieve_mature(version='CURRENT'):
    """Retrieve mature miRNAs"""
    mature_fa = urlopen('ftp://mirbase.org/pub/mirbase/{}/mature.fa.gz'.format(version))
    mature_fa = TextIOWrapper(gzip.GzipFile(fileobj=BytesIO(mature_fa.read())))
    return mature_fa

@_subst
def mature_accession_to_name(mature_accessions, mature_fa=None):
    """Convert mature miRNA accessions to their respective names"""
    if mature_fa is None:
        mature_fa = _retrieve_mature()
    mton = {record.description.split(' ')[1]: record.id
                for record in parse(mature_fa, 'fasta')}
    return mature_accessions, mton

@_subst
def mature_name_to_accession(mature_names, mature_fa=None):
    """Convert a miRNA names to their respective accessions"""
    if mature_fa is None:
        mature_fa = _retrieve_mature()
    ntom = {record.id: record.description.split(' ')[1] for record in parse(mature_fa, 'fasta')}
    return mature_names, ntom
