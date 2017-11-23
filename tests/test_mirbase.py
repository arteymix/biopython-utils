from bioutil.mirbase import mature_accession_to_name, mature_name_to_accession
import pandas as pd

def test_mature_accession_to_name():
    assert 'hsa-miR-155-5p' == mature_accession_to_name('MIMAT0000646')
    assert 'hsa-miR-155-5p' in mature_accession_to_name(['MIMAT0000646'])
    assert 'hsa-miR-155-5p' in mature_accession_to_name(pd.Index(['MIMAT0000646']))
    assert 'hsa-miR-155-5p' in mature_accession_to_name(pd.Series(['MIMAT0000646'])).values

def test_mature_accession_to_name():
    assert 'MIMAT0000646' == mature_name_to_accession('hsa-miR-155-5p')
    assert 'MIMAT0000646' in mature_name_to_accession(['hsa-miR-155-5p'])
    assert 'MIMAT0000646' in mature_name_to_accession(pd.Index(['hsa-miR-155-5p']))
    assert 'MIMAT0000646' in mature_name_to_accession(pd.Series(['hsa-miR-155-5p'])).values
