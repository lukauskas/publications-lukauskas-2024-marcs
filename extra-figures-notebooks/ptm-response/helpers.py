import pandas as pd
import numpy as np
from functools import lru_cache

from snapanalysis.models.ptm_response.main import OUTPUT_FILE as PTM_RESPONSE_FILE
from snapanalysis.models.ptm_response.main import PREDICTOR_ORDER

def load_limma_data(predictor=None, dropout_predictor=None):
    
    if dropout_predictor is None:
        df = pd.read_hdf(PTM_RESPONSE_FILE, '/ptm_stats/joint_limma_stats')
        if predictor is not None:
            return df.loc[predictor]
        else:
            return df
    else:
        return pd.read_hdf(PTM_RESPONSE_FILE, f'/ptm_stats/{predictor}/dropout/{dropout_predictor}/limma/stats')

def load_limma_design(predictor, dropout_predictor=None):
    
    if dropout_predictor is None:
        return pd.read_hdf(PTM_RESPONSE_FILE, f'/ptm_stats/{predictor}/limma/design')
    else:
        return pd.read_hdf(PTM_RESPONSE_FILE, f'/ptm_stats/{predictor}/dropout/{dropout_predictor}/limma/design')

def load_limma_matrix(predictor, dropout_predictor=None):
    
    if dropout_predictor is None:
        return pd.read_hdf(PTM_RESPONSE_FILE, f'/ptm_stats/{predictor}/limma/matrix')
    else:
        return pd.read_hdf(PTM_RESPONSE_FILE, f'/ptm_stats/{predictor}/dropout/{dropout_predictor}/limma/matrix')

def load_limma_long_matrix(predictor):
    return pd.read_hdf(PTM_RESPONSE_FILE, f'/ptm_stats/{predictor}/long_matrix')
    

# --- complexes --- 
from snapanalysis.preprocessing.protein_metadata import get_complex_memberships, get_complexes_as_dict

_COMPLEX_MEMBERSHIPS = get_complex_memberships()
_COMPLEX_MEMBERSHIPS_EXCLUSIVE_DICT = get_complexes_as_dict(exclusive_subunits_only=True)

def members_of(complexes):
    if isinstance(complexes, str):
        complexes = [complexes]
    ans = set()
    for complex_ in complexes:
        
        if complex_.endswith(' (exclusive subunits)'):
            real_complex = complex_[:-len(' (exclusive subunits)')]
            members = _COMPLEX_MEMBERSHIPS_EXCLUSIVE_DICT[real_complex]
        else:
            members = _COMPLEX_MEMBERSHIPS[_COMPLEX_MEMBERSHIPS==complex_].reset_index()['Gene label']
            
        if not len(members):
            raise Exception(f'Found no members for {complex_!r}')
        ans.update(members)
    
    return sorted(ans)

# --- URLs ---
from urllib.parse import quote

from snapanalysis.preprocessing.pulldown_metadata import OUTPUT_FILE as META_FILE
from snapanalysis.models.ptm_response.predictor_graph import EDGE_SEPARATOR, SPECIAL_PULLDOWN

def parse_edge(edge):
    left, __, right = edge.partition(EDGE_SEPARATOR)
    
    
    if right == SPECIAL_PULLDOWN:
        return [left]
    else:
        # Follow reverse convention here
        return [right, left]

@lru_cache(20)
def predictor_edges(predictor):
    
    pd_predictor_matrix = pd.read_hdf(META_FILE, '/meta/predictors')
    n_predictors = pd_predictor_matrix.sum(axis=1)

    lm = load_limma_long_matrix(predictor)
    edges = lm.reset_index()['edge'].unique()
    edges = [parse_edge(edge) for edge in edges]
    edges = sorted(edges, 
                   key=lambda edge: (len(edge), 
                                     min([n_predictors.loc[x] for x in edge])))

    return edges

def predictor_sorted_uri(predictor, proteins=None, extra_pds=None):        
    edges = predictor_edges(predictor)
    
    pds_flattened = []
    for e in edges:
        for pd in e:
            pds_flattened.append(pd)
    
    if extra_pds is not None:
        extra_pds = [p for p in extra_pds if p not in pds_flattened]
        pds_flattened.extend(extra_pds)
   
    pd_str = ','.join(pds_flattened)
    
    prot_str = ''
    if proteins:
        if isinstance(proteins, str):
            proteins = [proteins]
            
        prot_str = '&'.join([f'k=p:{x}' for x in map(quote, proteins)])
        prot_str = '&' + prot_str
    
    uri = f'http://ife-snap-data/proteins?pdorder={pd_str}{prot_str}&showsimilar=false&noclusterproteins=true'
    return uri

