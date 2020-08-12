import os

import pandas as pd
from snapanalysis.config import EXTERNAL_DATA_DIRECTORY

# Make sure below matches data/downloaded/mygene/populate_cache.py
# as well as data/downloaded makefile
CACHE_FILE = os.path.join(EXTERNAL_DATA_DIRECTORY, 'mygene/cache.h5cache')
_MYGENE_PROTEIN_ID_FILE = os.path.join(EXTERNAL_DATA_DIRECTORY, 'mygene/mygene.list')

# Make sure below matches data/downloaded/mygene/populate_cache.py
MYGENE_FIELDS = ['HGNC', 'entrezgene',
                 'symbol', 'alias', 'name', 'summary', 'ensembl.protein',
                 'interpro']

def query_mygene_from_cache(unique_protein_ids):

    protein_id_set = set(unique_protein_ids)

    with pd.HDFStore(CACHE_FILE, 'r') as store:
        query = store['query']

        query_set = set(query)

        difference = protein_id_set.difference(query_set)

        if difference:
            raise Exception('Protein IDs {!r} not in cache'.format(difference))

        cache = store['cache']

    subset = cache[[ix in protein_id_set for ix in cache.index]]
    return subset