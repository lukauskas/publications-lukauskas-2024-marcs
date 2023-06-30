import mygene
from tqdm import tqdm
import numpy as np
import pandas as pd
import os
from datetime import datetime

INPUT_FILE = os.path.join(os.path.dirname(__file__), 'mygene.list')
OUTPUT_FILE = os.path.join(os.path.dirname(__file__), 'cache.h5cache')

# Make sure below matches snapanalysis.external.mygene
MYGENE_FIELDS = ['HGNC', 'entrezgene',
                 'symbol', 'alias', 'name', 'summary', 'ensembl.protein',
                 'interpro']


def _query_mygeneinfo(unique_protein_ids, fields=None):
    mg = mygene.MyGeneInfo()
    if fields is None:
        fields = MYGENE_FIELDS

    query_kwargs = dict(fields=fields,
                        scopes='uniprot',
                        species='human')

    ans = mg.querymany(unique_protein_ids,
                       as_dataframe=True,
                       **query_kwargs)

    # their API sometimes returns different result for POST queries as compared to GET
    # due to this some proteins are classified as not-found when they are not
    # let's query all the missing ones again to fill-in these cases
    ids_to_query_again = set()

    if 'notfound' in ans.columns:
        ids_to_query_again.update(ans[ans['notfound'].fillna(False)].index)

    # Similarly, due to GET/POST weirdness requery  all the genes that have no entrez id
    entrez_counts = ans.groupby(level=0)['entrezgene'].count()
    no_entrez = entrez_counts[entrez_counts == 0].index
    ids_to_query_again.update(no_entrez)

    requeried = []

    for notfound in tqdm(ids_to_query_again, desc='Querying missing data again'):
        response = mg.query(notfound,
                            as_dataframe=True,
                            **query_kwargs)

        response.index = [notfound] * len(response)

        requeried.append(response)

    requeried = pd.concat(requeried)

    # Drop the not-found rows from answer and replace them with the requeried data
    ans = ans.drop(requeried.index)
    ans = pd.concat([ans, requeried])

    ans.index.name = 'Protein ID'
    return ans

def populate_cache():
    """
    This populates mygene cache.
    We now cheat a little and provide a list of gene ids to fetch data for before
    we theoretically obtain them in the pipeline.

    The reasoning for this is to pre-cache all of the mygene responses as they tend to
    change every once in a while when proteins are updated.

    While this is great for having up-to-date results, it is a bit annoying sometimes
    to keep looking for a new label of protein
    """

    with open(INPUT_FILE) as f:
        protein_ids = f.readlines()

    protein_ids = [x.strip() for x in protein_ids if not x.startswith('#') and x.strip() != '']
    protein_ids = np.unique(protein_ids)
    mygene_ans = _query_mygeneinfo(protein_ids)

    with pd.HDFStore(OUTPUT_FILE, 'w') as store:
        store['query'] = pd.Series(protein_ids, name='Protein ID')
        store['cache'] = mygene_ans

        store['datetime'] = pd.Series([datetime.now().isoformat()])
        store['mygene_version'] = pd.Series([mygene.__version__])

if __name__ == '__main__':
    populate_cache()