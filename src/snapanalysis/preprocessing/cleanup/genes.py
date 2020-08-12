import pandas as pd
from snapanalysis.config import get_logger
from snapanalysis.external.metadata.mygene import MYGENE_FIELDS, query_mygene_from_cache
from snapanalysis.preprocessing.raw.extract import INDEX_COL

def fetch_gene_meta(protein_id_map, sep=';'):
    logger = get_logger(__name__)

    unique_ids = protein_id_map['Protein ID'].unique()
    mg_response = query_mygene_from_cache(unique_ids)

    df = pd.merge(protein_id_map, mg_response,
                  left_on='Protein ID', right_index=True)

    columns_to_merge = [x for x in MYGENE_FIELDS if x not in ['alias', 'ensembl.protein',
                                                              'interpro']]

    long_form_metadata_fields = ['interpro']

    new_df = []
    new_ix = []

    long_form_datasets = {}

    for protein_id_group, subdata in df.groupby(INDEX_COL):
        new_ix.append(protein_id_group)
        row = []

        other_names = set()
        for aliases in subdata['alias'].dropna():
            if isinstance(aliases, str):
                other_names.add(aliases)
            else:
                other_names.update(aliases)

        row.append(other_names)

        ensembl_protein_ids = set()
        for ensembl_data in subdata['ensembl'].dropna():
            if isinstance(ensembl_data, dict):
                # sometimes they return a list of dicts, sometimes just a dict...
                ensembl_data = [ensembl_data]

            for ensembl_subdata in ensembl_data:
                ensembl_subdata = ensembl_subdata['protein']
                if isinstance(ensembl_subdata, str):
                    # Sometimes they return just a string..
                    ensembl_protein_ids.add(ensembl_subdata)
                else:
                    ensembl_protein_ids.update(ensembl_subdata)

        row.append(ensembl_protein_ids)

        for col in columns_to_merge:
            series = subdata[col].dropna()
            if col == 'entrezgene':
                # Convert to int, then to str
                # this has to be done after NAs are removed though
                series = series.astype(int).astype(str)

            unique = series.unique()
            row.append(list(unique))

        new_df.append(row)

        for col in long_form_metadata_fields:
            if col not in long_form_datasets:
                long_form_datasets[col] = []

            dict_dataset = subdata[col].dropna()
            for dicts in dict_dataset.values:
                if isinstance(dicts, dict):
                    dicts = [dicts]

                for d in dicts:
                    d = pd.Series(d)
                    d[INDEX_COL] = protein_id_group
                    long_form_datasets[col].append(d)

    for col in long_form_datasets:
        long_form_datasets[col] = pd.DataFrame(long_form_datasets[col]).drop_duplicates()

        if col == 'interpro':
            long_form_datasets[col] = long_form_datasets[col].set_index([INDEX_COL, 'id'])
        else:
            raise NotImplementedError('Specify index column for {!r} please'.format(col))

    new_df = pd.DataFrame(new_df,
                          columns=['alias'] + ['ensembl_protein_id'] + columns_to_merge,
                          index=pd.Index(new_ix, name=INDEX_COL))
    # sort aliases & ensembl protein ids
    new_df['alias'] = new_df['alias'].apply(lambda x: sep.join(sorted(x)))
    new_df['ensembl_protein_id'] = new_df['ensembl_protein_id'].apply(lambda x: sep.join(sorted(x)))

    def _join(x):
        try:
            return sep.join(x)
        except TypeError:
            logger.debug('Error joining {!r}'.format(x))
            raise

    # Preserve order for other names
    for col in columns_to_merge:
        new_df[col] = new_df[col].apply(_join)

    return new_df, long_form_datasets, mg_response
