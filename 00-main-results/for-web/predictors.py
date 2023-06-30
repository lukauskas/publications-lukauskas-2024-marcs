import pandas as pd
from snapanalysis.preprocessing.pulldown_metadata import OUTPUT_FILE as META_FILE

import json
import os

from snapanalysis.config import OUTPUT_DIRECTORY as SNAPANALYSIS_OUTPUT_DIRECTORY

OUTPUT_DIR = os.path.join(SNAPANALYSIS_OUTPUT_DIRECTORY, 'web/precompiled/')
OUTPUT_FILE = os.path.join(OUTPUT_DIR, 'pull_downs.js')

def _generate_json(names_and_types, predictors_web):
    """
    Creates a cached JSON to feed the pulldown list on interactive interface

    """

    payload = {}

    for id_, row in names_and_types.iterrows():

        d = {
            'key': id_,
            'name': row['Pull-Down name'],
            'type': row['Type'],
            'subtype': row['Subtype']
        }

        predictors_d = {}
        predictors = predictors_web[predictors_web['Pull-Down ID'] == id_]
        for __, row in predictors.iterrows():
             predictors_d[row['predictor']] = row['predictor_value']

        d['predictors'] = predictors_d

        payload[id_] = d

    return json.dumps(payload)

def main():

    names_and_types = pd.read_hdf(META_FILE,
                                 '/meta/names_and_types')
    predictors_web = pd.read_hdf(META_FILE, '/meta/predictors_web')

    pulldowns_json = _generate_json(names_and_types, predictors_web)

    if not os.path.isdir(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    TEMPLATE = """
    export const PULL_DOWNS = {pulldowns_json};
    """

    with open(OUTPUT_FILE, 'w') as f:
        f.write(TEMPLATE.format(pulldowns_json=pulldowns_json))


if __name__ == '__main__':
    main()

