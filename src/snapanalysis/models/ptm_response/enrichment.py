from snapanalysis.models.ptm_response.limma import limma_camera
from snapanalysis.preprocessing.protein_metadata import get_complexes_as_dict

import numpy as np

def limma_camera_complexes(matrix, design, weights,
                           limma_stats,
                           coef,
                           min_size=2, max_size=np.inf,
                           include_exclusive_subunits=True):
    gene_ids = matrix.index

    complexes = get_complexes_as_dict(gene_ids, min_size=min_size, max_size=max_size)

    if include_exclusive_subunits:
        complexes_exclusive_subunits = get_complexes_as_dict(gene_ids, min_size=min_size,
                                                             max_size=max_size,
                                                             exclusive_subunits_only=True)

        joint = {}

        for complex_, subunits in complexes.items():
            exclusive_subunits = complexes_exclusive_subunits.get(complex_, [])

            joint[complex_] = subunits
            if len(exclusive_subunits) > 0 and len(exclusive_subunits) != len(subunits):
                joint[f'{complex_} (exclusive subunits)'] = exclusive_subunits

        complexes = joint

    return limma_camera(matrix, design, weights,
                        limma_stats,
                        groups=complexes,
                        coef=coef,
                        group_name='Complex')
