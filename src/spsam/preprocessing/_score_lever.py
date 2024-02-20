import numpy as np
from anndata import AnnData
from ..get import check_indices


def score_lever(
        adata: AnnData,
        score_key,
        lever1: float = 0.95,
        lever2: float = 0.5
):
    check_indices(adata, score_key)
    score_type_lt = []
    med = np.median(adata.obs[score_key])
    mean = np.mean(adata.obs[score_key])
    std = np.std(adata.obs[score_key])
    print(f'{score_key} median:{med}; mean:{mean}; std:{std}')
    l1 = adata.obs[score_key].quantile(lever1)
    l2 = adata.obs[score_key].quantile(lever2)
    print(f'{score_key} lever1 value:{l1}; lever2 value:{l2}')
    lever1, lever2 = l1, l2
    for index, row in adata.obs.iterrows():
        score = row[score_key]
        if score >= lever1:
            score_type = 'lever1'
        elif score >= lever2:
            score_type = 'lever2'
        else:
            score_type = 'background'
        score_type_lt.append(score_type)
    adata.obs['score_type'] = score_type_lt
    print('score lever definition finished, score_type key is added to adata.obs, use adata.obs_keys() to check')
    