"""Calculate scores based on the expression of gene lists."""
from typing import Sequence, Optional, Union
from anndata import AnnData
from scipy.sparse import issparse
import numpy as np
import pandas as pd
from numpy import random
import logging


def _check_use_raw(adata: AnnData, use_raw: Union[None, bool]) -> bool:
    """
    Normalize checking `use_raw`.

    My intentention here is to also provide a single place to throw a deprecation warning from in future.
    """
    if use_raw is not None:
        return use_raw
    else:
        if adata.raw is not None:
            return True
        else:
            return False


def _sparse_nanmean(X, axis):
    """
    np.nanmean equivalent for sparse matrices
    """
    if not issparse(X):
        raise TypeError("X must be a sparse matrix")

    # count the number of nan elements per row/column (dep. on axis)
    Z = X.copy()
    Z.data = np.isnan(Z.data)
    Z.eliminate_zeros()
    n_elements = Z.shape[axis] - Z.sum(axis)

    # set the nans to 0, so that a normal .sum() works
    Y = X.copy()
    Y.data[np.isnan(Y.data)] = 0
    Y.eliminate_zeros()

    # the average
    s = Y.sum(axis, dtype='float64')  # float64 for score_genes function compatibility)
    m = s / n_elements

    return m


def score_genes(
        adata: AnnData,
        gene_list: Sequence[str],
        ctrl_size: int = 50,
        gene_pool: Optional[Sequence[str]] = None,
        n_bins: int = 25,
        score_name: str = 'score',
        random_state: Union[None, int, random.RandomState] = 0,
        copy: bool = False,
        use_raw: Optional[bool] = None,
) -> Optional[AnnData]:
    """\
    Score a set of genes [Satija15]_.

    The score is the average expression of a set of genes subtracted with the
    average expression of a reference set of genes. The reference set is
    randomly sampled from the `gene_pool` for each binned expression value.

    This reproduces the approach in Seurat [Satija15]_ and has been implemented
    for Scanpy by Davide Cittaro.

    Parameters
    ----------
    adata
        The annotated data matrix.
    gene_list
        The list of gene names used for score calculation.
    ctrl_size
        Number of reference genes to be sampled from each bin. If `len(gene_list)` is not too
        low, you can set `ctrl_size=len(gene_list)`.
    gene_pool
        Genes for sampling the reference set. Default is all genes.
    n_bins
        Number of expression level bins for sampling.
    score_name
        Name of the field to be added in `.obs`.
    random_state
        The random seed for sampling.
    copy
        Copy `adata` or modify it inplace.
    use_raw
        Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present.

        .. versionchanged:: 1.4.5
           Default value changed from `False` to `None`.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with an additional field
    `score_name`.
    """
    print(f'computing score {score_name!r}')
    adata = adata.copy() if copy else adata
    use_raw = _check_use_raw(adata, use_raw)

    if random_state is not None:
        np.random.seed(random_state)

    gene_list_in_var = []
    var_names = adata.raw.var_names if use_raw else adata.var_names
    genes_to_ignore = []
    for gene in gene_list:
        if gene in var_names:
            gene_list_in_var.append(gene)
        else:
            genes_to_ignore.append(gene)
    if len(genes_to_ignore) > 0:
        logging.warning(f'genes are not in var_names and ignored: {genes_to_ignore}')
    gene_list = set(gene_list_in_var[:])

    if len(gene_list) == 0:
        raise ValueError("No valid genes were passed for scoring.")

    if gene_pool is None:
        gene_pool = list(var_names)
    else:
        gene_pool = [x for x in gene_pool if x in var_names]
    if not gene_pool:
        raise ValueError("No valid genes were passed for reference set.")

    # Trying here to match the Seurat approach in scoring cells.
    # Basically we need to compare genes against random genes in a matched
    # interval of expression.

    _adata = adata.raw if use_raw else adata
    _adata_subset = (
        _adata[:, gene_pool] if len(gene_pool) < len(_adata.var_names) else _adata
    )
    if issparse(_adata_subset.X):
        obs_avg = pd.Series(
            np.array(_sparse_nanmean(_adata_subset.X, axis=0)).flatten(),
            index=gene_pool,
        )  # average expression of genes
    else:
        obs_avg = pd.Series(
            np.nanmean(_adata_subset.X, axis=0), index=gene_pool
        )  # average expression of genes

    obs_avg = obs_avg[
        np.isfinite(obs_avg)
    ]  # Sometimes (and I don't know how) missing data may be there, with nansfor

    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method='min') // n_items
    control_genes = set()

    # now pick `ctrl_size` genes from every cut
    for cut in np.unique(obs_cut.loc[list(gene_list)]):
        r_genes = np.array(obs_cut[obs_cut == cut].index)
        np.random.shuffle(r_genes)
        # uses full r_genes if ctrl_size > len(r_genes)
        control_genes.update(set(r_genes[:ctrl_size]))

    # To index, we need a list – indexing implies an order.
    control_genes = list(control_genes - gene_list)
    gene_list = list(gene_list)

    X_list = _adata[:, gene_list].X
    if issparse(X_list):
        X_list = np.array(_sparse_nanmean(X_list, axis=1)).flatten()
    else:
        X_list = np.nanmean(X_list, axis=1, dtype='float64')

    X_control = _adata[:, control_genes].X
    if issparse(X_control):
        X_control = np.array(_sparse_nanmean(X_control, axis=1)).flatten()
    else:
        X_control = np.nanmean(X_control, axis=1, dtype='float64')

    score = X_list - X_control

    adata.obs[score_name] = pd.Series(
        np.array(score).ravel(), index=adata.obs_names, dtype='float64'
    )
    print('finished\n'
          f'{score_name}, score of gene set (adata.obs).\n'
          f'{len(control_genes)} total control genes are used.'
          )
    return adata if copy else None
