from pandas import DataFrame
from anndata import AnnData
import numpy as np
from ..get import check_indices


def get_cell_abundance(
        adata: AnnData,
        cluster_df: DataFrame,
        cell_type_lt,
):
    check_indices(adata, cell_type_lt)  # check if cell type exists in adata.obs_keys()
    for cell_type in cell_type_lt:
        abundance_value_lt = []
        for idx, info in cluster_df.iterrows():
            row = info['row']
            col = info['col']
            cycle = info['cycle']
            if np.isnan(cycle):
                cell_abundance = np.nan
            else:
                cell_abundance = \
                    adata.obs.loc[(adata.obs['array_row'] == row) & (adata.obs['array_col'] == col), cell_type][0]
            abundance_value_lt.append(cell_abundance)
        cluster_df[cell_type] = abundance_value_lt
    # return cluster_df
