from pandas import DataFrame
import pandas as pd
from sklearn.preprocessing import MinMaxScaler


def minmax_scaler(
        cluster_df: DataFrame,
        cell_type_lt,
):
    # Map the cell abundance values to the range of [0,1].
    cluster_df1 = cluster_df.loc[:, 'row':'cycle']
    cluster_df2 = cluster_df.loc[:, cell_type_lt[0]:cell_type_lt[-1]]

    # initialize MinMaxScaler
    scaler = MinMaxScaler()
    scaled_df = pd.DataFrame(scaler.fit_transform(cluster_df2), columns=cluster_df2.columns)
    # scaled_pct_df = scaled_df.div(scaled_df.sum(axis=1), axis=0)
    cluster_scale_df = pd.concat([cluster_df1, scaled_df], axis=1)

    # calculate mean cell abundance, column "mean_abundance" added into dataframe
    cluster_scale_df['mean_abundance'] = cluster_scale_df[cell_type_lt].mean(axis=1)

    # score divided into 10 equal parts, column "score_bin10" added into dataframe
    cluster_scale_df['score_bin10'] = pd.cut(cluster_scale_df['score'], bins=10,
                                             labels=[f'{i}0%' for i in range(1, 11)])
    return cluster_scale_df
