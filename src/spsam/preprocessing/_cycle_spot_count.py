"""count spot number in each cycle around lever1"""

from pandas import DataFrame
import numpy as np
from collections import defaultdict


def cycle_spot_count(
        cluster_df: DataFrame,
):

    def _get_spot_distance(orange_df, red_df):
        for idx, orange_line in orange_df.iterrows():
            tmp = []
            orange_row = orange_line['row']
            orange_col = orange_line['col']
            for idx2, red_line in red_df.iterrows():
                red_row = red_line['row']
                red_col = red_line['col']
                tmp.append([abs(orange_row - red_row), abs(orange_col - red_col)])
            tmp.sort(key=lambda x: (x[0] + x[1], x[0], x[1]))
            picked_row, picked_col = tmp[0]
            cycle = round(picked_col / 2 + 0.1) if round(picked_col / 2 + 0.1) > picked_row else picked_row
            cluster_df.loc[(cluster_df['row'] == orange_row) & (cluster_df['col'] == orange_col), 'cycle'] = cycle

    # add "cycle" column to dataframe: all lever1 spot belong to cycle 0, lever2 belong to cycle 1-N
    cluster_df['cycle'] = np.nan
    cluster_df.loc[cluster_df['type'] == 'lever1', 'cycle'] = 0
    for cls in ['independent', 'adjacent']:
        df = cluster_df[cluster_df['class'] == cls]
        idx_lt = [i for i, j in df['cluster_idx'].value_counts().items() if j > 3]
        for idx in idx_lt:
            red_df = df[(df['cluster_idx'] == idx) & (df['type'] == 'lever1')]
            orange_df = df[(df['cluster_idx'] == idx) & (df['type'] == 'lever2')]
            _get_spot_distance(orange_df, red_df)
    # return cluster_df
