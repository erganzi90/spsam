"""use depth-first algorithm(DFS) to cluster visium spot"""

import pandas as pd
from anndata import AnnData


def _get_spot_val(row, col):
    # give one spot location, return six spot information around it
    row1, col1 = row, col + 2
    row2, col2 = row, col - 2
    row3, col3 = row - 1, col + 1
    row4, col4 = row - 1, col - 1
    row5, col5 = row + 1, col + 1
    row6, col6 = row + 1, col - 1
    surround_spot_lt = []
    for i in range(1, 7):
        this_row, this_col = eval(f'row{i}'), eval(f'col{i}')
        type_series = adata_obs_df[(adata_obs_df['array_row'] == this_row) & (adata_obs_df['array_col'] == this_col)][
            'score_type']
        score_type = type_series[0] if not type_series.empty else None
        if score_type:
            surround_spot_lt.append((this_row, this_col, score_type))
    return surround_spot_lt


def _find_connected_points(x, y, color, connected_color_lt, used_lever1_lt, connected_points):
    if (x, y, color) not in used_lever1_lt and color in connected_color_lt:
        used_lever1_lt.append((x, y, color))
        connected_points.append((x, y, color))
        # print(x,y, connected_points)
        surround_spot_lt = _get_spot_val(x, y)
        for spot in surround_spot_lt:
            new_x, new_y, new_color = spot
            _find_connected_points(new_x, new_y, new_color, connected_color_lt, used_lever1_lt, connected_points)


def _connect_color_spot(connected_color_lt):
    # connect spot according to the given lever [lever1, lever2] through DFS algorithm
    connected_groups = []
    used_lever1_lt = []
    for index, row in adata_obs_df[adata_obs_df['score_type'] == 'lever1'].iterrows():  # 遍历每一个红点
        spot_row = row['array_row']
        spot_col = row['array_col']
        color = 'lever1'
        if (spot_row, spot_col, color) in used_lever1_lt:
            continue
        connected_points = []
        _find_connected_points(spot_row, spot_col, color, connected_color_lt, used_lever1_lt, connected_points)
        connected_groups.append(connected_points)
    return connected_groups


def _get_four_direction(spot_cluster):
    left = min([spot[0] for spot in spot_cluster])
    right = max([spot[0] for spot in spot_cluster])
    lower = min([spot[1] for spot in spot_cluster])
    upper = max([spot[1] for spot in spot_cluster])
    return left, right, lower, upper


def _get_border(red_cluster, red_row_lt, red_col_lt):
    # give lever1 cluster split-coordinate according other lever1 cluster around it
    core_spot_left, core_spot_right, core_spot_lower, core_spot_upper = _get_four_direction(red_cluster)
    border_left, border_lower, border_right, border_upper = 0, 0, 1000, 1000  # 定义好初始边界值
    if core_spot_left != red_row_lt[0]:
        border_left = red_row_lt[red_row_lt.index(core_spot_left) - 1]  # 边界就是它前一个点
    if core_spot_right != red_row_lt[-1]:
        border_right = red_row_lt[
            len(red_row_lt) - red_row_lt[::-1].index(core_spot_right) - 1 + 1]  # 找元素最后出现的index号后一个元素
    if core_spot_lower != red_col_lt[0]:
        border_lower = red_col_lt[red_col_lt.index(core_spot_lower) - 1]
    if core_spot_upper != red_col_lt[-1]:
        border_upper = red_col_lt[
            len(red_col_lt) - red_col_lt[::-1].index(core_spot_upper) - 1 + 1]  # 找元素最后出现的index号后一个元素
    return border_left, border_lower, border_right, border_upper


def find_lever_core(
        adata: AnnData,
        score_key,
):
    global adata_obs_df
    adata_obs_df = adata.obs.copy()
    global input_key
    input_key = score_key
    only_red_groups_lt = _connect_color_spot(['lever1'])  # connect lever1 spot
    red_yellow_groups_lt = _connect_color_spot(['lever1', 'lever2'])  # connect lever1 and lever2 spot
    only_red_dt = {}
    for red_group in only_red_groups_lt:
        for spot in red_group:
            only_red_dt[spot] = sorted(red_group)
    cluster_df = pd.DataFrame(columns=['row', 'col', 'score', 'type', 'cluster_idx', 'class'])
    independent_idx = 0
    adjacent_idx = 0
    for spot_cluster in red_yellow_groups_lt:
        red_cluster_lt = []
        for spot in spot_cluster:
            row, col, tpe = spot
            if tpe == 'lever1' and spot not in sum(red_cluster_lt, []):
                connected_lever1 = only_red_dt[spot]
                red_cluster_lt.append(connected_lever1)
        if len(red_cluster_lt) == 1:  # independent type
            for spot in spot_cluster:
                this_row, this_col, this_tpe = spot
                score_series = \
                adata_obs_df[(adata_obs_df['array_row'] == this_row) & (adata_obs_df['array_col'] == this_col)][
                    input_key]
                hypoxia_score = score_series[0] if not score_series.empty else None
                row_info1 = [this_row, this_col, hypoxia_score, this_tpe, independent_idx, 'independent']
                cluster_df.loc[len(cluster_df)] = row_info1
            independent_idx += 1
        elif len(red_cluster_lt) > 1:  # adjacent type
            red_row_lt = sorted([spot[0] for spot in sum(red_cluster_lt, [])])
            red_col_lt = sorted([spot[1] for spot in sum(red_cluster_lt, [])])
            for red_cluster in red_cluster_lt:
                border_left, border_lower, border_right, border_upper = _get_border(red_cluster, red_row_lt, red_col_lt)
                sub_cluster_lt = []
                for spot in spot_cluster:
                    if (border_left < spot[0] < border_right) and (border_lower < spot[1] < border_upper):
                        sub_cluster_lt.append(spot)
                for sub_spot in sub_cluster_lt:
                    this_row, this_col, this_tpe = sub_spot
                    score_series = \
                    adata_obs_df[(adata_obs_df['array_row'] == this_row) & (adata_obs_df['array_col'] == this_col)][
                        input_key]
                    hypoxia_score = score_series[0] if not score_series.empty else None
                    row_info2 = [this_row, this_col, hypoxia_score, this_tpe, adjacent_idx, 'adjacent']
                    cluster_df.loc[len(cluster_df)] = row_info2
                adjacent_idx += 1
    return cluster_df
