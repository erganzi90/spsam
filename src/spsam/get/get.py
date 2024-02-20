from anndata import AnnData


def check_indices(
        adata: AnnData,
        obs_key,
):
    if type(obs_key) == str:
        check_lt = [obs_key]
    elif type(obs_key) == list:
        check_lt = obs_key
    else:
        raise TypeError(f'check indices must be str or list, unknown type found:{obs_key}')
    not_found_lt = []
    for check_args in check_lt:
        if check_args not in adata.obs_keys():
            not_found_lt.append(check_args)
    if len(not_found_lt):
        raise KeyError(
            f'{obs_key} not in  adata.obs: {",".join(not_found_lt)}'
        )
