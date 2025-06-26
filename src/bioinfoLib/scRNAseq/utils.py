import numpy as np
import pandas as pd


def grouped_obs_mean(adata, group_key, use_raw=True, layer=None, genes=None):
    if genes:
        if use_raw:
            genes = adata.raw.var_names[
                adata.raw.var_names.str.contains("|".join(genes))
            ]
        else:
            genes = adata.var_names[adata.var_names.str.contains("|".join(genes))]
    else:
        if use_raw:
            genes = adata.raw.var_names
        else:
            genes = adata.var_names

    def getX(x, genes):
        if use_raw:
            return x.raw[:, genes].X
        else:
            if layer is not None:
                return x[:, genes].layers[layer]
            else:
                return x[:, genes].X

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((genes.__len__(), len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=genes,
    )

    for group, idx in grouped.indices.items():
        X = getX(adata, genes)
        out[group] = np.ravel(X[idx, :].mean(axis=0, dtype=np.float64))
    return out
