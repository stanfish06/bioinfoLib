import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from pydantic import validate_call


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


@validate_call(config={"arbitrary_types_allowed": True})
def pairwise_differential_expression(
    adata: AnnData,
    label_key: str,
    targets: list[str],
    reference: list[str] = ["others"],
    max_pval: float = 0.05,
    min_lfc: float = 1.0,
    min_exp_target: float = 0.0,
    mode: str = "balanced",
):
    group_column = np.full(adata.shape[0], "others", dtype=object)
    group_column_subgroups = np.full(adata.shape[0], "others", dtype=object)
    reference_group = "others"
    if reference:
        if "others" not in reference:
            refernce_group = "|".join(f"({r})" for r in reference)
            for r in reference:
                group_column[adata.obs[label_key] == r] = reference_group
                group_column_subgroups[adata.obs[label_key] == r] = r
    if targets:
        target_group = "|".join(f"({t})" for t in targets)
        for t in targets:
            group_column[adata.obs[label_key] == t] = target_group
            group_column_subgroups[adata.obs[label_key] == t] = t
    else:
        raise ValueError("targets must contain at least one group")

    if reference_group == "others":
        group_column_subgroups[
            ~adata.obs[label_key].str.contains("|".join(targets))
        ] = adata.obs[label_key][~adata.obs[label_key].str.contains("|".join(targets))]
        reference = np.unique(
            adata.obs[label_key][~adata.obs[label_key].str.contains("|".join(targets))]
        )

    groupby = f"de_{target_group}_to_{reference_group}"
    adata.obs[groupby] = group_column
    adata.obs[groupby + "_sub"] = group_column_subgroups
    if mode == "balanced":
        group_exp = grouped_obs_mean(adata, groupby + "_sub")
        target_means = group_exp.loc[:, targets].mean(axis=1)
        reference_means = group_exp.loc[:, reference].mean(axis=1)
        lfc_balanced = target_means - reference_means
        # TODO: potentially make DE balanced as well
        sc.tl.rank_genes_groups(
            adata,
            groups=[target_group],
            groupby=groupby,
            reference=reference_group,
            key_added=groupby,
        )
        de_df = sc.get.rank_genes_groups_df(adata, group=target_group, key=groupby)
        de_df[f"{target_group}_balanced_mean"] = target_means[
            de_df["names"]
        ].__array__()
        de_df[f"{reference_group}_balanced_mean"] = reference_means[
            de_df["names"]
        ].__array__()
        de_df["balanced_lfc"] = lfc_balanced[de_df["names"]].__array__()
    elif mode == "pool":
        group_exp = grouped_obs_mean(adata, groupby)
    else:
        raise ValueError("mean mode unsupported")

    adata.uns[groupby] = de_df
