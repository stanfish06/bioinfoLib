from typing import Literal

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from fa2_modified import ForceAtlas2
from pydantic import validate_call
from umap.spectral import spectral_layout

from bioinfoLib.topology.utils import disk_2d_iso


def rotate_embedding_2d(adata, embedding_key: str, angle):
    rotation_matrix = np.array(
        [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    )
    adata.obsm[f"{embedding_key}_rotate"] = (
        adata.obsm[embedding_key] @ rotation_matrix.T
    )


def force_layout(
    adata,
    use_rep: str,
    initial_position=None,
    spectral_dim=10,
    spectral_components=(0, 1),
    distance_key="connectivities",
    num_iterations=10,
    scalingRatio=2.0,
    gravity=1.0,
    seed=0,
):
    assert distance_key in adata.obsp, "distance key not found"
    knn_dist_mat = adata.obsp[distance_key]
    data = adata.X if use_rep == "X" else adata.obsm[use_rep]
    # not sure why spectral embedding needs data. Check later
    # spectral initialization does not perform good, I wonder why
    if initial_position == "spectral":
        initial_position = spectral_layout(
            data=data, graph=knn_dist_mat, dim=spectral_dim, random_state=seed
        )[:, spectral_components]
    elif initial_position == "random_disk":
        initial_position = np.column_stack(
            disk_2d_iso(r=1, n_points=data.shape[0], noise=0, seed=seed)
        )
    else:
        initial_position = None
    # use spectral embedding as initial position
    fa = ForceAtlas2(verbose=True, scalingRatio=scalingRatio, gravity=gravity)
    adata.obsm["X_fl"] = np.array(
        fa.forceatlas2(G=knn_dist_mat, pos=initial_position, iterations=num_iterations)
    )


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
    mode: Literal["balanced", "pool"] = "balanced",
):
    """Differnetial expression

    Args:
        adata (AnnData): AnnData object
        label_key (str): Categorical key in obs
        targets (list[str]): Target groups, can have multiple groups
        reference (list[str], optional): Reference groups, can have multiple groups. If others, groups except target groups are included. Defaults to ["others"].
        mode (str, optional): Either balanced or pool. Defaults to "balanced".

    Raises:
        ValueError: _description_
        ValueError: _description_
    """
    group_column = np.full(adata.shape[0], "others", dtype=object)
    group_column_subgroups = np.full(adata.shape[0], "others", dtype=object)
    reference_group = "others"
    if reference:
        if "others" not in reference:
            reference_group = "|".join(f"({r})" for r in reference)
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
    # TODO: add separate mode that does rank gene for each target group and get shared markers
    if mode == "balanced":
        group_exp = grouped_obs_mean(adata, groupby + "_sub")
        target_means = group_exp.loc[:, targets].mean(axis=1)
        reference_means = group_exp.loc[:, reference].mean(axis=1)
        joint_means = (target_means + reference_means) / 2
        lfc_balanced = np.log2(
            (np.expm1(target_means) + 1e-9) / (np.expm1(reference_means) + 1e-9)
        )
        # TODO: potentially make DE balanced as well
        sc.tl.rank_genes_groups(
            adata,
            groups=[target_group],
            groupby=groupby,
            reference=reference_group,
            key_added=groupby,
        )
        de_df = sc.get.rank_genes_groups_df(adata, group=target_group, key=groupby)
        de_df["target_mean"] = target_means[de_df["names"]].__array__()
        de_df["reference_mean"] = reference_means[de_df["names"]].__array__()
        de_df["joint_mean"] = joint_means[de_df["names"]].__array__()
        de_df["lfc"] = lfc_balanced[de_df["names"]].__array__()
    elif mode == "pool":
        group_exp = grouped_obs_mean(adata, groupby)
        target_means = group_exp.loc[:, target_group]
        reference_means = group_exp.loc[:, reference_group]
        joint_means = (target_means + reference_means) / 2
        lfc = np.log2(
            (np.expm1(target_means) + 1e-9) / (np.expm1(reference_means) + 1e-9)
        )
        sc.tl.rank_genes_groups(
            adata,
            groups=[target_group],
            groupby=groupby,
            reference=reference_group,
            key_added=groupby,
        )
        de_df = sc.get.rank_genes_groups_df(adata, group=target_group, key=groupby)
        de_df["target_mean"] = target_means[de_df["names"]].__array__()
        de_df["reference_mean"] = reference_means[de_df["names"]].__array__()
        de_df["joint_mean"] = joint_means[de_df["names"]].__array__()
        de_df["lfc"] = lfc[de_df["names"]].__array__()
    else:
        raise ValueError("mean mode unsupported")

    adata.uns[groupby] = {"mode": mode, "result": de_df}
