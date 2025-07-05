from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from adjustText import adjust_text
from anndata import AnnData
from matplotlib.axes import Axes
from pydantic import Field, validate_call


def volcano_plot_differential_expression(adata, group: str):
    df = sc.get.rank_genes_groups_df(adata, group=group)


@validate_call(config={"arbitrary_types_allowed": True})
def embedding_label_repl(
    adata: AnnData,
    basis: tuple,
    groupby: str,
    exclude: list[str],
    ax: Optional[Axes] = None,
    components: tuple = (0, 1),
    adjust_kwargs: dict = Field(default_factory=dict),
    text_kwargs: dict = Field(default_factory=dict),
    color_by_group: bool = False,
    palette: dict = Field(default_factory=dict),
    text_path_effect: dict = Field(default_factory=dict),
):
    if adjust_kwargs is None:
        adjust_kwargs = {"text_from_points": False}
    if text_kwargs is None:
        text_kwargs = {}

    medians = {}

    for g, g_idx in adata.obs.groupby(groupby).groups.items():
        if g in exclude:
            continue
        medians[g] = np.median(adata[g_idx].obsm[basis][:, components], axis=0)

    # Fill the text colors dictionary
    text_colors = {group: None for group in adata.obs[groupby].cat.categories}

    if color_by_group:
        if palette is not None:
            for i, group in enumerate(adata.obs[groupby].cat.categories):
                if group in exclude:
                    continue
                text_colors[group] = palette[group]
        elif groupby + "_colors" in adata.uns:
            for i, group in enumerate(adata.obs[groupby].cat.categories):
                if group in exclude:
                    continue
                text_colors[group] = adata.uns[groupby + "_colors"][i]

    if ax is None:
        texts = [
            plt.text(x=x, y=y, s=k, color=text_colors[k], **text_kwargs)
            for k, (x, y) in medians.items()
        ]
    else:
        texts = [
            ax.text(x=x, y=y, s=k, color=text_colors[k], **text_kwargs)
            for k, (x, y) in medians.items()
        ]

    if text_path_effect is not None:
        for text in texts:
            text.set_path_effects(text_path_effect)

    adjust_text(texts, **adjust_kwargs)
