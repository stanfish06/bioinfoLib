from typing import Literal, Optional

import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
from anndata import AnnData
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from pydantic import Field, validate_call


def volcano_plot_differential_expression(
    adata: AnnData,
    key: str,
    max_pval: float = 0.01,
    min_lfc: float = 1.0,
    min_exp: float = 0.0,
    x_mode: Literal["target", "joint"] = "joint",
    y_annot_up: str = "",
    y_annot_down: str = "",
    y_annot_xpos: float = 0.75,
    y_annot_ypos: float = 0.75,
    y_max: float = np.inf,
    legend_loc: str = "upper right",
    gene_domain: list[str] = None,
    ax: Optional[Axes] = None,
    **kwargs,
) -> Axes:
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (5, 5)))
        created_fig = True
    else:
        created_fig = False
    de_mode = adata.uns[key]["mode"]
    print(f"result was generated with the {de_mode} mode")
    result = adata.uns[key]["result"].copy()
    if x_mode == "joint":
        xval = "joint_mean"
    elif x_mode == "target":
        result = result.loc[result["lfc"] >= 0, :]
        xval = "target_mean"
    np.random.seed(42)
    result = result.iloc[np.random.permutation(result.shape[0]), :]
    result["scatter_color"] = "#000000"
    result.loc[result["pvals_adj"] < max_pval, "scatter_color"] = "#ADD8E6"
    result_sig = result.query(
        f"{xval} > {min_exp} & pvals_adj < {max_pval} & (lfc >= {min_lfc} | lfc <= -{min_lfc})"
    ).reset_index()

    ymax = max(np.max(result["lfc"]), np.max(-result["lfc"]))
    ymax = min(ymax, y_max)
    ymin = -0.5 if x_mode == "target" else -ymax * 1.1
    y_clip = result["lfc"].__array__()
    y_clip[y_clip > ymax] = ymax
    y_clip[y_clip < -ymax] = -ymax
    y_sig_clip = result_sig["lfc"].__array__()
    y_sig_clip[y_sig_clip > ymax] = ymax
    y_sig_clip[y_sig_clip < -ymax] = -ymax

    ax.scatter(
        result[xval],
        y_clip,
        c=result["scatter_color"],
        s=0.5,
        alpha=0.5,
    )
    ax.scatter(
        result_sig[xval],
        y_sig_clip,
        color="red",
        s=4,
        alpha=0.5,
    )
    texts = []

    for i, row in result_sig.iterrows():
        if gene_domain is not None:
            if row["names"].split()[0] not in gene_domain:
                continue
        texts.append(
            ax.text(
                x=row[xval],
                y=y_sig_clip[i],
                s=row["names"].split()[0],
                fontsize=8,
                alpha=1.0,
                ha="left",
                va="bottom",
                bbox=dict(
                    boxstyle="round,pad=0.1", facecolor="white", edgecolor="black"
                ),
            )
        )

    xlims = ax.get_xlim()
    xmin = -0.5
    xmax = xlims[1] * 1.1

    ax.text(
        x=xmax * y_annot_xpos,
        y=ymax * y_annot_ypos,
        s=y_annot_up,  # assuming gene names are in 'names' column
        fontsize=16,
        alpha=1.0,
        ha="left",  # horizontal alignment
        va="bottom",  # vertical alignment
    )
    ax.text(
        x=xmax * y_annot_xpos,
        y=-ymax * y_annot_ypos,
        s=y_annot_down,  # assuming gene names are in 'names' column
        fontsize=16,
        alpha=1.0,
        ha="left",  # horizontal alignment
        va="top",  # vertical alignment
    )
    ax.annotate(
        "",
        xy=(xmax * y_annot_xpos, ymax * y_annot_ypos * 0.95),
        xytext=(xmax * y_annot_xpos, ymax * y_annot_ypos * 0.6),
        arrowprops=dict(arrowstyle="->", color="black", lw=2),
    )
    ax.annotate(
        "",
        xy=(xmax * y_annot_xpos, -ymax * y_annot_ypos * 0.95),
        xytext=(xmax * y_annot_xpos, -ymax * y_annot_ypos * 0.6),
        arrowprops=dict(arrowstyle="->", color="black", lw=2),
    )
    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="black",
            markersize=6,
            label="not differentially expressed",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="#ADD8E6",
            markersize=6,
            label="differentially expressed",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="red",
            markersize=6,
            label="unique markers",
        ),
    ]
    ax.legend(
        handles=legend_elements,
        loc=legend_loc,
        frameon=True,
        fancybox=True,
        shadow=True,
        fontsize=8,
    )
    ax.set_xticks(np.round(np.linspace(0, xmax, 5)))
    yticks = np.round(np.linspace(0, ymax, 4), 1)
    yticks = np.sort(np.unique(np.concatenate([-yticks, yticks])))
    yticks_lab = yticks.astype(str)
    yticks_lab[0] = f"<= -{ymax}"
    yticks_lab[-1] = f">= {ymax}"
    ax.set_yticks(yticks, labels=yticks_lab)
    ax.set_ylabel(r"Log$_{2}$ Fold Change")
    ax.set_xlabel("Mean Expression")

    ax.set_ylim([ymin, ymax * 1.1])
    ax.set_xlim([xmin, xmax])

    adjust_text(
        texts,
        arrowprops=dict(
            arrowstyle="-",
            color="black",
            lw=1,
            alpha=0.5,
            shrinkA=10,
            shrinkB=0,
        ),
        expand=(1.25, 1.25),
        avoid_text=True,
        avoid_point=True,
        prevent_crossings=True,
        min_arrow_len=1,
        max_move=[15, 15],
    )
    return ax


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
