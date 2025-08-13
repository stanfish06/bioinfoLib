from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from pandas import DataFrame


def volcano_plot_differential_expression(
    df: DataFrame,
    LFC_key: str,
    xval: str,
    annot_key: str,
    max_pval: float = 0.01,
    min_lfc: float = 1.0,
    min_exp: float = 0.0,
    y_annot_up: str = "",
    y_annot_down: str = "",
    y_annot_xpos: float = 0.75,
    y_annot_ypos: float = 0.75,
    y_max: float = np.inf,
    legend_loc: str = "upper right",
    legend_fontsize: int = 8,
    gene_domain: list[str] = None,
    ax: Optional[Axes] = None,
    **kwargs,
) -> Axes:
    if ax is None:
        fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (5, 5)))
        created_fig = True
    else:
        created_fig = False
    np.random.seed(42)
    df = df.iloc[np.random.permutation(df.shape[0]), :]
    df["scatter_color"] = "#000000"
    df_sig = df.query(
        f"{xval} > {min_exp} & ({LFC_key} >= {min_lfc} | {LFC_key} <= -{min_lfc})"
    ).reset_index()

    ymax = max(np.max(df[LFC_key]), np.max(-df[LFC_key]))
    ymax = max(ymax, y_max)
    ymin = -ymax * 1.1
    y_clip = df[LFC_key].__array__()
    y_clip[y_clip > ymax] = ymax
    y_clip[y_clip < -ymax] = -ymax
    y_sig_clip = df_sig[LFC_key].__array__()
    y_sig_clip[y_sig_clip > ymax] = ymax
    y_sig_clip[y_sig_clip < -ymax] = -ymax

    ax.scatter(
        df[xval],
        y_clip,
        c=df["scatter_color"],
        s=0.5,
        alpha=0.5,
    )
    ax.scatter(
        df_sig[xval],
        y_sig_clip,
        color="red",
        s=4,
        alpha=0.5,
    )
    texts = []

    for i, row in df_sig.iterrows():
        if gene_domain is not None:
            if row[annot_key].split()[0] not in gene_domain:
                continue
        texts.append(
            ax.text(
                x=row[xval],
                y=y_sig_clip[i],
                s=row[annot_key].split()[0],
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
        fontsize=legend_fontsize,
    )
    ax.set_xticks(np.round(np.linspace(0, xmax, 5)))
    yticks = np.round(np.linspace(0, ymax, 4), 1)
    yticks = np.sort(np.unique(np.concatenate([-yticks, yticks])))
    yticks_lab = yticks.astype(str)
    yticks_lab[0] = f"<= -{np.round(ymax, 1)}"
    yticks_lab[-1] = f">= {np.round(ymax, 1)}"
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
