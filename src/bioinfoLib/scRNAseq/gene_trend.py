# Copyright 2025 Zhiyuan Yu (Heemskerk's lab, University of Michigan)

import numpy as np
from scipy.sparse import issparse


def rolling_mean(t, y, w, s, min_n=50):
    max_t = np.max(t)
    min_t = np.min(t)
    times = []
    means = []
    sds = []
    ts = np.arange(min_t, max_t + s, s)
    for i in range(len(ts)):
        up_bound = ts[i] + w / 2
        low_bound = ts[i] - w / 2
        y_sub = y[np.logical_and(t >= low_bound, t < up_bound)]
        if len(y_sub) < min_n:
            continue
        times.append(ts[i])
        means.append(np.mean(y_sub))
        sds.append(np.std(y_sub))
    return np.array(times), np.array(means), np.array(sds)


# TODO: make it easier to use
def compute_gene_trends(
    adata,
    genes,
    branches,
    branch_key="branch",
    covariate_key="dpt_pseudotime",
    lower_bound=2,
    upper_bound=99,
    window_size=0.03,
    step_size=0.001,
    min_n=50,
):
    trends = {}
    for lineage, branch_labels in branches.items():
        trends[lineage] = {}
        for gene in genes:
            if any(adata.var_names.str.contains(gene)):
                gene_name = adata.var_names[adata.var_names.str.contains(gene)][0]
                expr = adata[:, gene_name].X
                if issparse(expr):
                    expr = expr.toarray()
                else:
                    expr = np.array(expr)
                branch_mask = adata.obs[branch_key].isin(branch_labels)
                covariate = adata.obs[covariate_key][branch_mask]
                lower_bound = np.percentile(covariate, lower_bound)
                upper_bound = np.percentile(covariate, upper_bound)
                dpt = np.clip(covariate, lower_bound, upper_bound)
                times, mean, sd = rolling_mean(
                    dpt, expr[branch_mask], window_size, step_size, min_n
                )
                trends[lineage][gene] = {
                    "covariate": times,
                    "response_mean": mean,
                    "response_sd": sd,
                }
    return trends
