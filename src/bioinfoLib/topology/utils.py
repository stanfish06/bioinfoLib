import itertools
from itertools import product

import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.spatial.distance import cdist, hamming
from scipy.stats import rankdata
from sksparse.cholmod import cholesky

from bioinfoLib.R.utils import trajectory_distance


# TODO: grid search, adjust radius based on gaussian prior
# TODO: generate a set of representative loops through custom distance function
# TODO: pair up loops from two sets though frechet distance and do regression based matching
# TODO: if I generate a set of homologous loops, can I statistically test if a query loop is different or the same as this set?
def sp_ridge_regression_mod2(
    one_ridx_A,
    one_cidx_A,
    nrow_A,
    ncol_A,
    b,
    ridge_coef_a=0.1,
    ridge_coef_b=1,
    n_neighbors=1,
    do_approximation=False,
    n_search_cutoff=100,
    max_bits_flip=10,
    max_bits_comb=10,
    n_post_process=2,
):
    # simplify A by only considering the neighborhood of the loop
    if do_approximation:
        columns_kept = []
        one_idx_buff = np.where(b == 1)[0]
        for _ in range(n_neighbors):
            incident_triangles = np.unique(
                one_cidx_A[np.isin(one_ridx_A, one_idx_buff)]
            )
            columns_kept.extend(incident_triangles)
            one_idx_buff = np.unique(
                one_ridx_A[np.isin(one_cidx_A, incident_triangles)]
            )
        columns_kept = np.unique(columns_kept)
        mask = np.isin(one_cidx_A, columns_kept)
        one_ridx_A = one_ridx_A[mask]
        one_cidx_A = one_cidx_A[mask]
        if len(one_cidx_A) == 0:
            return None
        # reduce number of columns
        one_cidx_A = rankdata(one_cidx_A, method="min") - 1
        ncol_A = np.max(one_cidx_A) + 1

    # figure out number of ones in each row
    one_ridx_A_uniq, n_ones_per_row = np.unique(one_ridx_A, return_counts=True)
    one_ridx_A_uniq = one_ridx_A_uniq[n_ones_per_row > 1]
    n_ones_per_row = n_ones_per_row[n_ones_per_row > 1]
    nrow_A_new = len(one_ridx_A_uniq)
    # the new columns are used to perform 1 + 1 - 2 = 0
    data = np.concatenate(
        [
            np.ones(len(one_ridx_A)),
            np.ones(nrow_A_new) * (n_ones_per_row // 2) * -2,
        ]
    )
    A = csr_matrix(
        (
            data,
            (
                np.concatenate([one_ridx_A, one_ridx_A_uniq]),
                np.concatenate([one_cidx_A, np.arange(nrow_A_new) + ncol_A]),
            ),
        ),
        shape=(nrow_A, ncol_A + nrow_A_new),
    )
    ridge_coef_emp = ridge_coef_a * ridge_coef_b / (ridge_coef_a + ridge_coef_b)
    B = A.transpose().dot(A) + diags(np.repeat(ridge_coef_emp, ncol_A + nrow_A_new))
    factor_emp = cholesky(B)
    B = A.transpose().dot(A) + diags(np.repeat(ridge_coef_a, ncol_A + nrow_A_new))
    factor = cholesky(B)
    best_diff = np.Inf
    best_s = None

    x = factor_emp(
        A.transpose().dot(b) + np.repeat(ridge_coef_emp * 0.5, ncol_A + nrow_A_new)
    )
    c = (ridge_coef_a * x + ridge_coef_b * 0.5) / (ridge_coef_a + ridge_coef_b)
    min_cut = np.min(c)
    max_cut = np.max(c)
    for cut_i in np.linspace(min_cut, max_cut, n_search_cutoff):
        c_scale = (c - cut_i) / (2 * (np.max(np.abs(c - cut_i))) + 1e-6) + 0.5
        s = factor(A.transpose().dot(b) + ridge_coef_a * c_scale)
        min_cut_s = np.min(s)
        max_cut_s = np.max(s)
        for cut_j in np.linspace(min_cut_s, max_cut_s, n_search_cutoff):
            s_bin = s.copy()
            s_bin[s >= cut_j] = 1
            s_bin[s < cut_j] = 0
            pred = A[:, :ncol_A].dot(s_bin[:ncol_A]) % 2
            diff = hamming(pred, b)
            if diff < best_diff:
                best_diff = diff
                best_s = s_bin

    for _ in range(n_post_process):
        if best_diff == 0:
            break
        # try bit flip to improve solution
        grad_vec = np.zeros(ncol_A)
        for i in range(ncol_A):
            best_s_tmp = best_s.copy()
            best_s_tmp[i] = 1 - best_s_tmp[i]
            pred = A[:, :ncol_A].dot(best_s_tmp[:ncol_A]) % 2
            diff = hamming(pred, b)
            grad_vec[i] = diff - best_diff
        # find bits that improve the result
        good_bits = np.where(grad_vec <= 0)[0]
        if len(good_bits) > 0:
            best_bits = np.argsort(grad_vec[good_bits])[
                : min(max_bits_flip, len(good_bits))
            ]
            best_bits = good_bits[best_bits]
            best_s_k = best_s.copy()
            best_diff_k = best_diff
            for n_flips in range(1, min(max_bits_comb, len(best_bits)) + 1):
                for bit_comb in itertools.combinations(best_bits, n_flips):
                    best_s_tmp = best_s.copy()
                    for bit_idx in bit_comb:
                        best_s_tmp[bit_idx] = 1 - best_s_tmp[bit_idx]
                    pred = A[:, :ncol_A].dot(best_s_tmp[:ncol_A]) % 2
                    diff = hamming(pred, b)
                    if diff < best_diff:
                        best_s_k = best_s_tmp
                        best_diff_k = diff
                        if diff == 0:
                            break
            best_s = best_s_k.copy()
            best_diff = best_diff_k

    return (A[:, :ncol_A], best_s[:ncol_A])


def edge_idx_encode(i, j):
    if i > j:
        i_tmp = i
        i = j
        j = i_tmp
    return ((i + j) * (i + j + 1)) // 2 + j


def evaluate_match(
    nrow_A,
    ncol_A,
    max_frechet_dist,
    ridge_coef_a,
    ridge_coef_b,
    n_search,
    source_loops_edges,
    target_loops_edges,
    source_loops_coords,
    target_loops_coords,
    one_ridx_A,
    one_cidx_A,
    do_regression,
    do_approximation,
    n_neighbors,
    bd_column_birth_t,
    source_loop_death_t,
    similarity_func,
    type,
):
    n_sources = len(source_loops_edges)
    n_targets = len(target_loops_edges)
    if not do_regression:
        dists = []
        for i, j in product(range(n_sources), range(n_targets)):
            sloop_coords = source_loops_coords[i]
            tloop_coords = target_loops_coords[j]
            # find closest pair as the starting point of frechet distance
            dist_mat = cdist(sloop_coords, tloop_coords)
            p, q = np.unravel_index(np.argmin(dist_mat), dist_mat.shape)
            try:
                d1 = trajectory_distance(
                    np.roll(sloop_coords, -p, axis=0),
                    np.roll(tloop_coords, -q, axis=0),
                    type,
                    similarity_func,
                )[0]
                d2 = trajectory_distance(
                    np.roll(sloop_coords, -p, axis=0),
                    np.roll(tloop_coords[::-1], q + 1, axis=0),
                    type,
                    similarity_func,
                )[0]
                # not sure what negative value means
                if d1 < 0:
                    d1 = np.inf
                if d2 < 0:
                    d2 = np.inf
                dist = min(d1, d2)
                if dist < np.inf:
                    dists.append(dist)
                else:
                    dists.append(np.nan)
            except ValueError:
                dists.append(np.nan)
        dist = np.nanmean(dists)
        if dist > max_frechet_dist:
            return None
        else:
            return dist
    else:
        dists = []
        for i, j in product(range(n_sources), range(n_targets)):
            sloop_eidx = source_loops_edges[i]
            tloop_eidx = target_loops_edges[j]
            if len(sloop_eidx) == 0 or len(tloop_eidx) == 0:
                dists.append(np.nan)
                continue
            b1 = np.zeros(nrow_A)
            b1[sloop_eidx] = 1
            b2 = np.zeros(nrow_A)
            b2[tloop_eidx] = 1
            b = np.logical_xor(b1, b2).astype(int)
            # remove columns with birth t larger than the death t of the loop
            columns_kept = np.where(bd_column_birth_t <= source_loop_death_t)[0]
            if columns_kept.size == 0:
                return np.nan
            mask = ~np.isin(one_cidx_A, columns_kept)
            one_ridx_A = one_ridx_A[~mask]
            one_cidx_A = one_cidx_A[~mask]
            ncol_A = np.max(one_cidx_A) + 1
            res = sp_ridge_regression_mod2(
                one_ridx_A,
                one_cidx_A,
                nrow_A,
                ncol_A,
                b,
                ridge_coef_a,
                ridge_coef_b,
                do_approximation=do_approximation,
                n_neighbors=n_neighbors,
            )
            if res is None:
                dists.append(np.nan)
            else:
                A, s = res
                pred = np.round(A.dot(s)) % 2
                if np.sum(np.logical_or(pred, b)) == 0:
                    dist = np.nan
                else:
                    tp = np.sum(np.logical_and(pred == 1, b == 1))
                    fp = np.sum(np.logical_and(pred == 1, b == 0))
                    fn = np.sum(np.logical_and(pred == 0, b == 1))
                    if tp + fp + fn == 0:
                        dist = 0.0
                    else:
                        f1 = 2 * tp / (2 * tp + fp + fn)
                        dist = 1 - f1
                dists.append(dist)
        dist = np.nanmean(dists)
        return dist


def evaluate_match_worker(args):
    (
        nrow_A,
        ncol_A,
        max_frechet_dist,
        ridge_coef_a,
        ridge_coef_b,
        n_search,
        source_loop_edges,
        target_loop_edges,
        source_loop_coords,
        target_loop_coords,
        one_ridx_A,
        one_cidx_A,
        do_regression,
        do_approximation,
        n_neighbors,
        bd_column_birth_t,
        source_loop_birth_t,
    ) = args
    return evaluate_match(
        nrow_A,
        ncol_A,
        max_frechet_dist,
        ridge_coef_a,
        ridge_coef_b,
        n_search,
        source_loop_edges,
        target_loop_edges,
        source_loop_coords,
        target_loop_coords,
        one_ridx_A,
        one_cidx_A,
        do_regression,
        do_approximation,
        n_neighbors,
        bd_column_birth_t,
        source_loop_birth_t,
    )


def donut_2d_iso(r1, r2, n_points, noise, seed):
    np.random.seed(seed)
    rho = (r1 - r2) * np.sqrt(np.random.rand(n_points)) + r2
    theta = 2 * np.pi * np.random.rand(n_points)
    x = rho * np.cos(theta) + np.random.normal(0, noise, n_points)
    y = rho * np.sin(theta) + np.random.normal(0, noise, n_points)
    return (x, y)


def disk_2d_iso(r, n_points, noise, seed):
    np.random.seed(seed)
    rho = r * np.sqrt(np.random.rand(n_points))
    theta = 2 * np.pi * np.random.rand(n_points)
    x = rho * np.cos(theta) + np.random.normal(0, noise, n_points)
    y = rho * np.sin(theta) + np.random.normal(0, noise, n_points)
    return (x, y)


def disk_2d_two_holes_iso(r, r1, r2, c1_x, c1_y, c2_x, c2_y, n_points, noise, seed):
    np.random.seed(seed)
    rho = r * np.sqrt(np.random.rand(n_points * 2))
    theta = 2 * np.pi * np.random.rand(n_points * 2)
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)

    x1 = x - c1_x
    y1 = y - c1_y
    mask1 = np.sqrt(np.power(x1, 2) + np.power(y1, 2)) < r1

    x2 = x - c2_x
    y2 = y - c2_y
    mask2 = np.sqrt(np.power(x2, 2) + np.power(y2, 2)) < r2

    x = x[~np.logical_or(mask1, mask2)]
    y = y[~np.logical_or(mask1, mask2)]

    idx_keep = np.random.choice(np.arange(len(x)), size=n_points, replace=False)
    x = x[idx_keep] + np.random.normal(0, noise, n_points)
    y = y[idx_keep] + np.random.normal(0, noise, n_points)
    return (x, y)
