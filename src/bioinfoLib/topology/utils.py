import itertools
import warnings
from itertools import product

import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.spatial.distance import cdist, hamming
from sksparse.cholmod import cholesky

from bioinfoLib.linearAlgebra.gauss_mod2_m4ri import mod2Solve_m4ri_py
from bioinfoLib.R.utils import trajectory_distance

warnings.filterwarnings("ignore", category=RuntimeWarning)


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
    n_search_cutoff=100,
    max_bits_flip=10,
    max_bits_comb=10,
    n_post_process=2,
):
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


def edge_idx_encode(i, j, n_vertices, self_edge=True):
    if i > j:
        i, j = j, i
    # // make sure it returns integer
    if self_edge:
        return i * n_vertices + j
    else:
        return i * n_vertices - i * (i + 1) // 2 + (j - i - 1)


def trig_idx_encode(i, j, k, n_vertices, self_edge=True):
    i, j, k = sorted([i, j, k])
    # // make sure it returns integer
    if self_edge:
        return i * np.power(n_vertices, 2) + j * n_vertices + k
    else:
        return None


def compute_geometric_similarity(
    source_loops_coords,
    target_loops_coords,
    max_frechet_dist,
    similarity_func,
    similarity_type,
):
    n_sources = len(source_loops_coords)
    n_targets = len(target_loops_coords)
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
                similarity_type,
                similarity_func,
            )[0]
            d2 = trajectory_distance(
                np.roll(sloop_coords, -p, axis=0),
                np.roll(tloop_coords[::-1], q + 1, axis=0),
                similarity_type,
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


def compute_homological_equivalence(
    source_loops_edges,
    target_loops_edges,
    one_ridx_A,
    one_cidx_A,
    nrow_A,
    ncol_A,
    ridge_coef_a,
    ridge_coef_b,
    bd_column_birth_t,
    source_loop_death_t,
    regression_mode,
    do_approximation=False,  # this does not perform good
    n_neighbors=1,
    do_downsample=True,  # this is better than nn approximation
    fraction_downsample=0.0,
):
    n_sources = len(source_loops_edges)
    n_targets = len(target_loops_edges)
    dists = []
    # remove columns with birth t larger than the death t of the loop
    columns_use = np.where(bd_column_birth_t <= source_loop_death_t)[0]
    if columns_use.size == 0:
        return np.nan
    mask = np.isin(one_cidx_A, columns_use)
    one_ridx_A_local = one_ridx_A[mask]
    one_cidx_A_local = one_cidx_A[mask]
    _, one_cidx_A_local = np.unique(one_cidx_A_local, return_inverse=True)
    ncol_A_local = np.max(one_cidx_A_local) + 1
    for i, j in product(range(n_sources), range(n_targets)):
        bd_column_birth_t_sub = bd_column_birth_t[columns_use].copy()
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
        if do_approximation:
            columns_kept = []
            one_idx_buff = np.where(b == 1)[0]
            for _ in range(n_neighbors):
                incident_triangles = np.unique(
                    one_cidx_A_local[np.isin(one_ridx_A_local, one_idx_buff)]
                )
                columns_kept.extend(incident_triangles)
                one_idx_buff = np.unique(
                    one_ridx_A_local[np.isin(one_cidx_A_local, incident_triangles)]
                )
            columns_kept = np.unique(columns_kept)
            if len(columns_kept) > 0:
                bd_column_birth_t_sub = bd_column_birth_t_sub[columns_kept]
                mask = np.isin(one_cidx_A_local, columns_kept)
                one_ridx_A_local = one_ridx_A_local[mask]
                one_cidx_A_local = one_cidx_A_local[mask]
            else:
                dists.append(np.nan)
            # reduce number of columns
            _, one_cidx_A_local = np.unique(one_cidx_A_local, return_inverse=True)
            ncol_A_local = np.max(one_cidx_A_local) + 1
        if ncol_A_local > nrow_A:
            mask = np.argsort(bd_column_birth_t_sub)[:nrow_A]
            bd_column_birth_t_sub = bd_column_birth_t_sub[mask]
            mask = np.isin(one_cidx_A_local, mask)
            one_ridx_A_local = one_ridx_A_local[mask]
            one_cidx_A_local = one_cidx_A_local[mask]
            _, one_cidx_A_local = np.unique(one_cidx_A_local, return_inverse=True)
            ncol_A_local = np.max(one_cidx_A_local) + 1
        if regression_mode == "exact":
            # if number of columns is greater than number of rows, remove triangles with large birth t
            # large triangles are less likely be an important component between two equivalent loops
            # print(f"A: {nrow_A} x {ncol_A_local}")
            A = np.zeros([nrow_A, ncol_A_local])
            A[one_ridx_A_local, one_cidx_A_local] = 1
            # if whole row and b are zeros, then ignore that row
            zero_rows = np.logical_and(np.all(A == 0, axis=1), b == 0)
            nrow_A_new = nrow_A
            nrow_A_new = nrow_A_new - np.sum(zero_rows)
            # fail immediately if b is one but A is all zeros
            unsolvable_rows = np.logical_and(np.all(A == 0, axis=1), b == 1)
            if np.any(unsolvable_rows):
                res = 0
            else:
                if ncol_A_local > nrow_A_new:
                    mask = np.argsort(bd_column_birth_t_sub)[:nrow_A_new]
                    if do_downsample:
                        mask = mask[
                            : int(np.ceil(nrow_A_new * (1 - fraction_downsample)))
                        ]
                    if nrow_A_new == 0:
                        res = 1
                    else:
                        res = mod2Solve_m4ri_py(
                            A[~zero_rows, :][:, mask], b[~zero_rows]
                        )
                else:
                    res = mod2Solve_m4ri_py(A[~zero_rows, :], b[~zero_rows])
            dists.append(1 - int(res))
        elif regression_mode == "approx":
            res = sp_ridge_regression_mod2(
                one_ridx_A_local,
                one_cidx_A_local,
                nrow_A,
                ncol_A_local,
                b,
                ridge_coef_a,
                ridge_coef_b,
            )
            if res is None:
                dists.append(np.nan)
            else:
                A, s = res
                pred = np.round(A.dot(s)) % 2
                if np.sum(np.logical_or(pred, b)) == 0:
                    dist = 0
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
