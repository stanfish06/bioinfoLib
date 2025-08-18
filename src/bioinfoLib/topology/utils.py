import itertools

import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.spatial.distance import hamming
from shapely import LineString, frechet_distance
from sksparse.cholmod import cholesky


# TODO: grid search, adjust radius based on gaussian prior
# TODO: generate a set of representative loops through custom distance function
# TODO: pair up loops from two sets though frechet distance and do regression based matching
def sp_ridge_regression_mod2(
    one_ridx_A,
    one_cidx_A,
    nrow_A,
    ncol_A,
    b,
    ridge_coef_a,
    ridge_coef_b,
    n_neighbors=10,
    do_approximation=False,
    n_search_cutoff=10,
    max_bits_flip=10,
    max_bits_comb=10,
    n_post_process=10,
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
):
    if not do_regression:
        dist = min(
            np.min(
                [
                    frechet_distance(
                        LineString(
                            np.concatenate(
                                [source_loop_coords[i:, :], source_loop_coords[0:i, :]]
                            )
                        ),
                        LineString(target_loop_coords),
                        densify=0.5,
                    )
                    for i in range(source_loop_coords.shape[0])
                ]
            ),
            np.min(
                [
                    frechet_distance(
                        LineString(
                            np.concatenate(
                                [source_loop_coords[i:, :], source_loop_coords[0:i, :]]
                            )[::-1, :]
                        ),
                        LineString(target_loop_coords),
                        densify=0.5,
                    )
                    for i in range(source_loop_coords.shape[0])
                ]
            ),
        )
        if dist > max_frechet_dist:
            return None
        else:
            return dist
    else:
        # start_time = time.time()
        b1 = np.zeros(nrow_A)
        b1[source_loop_edges] = 1
        b2 = np.zeros(nrow_A)
        b2[target_loop_edges] = 1
        b = np.logical_xor(b1, b2).astype(int)

        # remove columns with birth t larger than the birth t of the loop
        columns_kept = np.where(bd_column_birth_t <= source_loop_birth_t)[0]
        if columns_kept.size == 0:
            return np.nan
        mask = ~np.isin(one_cidx_A, columns_kept)
        one_ridx_A = one_ridx_A[~mask]
        one_cidx_A = one_cidx_A[~mask]
        ncol_A = np.max(one_cidx_A) + 1
        A, s = sp_ridge_regression_mod2(
            one_ridx_A,
            one_cidx_A,
            nrow_A,
            ncol_A,
            b,
            ridge_coef_a,
            ridge_coef_b,
            n_search,
            n_neighbors,
            do_approximation,
        )
        pred = np.round(A.dot(s)) % 2
        # ham_dist = hamming(pred, b) * nrow_A / np.sum(b)
        if np.sum(np.logical_or(pred, b)) == 0:
            return np.nan
        ham_dist = 1 - np.sum(np.logical_and(pred, b)) / np.sum(np.logical_or(pred, b))

        # end_time = time.time()
        # print(f"\ndone in {end_time - start_time}s")
        return ham_dist


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
