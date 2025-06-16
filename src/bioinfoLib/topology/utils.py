import warnings

import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.spatial.distance import hamming
from shapely import LineString, frechet_distance
from sksparse.cholmod import cholesky


def sp_ridge_regression_mod2(
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
    # weight matrix
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        W = diags(
            b * (len(b) / np.sum(b) - len(b) / (len(b) - np.sum(b)))
            + len(b) / (len(b) - np.sum(b))
        )
    # W = diags(np.repeat(1, nrow_A))
    # as suggested in sksparse's docment, make A csr so that its transpose is csc
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
    B = A.transpose().dot(W).dot(A) + diags(
        np.concatenate(
            [
                np.repeat(ridge_coef_a, ncol_A),
                np.repeat(ridge_coef_b, nrow_A_new),
            ]
        )
    )
    factor = cholesky(B)
    best_diff = np.Inf
    best_s = None
    for i in np.linspace(0, 1, n_search):
        for j in np.linspace(0, 1, n_search):
            s = factor(
                A.transpose().dot(W.dot(b))
                + np.concatenate(
                    [
                        np.repeat(ridge_coef_a * i, ncol_A),
                        np.repeat(ridge_coef_b * j, nrow_A_new),
                    ]
                )
            )
            s[s >= 0.5] = 1
            s[s < 0.5] = 0
            pred = A[:, :ncol_A].dot(s[:ncol_A]) % 2
            diff = hamming(pred, b)
            if diff < best_diff:
                best_diff = diff
                best_s = s
    return (A[:, :ncol_A], best_s[:ncol_A])


def edge_idx_encode(i, j):
    if i > j:
        i_tmp = i
        i = j
        j = i_tmp
    return ((i + j) * (i + j + 1)) // 2 + j


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
