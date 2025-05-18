import numpy as np
from scipy.sparse import csr_matrix, diags
from sksparse.cholmod import cholesky


def sp_ridge_regression_mod2(one_ridx_A, one_cidx_A, nrow_A, ncol_A, b, ridge_coef):
    # figure out number of ones in each row
    _, n_ones_per_row = np.unique(one_ridx_A, return_counts=True)
    # the new columns are used to perform 1 + 1 - 2 = 0
    data = np.concatenate(
        [np.ones(len(one_ridx_A)), np.ones(nrow_A) * (n_ones_per_row // 2 + 1) * -2]
    )
    A = csr_matrix(
        (
            data,
            (
                np.concatenate([one_ridx_A, np.arange(nrow_A)]),
                np.concatenate([one_cidx_A, np.arange(nrow_A) + ncol_A]),
            ),
        ),
        shape=(nrow_A, ncol_A + nrow_A),
    )
    factor = cholesky(
        A.transpose().dot(A)
        + diags(np.repeat(ridge_coef, ncol_A + nrow_A)) * ridge_coef
    )
    return (
        A,
        factor(
            A.transpose().dot(b)
            + np.concatenate(
                [
                    np.repeat(0.5, ncol_A),
                    np.repeat(0, nrow_A),
                ]
            )
            * ridge_coef
        ),
    )
