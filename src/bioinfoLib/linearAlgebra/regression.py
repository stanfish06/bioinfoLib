import numpy as np
from scipy.sparse import csr_matrix, diags
from sksparse.cholmod import cholesky_AAt


def sp_ridge_regression_mod2(one_ridx_A, one_cidx_A, nrow_A, ncol_A, b, ridge_coef):
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
    W = diags(
        b * (len(b) / np.sum(b) - len(b) / (len(b) - np.sum(b)))
        + len(b) / (len(b) - np.sum(b))
    )
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
    factor = cholesky_AAt(A.transpose().dot(np.sqrt(W)), beta=ridge_coef)
    return (
        A,
        factor(A.transpose().dot(W.dot(b))),
    )
