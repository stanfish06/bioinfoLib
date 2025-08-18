import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial.distance import hamming

from bioinfoLib.topology.utils import sp_ridge_regression_mod2


class Test_ridge_regression_mod2:
    @staticmethod
    def create_random_matrix(nrow, ncol, density=0.1, seed=1):
        np.random.seed(seed)
        n_ones = int(nrow * ncol * density)
        rows = np.random.choice(nrow, n_ones)
        cols = np.random.choice(ncol, n_ones)
        unique_pairs = np.unique(np.column_stack([rows, cols]), axis=0)
        rows, cols = unique_pairs[:, 0], unique_pairs[:, 1]
        return (rows, cols, nrow, ncol)

    @staticmethod
    def create_problems(
        one_ridx_A,
        one_cidx_A,
        nrows_A,
        ncols_A,
        solution_density=0.1,
        seed=1,
        n_problems=1,
    ):
        np.random.seed(seed)
        A = csr_matrix(
            (np.ones(len(one_ridx_A)), (one_ridx_A, one_cidx_A)),
            shape=(nrows_A, ncols_A),
        )
        n_ones_solution = int(nrows_A * solution_density)
        solutions = []
        targets = []
        for _ in range(n_problems):
            one_idx_solution = np.random.choice(ncols_A, n_ones_solution)
            s = np.zeros(ncols_A)
            s[one_idx_solution] = 1
            b = A.dot(s) % 2
            solutions.append(s)
            targets.append(b)
        return (solutions, targets)

    def test_small_cases(self):
        nrows_A = 6
        ncols_A = 6
        (one_ridx_A, one_cidx_A, _, _) = (
            Test_ridge_regression_mod2.create_random_matrix(
                nrows_A, ncols_A, density=0.25, seed=1
            )
        )
        assert len(one_ridx_A) == len(one_cidx_A)
        assert np.max(one_ridx_A) < nrows_A and np.min(one_ridx_A) >= 0
        assert np.max(one_cidx_A) < ncols_A and np.min(one_cidx_A) >= 0

        (solutions, targets) = Test_ridge_regression_mod2.create_problems(
            one_ridx_A,
            one_cidx_A,
            nrows_A,
            ncols_A,
            solution_density=0.5,
            seed=1,
            n_problems=10,
        )

        A_dense = np.zeros((nrows_A, ncols_A), dtype=int)
        A_dense[one_ridx_A, one_cidx_A] = 1
        print("Input matrix A:")
        print(A_dense)

        ridge_coef_a = 0.1
        ridge_coef_b = 0.1
        for s, b in zip(solutions, targets):
            print("Try to solve:")
            print(b)
            print("One solution is:")
            print(s)
            A_sparse, s_pred = sp_ridge_regression_mod2(
                one_ridx_A,
                one_cidx_A,
                nrows_A,
                ncols_A,
                b,
                ridge_coef_a,
                ridge_coef_b,
            )
            print(f"Solution vector: {s_pred.astype(int)}")
            pred = A_sparse.dot(s_pred) % 2
            dist = hamming(pred, b)
            print(f"Prediction A*x mod 2: {pred.astype(int)}")
            print(f"Match: {['✓' if p == t else '✗' for p, t in zip(pred, b)]}")
            print(f"Hamming distance: {dist:.3f}")
            assert dist <= 0.2, f"Accuracy is: {1 - dist}"

    def test_medium_cases(self):
        nrows_A = 40
        ncols_A = 40
        (one_ridx_A, one_cidx_A, _, _) = (
            Test_ridge_regression_mod2.create_random_matrix(
                nrows_A, ncols_A, density=0.25, seed=1
            )
        )
        assert len(one_ridx_A) == len(one_cidx_A)
        assert np.max(one_ridx_A) < nrows_A and np.min(one_ridx_A) >= 0
        assert np.max(one_cidx_A) < ncols_A and np.min(one_cidx_A) >= 0

        (solutions, targets) = Test_ridge_regression_mod2.create_problems(
            one_ridx_A,
            one_cidx_A,
            nrows_A,
            ncols_A,
            solution_density=0.1,
            seed=1,
            n_problems=10,
        )

        A_dense = np.zeros((nrows_A, ncols_A), dtype=int)
        A_dense[one_ridx_A, one_cidx_A] = 1
        print("Input matrix A:")
        print(A_dense)

        ridge_coef_a = 0.1
        ridge_coef_b = 10
        for s, b in zip(solutions, targets):
            print("Try to solve:")
            print(b)
            print("One solution is:")
            print(s)
            A_sparse, s_pred = sp_ridge_regression_mod2(
                one_ridx_A,
                one_cidx_A,
                nrows_A,
                ncols_A,
                b,
                ridge_coef_a,
                ridge_coef_b,
            )
            print(f"Solution vector: {s_pred.astype(int)}")
            pred = A_sparse.dot(s_pred) % 2
            dist = hamming(pred, b)
            print(f"Prediction A*x mod 2: {pred.astype(int)}")
            print(f"Match: {['✓' if p == t else '✗' for p, t in zip(pred, b)]}")
            print(f"Hamming distance: {dist:.3f}")
            assert dist <= 0.25, f"Accuracy is: {1 - dist}"
