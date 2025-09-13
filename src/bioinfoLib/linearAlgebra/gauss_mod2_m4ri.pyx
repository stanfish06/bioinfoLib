cdef extern from "m4ri/m4ri.h":
    ctypedef int rci_t
    ctypedef int wi_t
    ctypedef unsigned long long m4ri_word "word"
    ctypedef int BIT
    ctypedef struct mzd_t:
        rci_t nrows
        rci_t ncols
        wi_t width

    cdef mzd_t *mzd_init(rci_t, rci_t)
    cdef void mzd_free(mzd_t *)
    cdef void mzd_write_bit(mzd_t *m, rci_t row, rci_t col, BIT value)

    cdef void mzd_print(mzd_t *)
    cdef int mzd_solve_left(mzd_t *A, mzd_t *B, int cutoff, int inconsistency_check)


def mod2Solve_m4ri_py(A, b):
    nr = len(A)
    nc = len(A[0])

    assert nr >= nc, "number of rows must be greater than or equal to the number of columns" 

    cdef mzd_t *A2 = mzd_init(nr, nc);
    cdef mzd_t *b2 = mzd_init(nr, 1);
    
    try:
        for i in range(nr):
            for j in range(nc):
                mzd_write_bit(A2, i, j, A[i][j])
        for i in range(nr):
            mzd_write_bit(b2, i, 0, b[i])
        result = mzd_solve_left(A2, b2, 0, 1)
        return result == 0
    finally:
        mzd_free(A2)
        mzd_free(b2)

