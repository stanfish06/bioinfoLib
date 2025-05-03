ctypedef int size_t
        
cdef extern from "givaro/modular.h" namespace "Givaro":
    cdef cppclass Modular_uint64 "Givaro::Modular<uint64_t>":
        Modular_uint64(int modulus)
    cdef cppclass Modular_double "Givaro::Modular<double>":
        Modular_double(int modulus) except +

cdef extern from "linbox/matrix/sparse-matrix.h" namespace "LinBox":
    cdef cppclass SparseMatrix_F2 "LinBox::SparseMatrix<Givaro::Modular<double>>":
        ctypedef Modular_double Field
        SparseMatrix_F2(const Field& F, size_t m, size_t n) except +
    cdef cppclass SparseMatrix_Modular_uint64 "LinBox::SparseMatrix<Givaro::Modular<uint64_t>, LinBox::SparseMatrixFormat::SparseSeq>":
        ctypedef Modular_uint64 Field
        SparseMatrix_Modular_uint64(Field &F, size_t m, size_t n)

def mod2Solve_linbox():
    # pass
    # cdef Modular_double* F2 = new Modular_double(2)
    # cdef SparseMatrix_F2* A = new SparseMatrix_F2(F2[0], 4, 4)
    cdef Modular_uint64 * F2 = new Modular_uint64(2)
    cdef SparseMatrix_Modular_uint64* A = new SparseMatrix_Modular_uint64(F2[0], 4, 4)