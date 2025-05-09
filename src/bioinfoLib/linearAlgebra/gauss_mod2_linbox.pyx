ctypedef int size_t
        
cdef extern from "givaro/modular.h" namespace "Givaro":
    cdef cppclass Element_uint64 "Givaro::Modular<uint64_t>::Element":
        pass
    cdef cppclass Modular_uint64 "Givaro::Modular<uint64_t>":
        Modular_uint64(int modulus) except +
        const Element_uint64 one
        Element_uint64 assign(Element_uint64 &x, const Element_uint64 &y)

    cdef cppclass Element_double "Givaro::Modular<double>::Element":
        pass
    cdef cppclass Modular_double "Givaro::Modular<double>":
        Modular_double(int modulus) except +
        const Element_double one
        Element_double assign(Element_double &x, const Element_double &y)

cdef extern from "linbox/matrix/sparse-matrix.h" namespace "LinBox":
    cdef cppclass SparseMatrix_Modular_uint64 "LinBox::SparseMatrix<Givaro::Modular<uint64_t>, LinBox::SparseMatrixFormat::SparseSeq>":
        ctypedef Modular_uint64 Field
        ctypedef Element_uint64 Element
        SparseMatrix_Modular_uint64(Field &F, size_t m, size_t n) except +
        const Element& setEntry(size_t i, size_t j, const Element &value)
        
    cdef cppclass SparseMatrix_Modular_double "LinBox::SparseMatrix<Givaro::Modular<double>, LinBox::SparseMatrixFormat::SparseSeq>":
        ctypedef Modular_double Field
        ctypedef Element_double Element
        SparseMatrix_Modular_double(Field &F, size_t m, size_t n) except +
        const Element& setEntry(size_t i, size_t j, const Element &value)
        
cdef extern from "linbox/vector/vector.h" namespace "LinBox":
    cdef cppclass DenseVector_Modular_uint64 "LinBox::DenseVector<Givaro::Modular<uint64_t>>":
        ctypedef Modular_uint64 Field
        ctypedef Element_uint64 Element
        DenseVector_Modular_uint64(Field &F, size_t n) except +       
        Element& operator[](size_t i)

    cdef cppclass DenseVector_Modular_double "LinBox::DenseVector<Givaro::Modular<double>>":
        ctypedef Modular_double Field
        ctypedef Element_double Element
        DenseVector_Modular_double(Field &F, size_t n) except +       
        Element& operator[](size_t i)

cdef extern from "linbox/solutions/methods.h" namespace "LinBox::Method":
    cdef cppclass SparseElimination "LinBox::Method::SparseElimination":
        SparseElimination() except +
    cdef cppclass Wiedemann "LinBox::Method::Wiedemann":
        Wiedemann() except +

cdef extern from "linbox/solutions/solve.h" namespace "LinBox":
    # Declare the specific solve function template instantiation
    cdef cppclass ModularTag "LinBox::RingCategories::ModularTag":
        ModularTag() except +

    DenseVector_Modular_double& solve(
        DenseVector_Modular_uint64& x,
        const SparseMatrix_Modular_uint64& A,
        const DenseVector_Modular_uint64& b,
        const ModularTag& tag,
        const SparseElimination& m
    ) except +

    DenseVector_Modular_double& solve(
        DenseVector_Modular_double& x,
        const SparseMatrix_Modular_double& A,
        const DenseVector_Modular_double& b,
        const ModularTag& tag,
        const Wiedemann& m
    ) except +

def mod2Solve_linbox_py(one_ridx_A, one_cidx_A, nrow_A, ncol_A, one_i_b):
    cdef Modular_uint64 * F2 = new Modular_uint64(2)
    cdef SparseMatrix_Modular_uint64* A = new SparseMatrix_Modular_uint64(F2[0], nrow_A, ncol_A)

    for r, c in zip(one_ridx_A, one_cidx_A):
        A.setEntry(r, c, F2[0].one)

    cdef DenseVector_Modular_uint64 * b = new DenseVector_Modular_uint64(F2[0], nrow_A)

    cdef DenseVector_Modular_uint64 * x = new DenseVector_Modular_uint64(F2[0], ncol_A)

    for i in one_i_b:
        F2.assign(b[0][i], F2[0].one);

    cdef ModularTag tag;
    cdef SparseElimination method;
    # cdef Wiedemann method;
    solved = True
    try:
        solve(x[0], A[0], b[0], tag, method)
    except:
        solved = False
    # cython will translate del to delete
    del A, x, b, F2
    return solved
