# distutils: language = c++
from libcpp.vector cimport vector
from libcpp cimport bool

cdef bool mod2solve(vector[vector[bool]] Ab):
    nr = Ab.size()
    nc = Ab[0].size()

    # If there is a row of 0s and the rhs is 1
    # we know it is not solvable
    for i in range(0, nr):
        if Ab[i][nc - 1]:
            allzero = True;
            for j in range(0, nc - 1):
                if Ab[i][j]:
                    allzero = False;
                    break;
            if allzero:
                return False;

    # this stores the next pivot
    row = 1
    for i in range(0,nc-1):
        pivot = -1
        for j in range(row, nr):
            if Ab[j][i]:
                pivot = j
                break
        if pivot != -1:
            if pivot != row:
                Ab[row].swap(Ab[pivot])
            for j in range(0, nr):
                if j != row and Ab[j][i]:
                    for k in range(0, nc):
                        Ab[j][k] = Ab[j][k] ^ Ab[row][k];
            if row == nr:
                return True
            row = row + 1;

    for i in range(row - 1, nr):
        if Ab[i][nc - 1]:
            allzero = True;
            for j in range(0, nc - 1):
                if Ab[i][j]:
                    allzero = False;
                    break;
            if allzero:
                return False;
    return True;

# A should be 2d list of bool values
def solve_mod2(A, b):
    cdef vector[vector[bool]] Ab
    cdef vector[bool] row

    for ai, bi in zip(A, b):
        row = vector[bool]()
        for aij in ai:
            row.push_back(aij)
        row.push_back(bi)
        Ab.push_back(row)

    return mod2Solve(Ab)

