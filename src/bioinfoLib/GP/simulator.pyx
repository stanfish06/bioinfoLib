from libc.math cimport exp, pow, sqrt, log, sin, cos, M_PI
from libc.stdlib cimport rand, RAND_MAX
import numpy as np
from numpy.linalg import cholesky

# for fun, this is not the most efficient to sample normal random variable
cdef double box_muller_rnorm(): 
    cdef double u1 = rand() / RAND_MAX
    cdef double u2 = rand() / RAND_MAX
    return sqrt(-2*log(u1))*cos(2*M_PI*u2)

def gp_1d(mean_func, kernel_func, xs, l, sig):
    # this probably won't speed things up a lot
    # consider rewirte this to c code
    n = len(xs)
    cov_mat = np.zeros([n, n])
    for i in range(n):
        for j in range(i,n):
            val = kernel_func(xs[i], xs[j], l, sig)
            cov_mat[i,j] = val 
            cov_mat[j,i] = val
    A = cholesky(cov_mat)
    mu = np.array([mean_func(x) for x in xs])
    zs = np.zeros(n)
    for i in range(n):
        zs[i] = box_muller_rnorm()
    ys = mu + np.dot(A,zs)
    return ys


# 1d squared-exponential kernel
def sq_exp_kernel_1d(a, b, l, sig):
    return sig*np.exp(-np.power(a-b,2)/(2*np.power(l,2)))
