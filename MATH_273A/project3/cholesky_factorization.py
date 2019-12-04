import numpy as np
from math import sqrt


def incomplete_cholesky_factorization(A):
    """Compute the sparse Upper Right Triangular matrix R in R^T*R=A"""
    R = np.zeros(A.shape)
    for i in range(A.shape[0]):
        r = A[i,i]
        for k in range(i-1):
            r -= R[k,i]**2
        R[i,i] = sqrt(r)

        for j in range(i, A.shape[1]):
            # "Incomplete" factorization
            if A[i,j] == 0:
                continue

            r = A[i,j]
            for k in range(i-1):
                r -= R[k,i]*R[k,j]
            R[i,j] = r/R[i,i]

    return R
