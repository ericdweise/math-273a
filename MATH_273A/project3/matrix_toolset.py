import numpy as np


def is_sym_pos_def(A):
    """Test if a matrix A is symmetric positive definite"""
    for i in range(A.shape[0]):
        if not np.linalg.det(A[:i+1,:i+1])>0:
            return False

        for j in range(A.shape[1]):
            if i != j:
                if A[i,j] != A[j,i]:
                    return False

    return True


def iter_matrix(A):
    for row in range(A.shape[0]):
        for col in range(A.shape[1]):
            yield row, col
