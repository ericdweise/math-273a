import numpy as np


def create_mesh(gridstep=0.01, xmin=-1, xmax=1, ymin=-1, ymax=1):
    xlist = np.arange(xmin, xmax + 0.01*gridstep, gridstep)
    ylist = np.arange(ymin, ymax + 0.01*gridstep, gridstep)
    xgrid, ygrid = np.meshgrid(xlist, ylist)
    return xgrid, ygrid


def is_symmetric(A):
    for i in range(1, A.shape[0]):
        for j in range(i, A.shape[1]):
            if i != j:
                if A[i,j] != A[j,i]:
                    return False
    return True


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


def index_vec_to_mat(idx, mat_width):
    return idx // mat_width, idx % mat_width


def index_mat_to_vec(row, col, mat_width):
    return row*mat_width + col
