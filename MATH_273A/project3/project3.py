import numpy as np
from math import cos
from math import exp
from math import pi
from math import sin
from math import sqrt


def create_grids():
    xmin = -1
    ymin = -1
    xmax = 1
    ymax = 1
    gridstep = 0.01
    xlist = np.arange(xmin, xmax, gridstep)
    ylist = np.arange(ymin, ymax, gridstep)
    xgrid, ygrid = np.meshgrid(xlist, ylist)
    return xgrid, ygrid


def define_ellipse_func(xcent, a, ycent, b):
    """Create a function that defines the interior of an ellipse"""
    def func(x,y):
        """returns True if point (x,y) is inside of the ellipse. False otherwise."""
        if ((x-xcent)/a)**2 + ((y-ycent)/b)**2 < 1:
            return True
        else:
            return False
    return(func)


def func_f(x,y):
    """Apply function f to a point"""
    return(4*exp(y)*sin(pi*x))


def func_g(x,y):
    """Apply function g to a point"""
    return(exp(x)*cos(2*pi*y))


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


def build_soe_matrix(ellipse_func, xgrid, ygrid):
    """Will create a matrix that represents the system of equations"""
    pass


def add_ghost_fluid_point():
    pass


def compute_residual():
    pass


if __name__ == '__main__':

    # DEFINE ELLIPSE
    XCENTER = 0.5
    YCENTER = 0.5
    A = 0.1
    B = 0.25
    is_in_ellipse = define_ellipse_func(XCENTER, A, YCENTER, B)

    # CREATE MATRIX FOR SYSTEM OF EQUATIONS
    soe_matrix = build_soe_matrix(is_in_ellipse, func_f, func_g)
