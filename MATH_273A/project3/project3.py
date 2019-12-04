import numpy as np
from math import cos
from math import exp
from math import pi
from math import sin
from math import sqrt
from matrix_toolset import iter_matrix
from cholesky_factorization import incomplete_cholesky_factorization
from congugate_gradient import congugate_gradient


def create_mesh():
    xmin = -1
    ymin = -1
    xmax = 1
    ymax = 1
    gridstep = 0.01
    xlist = np.arange(xmin, xmax, gridstep)
    ylist = np.arange(ymin, ymax, gridstep)
    xgrid, ygrid = np.meshgrid(xlist, ylist)
    return xgrid, ygrid


def define_regions(xcent, a, ycent, b):
    """Create a function that defines the interior of an ellipse"""

    def func(x,y):
        """Given a point returns what region the point is in:
        OMEGA: The point is in the ellipse interior
        dOMEGA: The point is on the edge of the ellipse
        D: The point is outside of the ellipse but not on dD
        dD: The point is on the boundary of the computational domain
        """
        if (x==-1) or (x==1) or (y==-1) or (y==-1):
            return 'dD'

        v = ((x-xcent)/a)**2 + ((y-ycent)/b)**2
        if v < 1:
            return 'OMEGA'
        if v == 1:
            return 'dOMEGA'

        return 'D'

    return(func)


def func_f(x,y):
    """Apply function f to a point"""
    return(4*exp(y)*sin(pi*x))


def func_g(x,y):
    """Apply function g to a point"""
    return(exp(x)*cos(2*pi*y))


def build_soe_matrix(find_region, xgrid, ygrid):
    """Will create a matrix that represents the system of equations"""

    ashape = xgrid.shape[0]*xgrid.shape[1]
    A = np.zeros((ashape, ashape))

    for row, col in iter_matrix(xgrid):
        region = find_region(row, col)


def add_ghost_fluid_point():
    pass


def compute_residual(A, b, x):
    return np.linalg.norm(b-np.matmul(A,x), 2)/np.linalg.norm(b, 2)


if __name__ == '__main__':

    # Cutoff point
    residual = 10**(-7)

    # DEFINE ELLIPSE
    XCENTER = 0.5
    YCENTER = 0.5
    A = 0.1
    B = 0.25
    find_region = define_regions(XCENTER, A, YCENTER, B)

    # CREATE MATRIX FOR SYSTEM OF EQUATIONS
    A = build_soe_matrix(find_region, func_f, func_g)

    # Find Incomplete Cholesky Factorization
    R = incomplete_cholesky_factorization(A)
    R_inverse = np.linalg.inv(R)
    R_transpose_inverse = np.linalg.inv(np.matrix.transpose(R))

    # Precondition Matrix:
    A_hat = np.matmul(np.matmul(R_transpose_inverse, A), R_inverse)
    b_hat = np.matmul(R_transpose_inverse, b)
    x_hat = congugate_gradient()

    x = np.matmul(R_inverse, x_hat)
