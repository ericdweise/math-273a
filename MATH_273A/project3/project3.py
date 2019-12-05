import numpy as np
import sys

from math import acos
from math import cos
from math import exp
from math import pi
from math import sin
from math import sqrt

from matrix_toolset import create_mesh
from matrix_toolset import iter_matrix
from matrix_toolset import is_sym_pos_def
from matrix_toolset import is_symmetric
from plot_tools import plot_heatmap
from plot_tools import plot_surface

from cholesky_factorization import incomplete_cholesky_factorization

from conjugate_gradient import conjugate_gradient


np.set_printoptions(threshold=sys.maxsize)


class Ellipse(object):

    def __init__(self, xcent, a, ycent, b):
        self.xcent = xcent
        self.a = a
        self.ycent = ycent
        self.b = b

    def inside(self, x, y):
        """Will return True if (x,y) is inside the ellipse. False otherwise"""
        if ((x-self.xcent)/self.a)**2 + ((y-self.ycent)/self.b)**2 < 1:
            return True
        return False

    def find_alpha(self, x1, y1, x2, y2):
        """Find distance between (x1, y1) and boundary of elipse in direction of (x2,y2)"""
        n_steps = 4
        inside1 = self.inside(x1, y1)
        inside2 = self.inside(x2, y2)
        assert(inside1 != inside2)

        x_step = x2-x1
        y_step = y2-y1
        DIST = sqrt((x1-x2)**2 + (y1-y2)**2)

        for _ in range(n_steps):

            x_step = x_step/2.0
            y_step = y_step/2.0
            inside2 = self.inside(x2, y2)

            if inside2==inside1:
                # move away from x1,y1
                x2 = x2 + x_step
                y2 = y2 + y_step
            else:
                # move towards x1,y1
                x2 = x2 - x_step
                y2 = y2 - y_step

        alpha = sqrt((x2-x1)**2 + (y2-y1)**2)/DIST
        assert(alpha<1), f'alpha is greater than 1: ({x1},{y1}), ({x2},{y2})'

        return alpha


def func_f(x,y):
    """Apply function f to a point"""
    return(4*exp(y)*sin(pi*x))


def func_g(x,y):
    """Apply function g to a point"""
    return(exp(x)*cos(2*pi*y))


def build_soe_matrix(omega, xgrid, ygrid, gridstep):
    """Creates the matrix representation of the system of equations"""
    ashape = xgrid.shape[0]*xgrid.shape[1]
    A = np.zeros((ashape, ashape))
    b = np.zeros(ashape)

    # Gridpoints NOT on Boundary
    for i,j in iter_matrix(xgrid):

        # skip boundary points
        if i in (0,xgrid.shape[0]-1):
            continue
        if j in (0, xgrid.shape[1]-1):
            continue

        x_ij = xgrid[i,j]
        y_ij = ygrid[i,j]
        row = (i*xgrid.shape[0])+j # Row index for A, b

        # Points not on boundary of computational domain
        inside_ij = omega.inside(x_ij,y_ij)

        if inside_ij:
            diff_func = func_f
        else:
            diff_func = func_g

        # Set diagonal value
        A[row,row] = -4*gridstep**(-2)

        # Set value in b vector
        b[row] = diff_func(x_ij,y_ij)

        # Add neighbor dependencies
        for (k,l) in [(i+1,j),(i,j+1),(i-1,j),(i,j-1)]:
            col = (k*xgrid.shape[0])+l # Column index for A
            x_kl = xgrid[k,l]
            y_kl = ygrid[k,l]
            inside_kl = omega.inside(x_kl,y_kl)

            # Neighbor in same region as center point
            if inside_ij==inside_kl:
                A[row,col] += gridstep**(-2)

            # Neighbor in different region from center point
            else:
                if inside_ij:
                    boundary_value = 2
                else:
                    boundary_value = 1

                alpha = omega.find_alpha(x_ij, y_ij, x_kl, y_kl)

                A[row,row] -= gridstep**(-2)/alpha
                b[row] -= boundary_value*gridstep**(-2)/alpha

    # Grid Points on BOUNDARY
    bdy = []
    bdy.extend([(0, j) for j in range(xgrid.shape[1])])
    bdy.extend([(i, 0) for i in range(xgrid.shape[0])])
    bdy.extend([(xgrid.shape[0]-1, j) for j in range(xgrid.shape[1])])
    bdy.extend([(i, xgrid.shape[1]-1) for i in range(xgrid.shape[0])])

    for i,j in bdy:
        row = (i*xgrid.shape[0])+j # Row index for A, b
        A[row,:] = 0 # Remove dependence of other points on this point
        A[:,row] = 0
        A[row,row] = -1 # Direct dependence on 
        b[row] = 0 # Negative u value at boundary

    return A, b


if __name__ == '__main__':

    # distance between grid points
    GRIDSTEP = 0.05

    # Residual cutoff point
    cutoff = 10**(-7)

    # DEFINE ELLIPSE
    XCENTER = 0.
    YCENTER = 0.1
    ELLIPSE_A = 0.5
    ELLIPSE_B = 0.25
    ellipse = Ellipse(XCENTER, ELLIPSE_A, YCENTER, ELLIPSE_B)

    # Define grids
    xgrid, ygrid = create_mesh(gridstep=GRIDSTEP )

    # CREATE MATRIX FOR SYSTEM OF EQUATIONS
    A, b = build_soe_matrix(ellipse, xgrid, ygrid, GRIDSTEP)
    with open('results/project3/A.np', 'w') as f:
        A.tofile(f)
    with open('results/project3/b.np', 'w') as f:
        b.tofile(f)

    assert(is_symmetric(A))

    # Find Incomplete Cholesky Factorization
    A_neg = -1*A
    b_neg = -1*b
    # assert(is_sym_pos_def(A))
    print('Find Incomplete Cholesky Factor, R')
    R = incomplete_cholesky_factorization(A_neg)
    print('Finding R inverse')
    R_inverse = np.linalg.inv(R)
    print('Finding Inverse of Transpose of R')
    R_transpose_inverse = np.linalg.inv(np.matrix.transpose(R))

    # Precondition Matrix:
    print('Preconditioning Matrix')
    A_hat = np.matmul(np.matmul(R_transpose_inverse, A_neg), R_inverse)
    b_hat = np.matmul(R_transpose_inverse, b_neg)

    print('Performing conjugate gradient')
    x_hat = conjugate_gradient(A_hat, b_hat, cutoff)

    x = -1*np.matmul(R_inverse, x_hat)
    with open('results/project3/x.np', 'w') as f:
        x.tofile(f)

    u = np.zeros(xgrid.shape)
    for k in range(x.shape[0]):
        i = k % xgrid.shape[0]
        j = k // xgrid.shape[0]
        u[i,j] = x[k]
    
    plot_heatmap(xgrid, ygrid, u, 'results/project3/heatmap.png')
    plot_surface(xgrid, ygrid, u, 'results/project3/surface.png')
