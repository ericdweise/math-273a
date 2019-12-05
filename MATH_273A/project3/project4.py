import numpy as np

from math import pi
from math import cos
from math import sin
from math import sqrt
from math import log

from matrix_toolset import create_mesh
from matrix_toolset import iter_matrix
from plot_tools import plot_heatmap
from plot_tools import plot_surface

# CONSTANTS
PARTICLE_H = pi/250
GRIDSTEP = 0.02


def particle(j):
    theta = j*PARTICLE_H
    x = pi/30 + 0.7*cos(theta)
    y = sqrt(2)/5 + 0.4*sin(theta)
    q = x + 1
    return x, y, q


def direct_sum_contribution(x,y,xq,yq,q):
    return -1*q*log(sqrt((x-xq)**2 + (y-yq)**2))


if __name__ == '__main__':

    # define mesh
    xgrid, ygrid = create_mesh(gridstep=GRIDSTEP)
    output = np.zeros(xgrid.shape)

    # define atom locations
    atoms = []
    for j in range(1,501):
        atoms.append(particle(j))

    # iterate through grid points in solution mesh
    for i,j in iter_matrix(xgrid):
        x = xgrid[i,j]
        y = ygrid[i,j]

        for xq,yq,q in atoms:
            output[i,j] += direct_sum_contribution(x,y,xq,yq,q)

    plot_heatmap(xgrid, ygrid, output, 'results/project4/heatmap.png')
    plot_surface(xgrid, ygrid, output, 'results/project4/surface.png')
