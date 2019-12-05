from math import pi
from math import cos
from math import sin
from math import sqrt

from matrix_toolset import create_mesh

# CONSTANTS
PARTICLE_H = pi/250


def particle(j):
    theta = j*PARTICLE_H
    x = pi/30 + 0.7*cos(theta)
    y = sqrt(2)/5 + 0.4*sin(theta)
    q = x + 1
    return x, y, q


def direct_sum_contribution(x,y,q):
    return 


if __name__ == '__main__':

    # define mesh
    xgrid, ygrid = create_mesh(gridstep=0.02)

    # define atom locations
    # for each atom add contribution
    # iterate through grid points in solution mesh
    # compute contribution from each point in direct contribution region
