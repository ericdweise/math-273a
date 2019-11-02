"""
Project 2 for Math 273A: Level set functions.

To run this code use the following command:
    python3 project2.py

Package Requirements:
    numpy: pip install numpy
    matplotlib: pip install matplotlib

Instructor: Li-Tien Cheng
Academic Quarter: Fall 2019
Author/student: Eric Weise
"""

import numpy as np
from math import floor, sqrt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

TEST = True


class Phi(object):

    def __init__(self, xmin=-1, xmax=1, ymin=-1, ymax=1, grid_step=0.02, eye_x=0, eye_y=0):
        """Initialize grid with values of float('inf')"""
        self.eye_x = eye_x
        self.eye_y = eye_y

        self.xmin = xmin
        self.xmax = xmax

        self.ymin = ymin
        self.ymax = ymax

        self.grid_step = grid_step

        # 1D coordinate axes
        xlist = np.arange(xmin, xmax + grid_step/2, grid_step)
        ylist = np.arange(ymin, ymax + grid_step/2, grid_step)

        # 2D grid to hold x and y values
        self.xgrid, self.ygrid = np.meshgrid(xlist, ylist)

        # 2D grid to hold Phi values
        # Initialize with Inf to accommodate min() for union
        self.phigrid = np.full(self.xgrid.shape, np.inf)

        # Set flow field
        denom = np.sqrt((self.xgrid - eye_x)**2 + (self.ygrid - eye_y)**2)
        with np.nditer(denom, op_flags=['readwrite']) as it:
            for x in it:
                if x == 0:
                    x[...] = np.inf

        self.flowx = np.divide(self.xgrid, denom)
        self.flowy = np.divide(self.ygrid, denom)


    def plot_phi(self, path):
        """Will plot the values stored in the grid.
        Args:
            path: File where the plot will be saved.
        """
        # https://matplotlib.org/mpl_toolkits/mplot3d/tutorial.html#surface-plots
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(self.xgrid, self.ygrid, self.phigrid)

        # Set axis limits
        ax.set_xlim(self.xmin, self.xmax)
        ax.set_ylim(self.ymin, self.ymax)
        ax.set_zlim(-1, 0)

        # Set number of ticks along axes
        plt.locator_params(axis='x', nbins=3)
        plt.locator_params(axis='y', nbins=3)
        plt.locator_params(axis='z', nbins=3)

        # set axis labels
        plt.xlabel('X')
        plt.ylabel('Y')

        plt.savefig(path)

        plt.cla()
        plt.clf()
        plt.close()


    def plot_level_set(self, path, value=0):
        """Create a 2D plot of Phi(x,y)=value.
        Args:
            value: The value at which the level set is comuputed.
        """
        fig = plt.figure()
        ax = fig.gca()
        contour = ax.contour(self.xgrid, self.ygrid, self.phigrid, [value,])

        # set axis limits
        ax.set_xlim(self.xmin, self.xmax)
        ax.set_ylim(self.ymin, self.ymax)

        # set axis labels
        plt.xlabel('X')
        plt.ylabel('Y')

        # Set axes to be same length
        plt.gca().set_aspect('equal')

        # plot eye point
        ax.plot(self.eye_x, self.eye_y, 'ro')

        plt.savefig(path)

        plt.cla()
        plt.clf()
        plt.close()


    def add_circle(self, x, y, r):
        """
        Will add a circle to Phi.
        """
        newgrid = (self.xgrid - x)**2 + (self.ygrid - y)**2 - r**2
        self.phigrid = np.minimum(self.phigrid, newgrid)

    def add_ellipse(self, x, a, y, b):
        """
        Will add an elipse to Phi. The major axis will be aligned with
        either the x or y coordinate axes.
        """
        newgrid = ((self.xgrid - x)/a)**2 + ((self.ygrid - y)/b)**2 - 1
        self.phigrid = np.minimum(self.phigrid, newgrid)


    def add_rectangle(self, xmin, xmax, ymin, ymax):
        """
        Will add a rectangle to Phi. The rectangle will be aligned with
        the coordinate axes.
        """
        newgrid = self.xgrid - xmax

        newgrid1 = self.ygrid - ymax
        newgrid = np.maximum(newgrid, newgrid1)

        newgrid1 = xmin - self.xgrid
        newgrid = np.maximum(newgrid, newgrid1)

        newgrid1 = ymin - self.ygrid
        newgrid = np.maximum(newgrid, newgrid1)

        self.phigrid = np.minimum(self.phigrid, newgrid)

    def transport(self, k=0.01):
        """
        Will evolve figures according to the transport equation and the eye-point.

        args:
            k: A positive float indicating the timestep. If no timestep
                is provided then k=grid_step_size/2.
        """
        if k is None:
            k = self.grid_step/2

        C = k/self.grid_step

        # TODO: This is not right!!! 
        # newgrid = self.phigrid + C*(2*self.phigrid - self.phigrid - self.phigrid)
        # TODO: Try this instead
        newgrid = np.empty_like(self.phigrid)
        for i in range(newgrid.shape[0]-1):
            for j in range(newgrid.shape[1]-1):
                newgrid[i,j] = self.phigrid[i,j] + C*(
                    2*self.phigrid[i,j] - 
                    self.flowy[i,j] * self.phigrid[i,j+1] - 
                    self.flowx[i,j] * self.phigrid[i+1,j])

        self.phigrid = newgrid


def test():
    # TESTING
    # one circle
    phi_circ = Phi()
    phi_circ.add_circle(0,0,0.5)
    phi_circ.plot_phi('circle-1-3d.png')
    phi_circ.plot_level_set('circle-1-levelset.png')
    phi_circ.plot_level_set('circle-1-levelset-25.png', 0.25)

    # two circles
    phi_2circ = Phi(eye_x=0.5, eye_y=0.5)
    phi_2circ.add_circle(0.5, 0.5, 0.5)
    phi_2circ.add_circle(-0.5, -0.5, 0.25)
    phi_2circ.plot_phi('circle-2-3d.png')
    phi_2circ.plot_level_set('circle-2--levelset.png')

    # one rectangle
    phi_rect = Phi(eye_x=0, eye_y=-0.75)
    phi_rect.add_rectangle(-0.5, 0.5, -0.75, 0.75)
    phi_rect.plot_phi('rectangle-1-3d.png')
    phi_rect.plot_level_set('rectangle-1--levelset.png')

    # two rectangles
    phi_2rect = Phi()
    phi_2rect.add_rectangle(-0.75, 0, -0.75, 0)
    phi_2rect.add_rectangle(0, 0.25, 0, 0.5)
    phi_2rect.plot_phi('rectangle-2-3d.png')
    phi_2rect.plot_level_set('rectangle-2-levelset.png')
    phi_2rect.plot_level_set('rectangle-2-levelset-50.png', 0.5)

    # one ellipse
    phi_parab = Phi()
    phi_parab.add_ellipse(0,2,0,1)
    phi_parab.plot_phi('ellipse-long-x-3d.png')
    phi_parab.plot_level_set('ellipse-long-x-levelset.png')
    phi_parab2 = Phi()
    phi_parab2.add_ellipse(0.5,1,0,2)
    phi_parab2.plot_phi('ellipse-long-y-3d.png')
    phi_parab2.plot_level_set('ellipse-long-y-levelset.png')

    phit = Phi()
    phit.add_circle(0, 0.5, -0.1)
    phit.plot_level_set('test-transport-levelset-0.png')
    phit.plot_phi('test-transport-surface-0.png')
    for k in range(11):
        phit.transport()
        phit.plot_level_set('test-transport-levelset-{}.png'.format(k))
        phit.plot_phi('test-transport-surface-{}.png'.format(k))

    
def run_transport(eye_x, eye_y, plotname):
    # parameters for running and saving
    n_iter = 1000
    savepoints = [10,50,100,500,1000]

    # initialize phi
    phi = Phi(eye_x=eye_x, eye_y=eye_y)

    # circles
    phi.add_circle(-0.85, 0, 0.05)
    phi.add_circle(-0.25, 0, 0.1)
    phi.add_circle(0.5, 0.5, 0.15)
    # ellipses
    phi.add_ellipse(0.75, 0.1, 0.05, 0.15)
    phi.add_ellipse(0.25, 0.3, -0.5, 0.1)
    # rectangle
    phi.add_rectangle(-0.7, -0.5, -0.5, 0.5)

    # save figures of initial states:
    phi.plot_level_set('{}-levelset-0.png'.format(plotname))
    phi.plot_phi('{}-surface-0.png'.format(plotname))

    for it in range(1, n_iter+1):
        phi.transport()
        if it in savepoints:
            phi.plot_level_set('{}-levelset-{}.png'.format(plotname, it))
            phi.plot_phi('{}-surface-{}.png'.format(plotname, it))


if __name__ == '__main__':
    if TEST:
        test()

    else:
        run_transport(0, 0, 'PartA')
