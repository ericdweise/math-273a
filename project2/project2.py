"""
Project 2 for Math 273A: Level set functions.

To run this code use the following command:
    python3 project2.py

Package Requirements:
    numpy: pip install numpy

Instructor: Li-Tien Cheng
Academic Quarter: Fall 2019
Author/student: Eric Weise
"""

import numpy as np
from math import floor
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


class Phi(object):

    # Resolution of the grid
    grid_step_size = 0.02

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
        self.x_list = np.arange(xmin, xmax, grid_step)
        self.y_list = np.arange(ymin, ymax, grid_step)

        # 2D grid to hold x and y values
        self.xgrid, self.ygrid = np.meshgrid(self.x_list, self.y_list)

        # 2D grid to hold Phi values
        # Initialize with Inf
        self.phigrid = np.full(self.xgrid.shape, np.inf)


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
        plt.locator_params(axis='z', nbins=2)

        # set axis labels
        plt.xlabel('X')
        plt.ylabel('Y')

        plt.savefig(path)


    def plot_level_set(self, value=0):
        """Create a 2D plot of Phi(x,y)=value.
        Args:
            value: The value at which the level set is comuputed.
        """
        pass

    def add_circle(self, x, y, r):
        """
        Will add a circle to Phi.
        """
        newgrid = (self.xgrid - x)**2 + (self.ygrid - y)**2 - r**2
        self.phigrid = np.minimum(self.phigrid, newgrid)

    def add_ellipse(self, x, y, r):
        """
        Will add an elipse to Phi. The major axis will be aligned with
        either the x or y coordinate axes.
        """
        pass

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
            k = self.grid_step_size/2
        pass


if __name__ == '__main__':
    phi = Phi()
    phi.add_rectangle(-1, 0, -1, 0)
    phi.add_rectangle(0, 1, 0, 1)
    phi.plot_phi('two-rectangels.png')
