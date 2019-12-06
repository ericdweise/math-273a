import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.mplot3d import axes3d, Axes3D


def plot_surface(xgrid, ygrid, graph, path):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(xgrid, ygrid, graph)

    # Set axis limits
    ax.set_xlim(np.min(xgrid), np.max(xgrid))
    ax.set_ylim(np.min(ygrid), np.max(ygrid))
    ax.set_zlim(np.min(graph), np.max(graph))

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


def plot_heatmap(xgrid, ygrid, graph, path):
    fig = plt.figure()
    ax = fig.gca()
    extent = [np.min(xgrid), np.max(xgrid), np.min(ygrid), np.max(ygrid)]
    plt.imshow(graph, cmap='cividis', interpolation='nearest', extent=extent)

    # Set number of ticks along axes
    plt.locator_params(axis='x', nbins=3)
    plt.locator_params(axis='y', nbins=3)

    # set axis labels
    plt.xlabel('X')
    plt.ylabel('Y')

    # plt.show()
    plt.savefig(path)

    plt.cla()
    plt.clf()
    plt.close()


def plot_level_set(xgrid, ygrid, graph, level, path):
    fig = plt.figure()
    ax = fig.gca()
    _ = plt.contour(self.xgrid, self.ygrid, self.phigrid, [level,])

    # set axis limits
    ax.set_xlim(np.min(xgrid), np.max(xgrid))
    ax.set_ylim(np.min(ygrid), np.max(ygrid))

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
