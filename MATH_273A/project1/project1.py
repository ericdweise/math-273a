"""
Project 1 for Math 273A. 

To run this code use the following command:
    python3 project1.py

File requirements:
    Ensure that the file dolphinorig.dat exists in the same directory as this file.

Package Requirements:
    numpy: pip install numpy

Instructor: Li-Tien Cheng
Academic Quarter: Fall 2019
Author/student: Eric Weise
"""

import numpy as np
from math import floor, ceil
import cv2


def read_image(path):
    """
    Will read in the image and return a numpy array.
    Args:
        path: Path to the image file.
    Returns:
        A Numpy array representation of the file.
    """
    img = np.fromfile(path, dtype=float, sep=' ')

    with open(path, 'r') as f:
        line = f.readline()
        row = np.fromstring(line.strip(), dtype=float, sep=' ')
        N = row.shape[0]

    img = np.reshape(img, (-1, N))

    rotated_img = np.rot90(img)

    return rotated_img


def save_image(array, path, save_png=False):
    """
    Will save a 2D Numpy array to a file.

    Args:
        array: A 2D numpy array.
        path: The path where the image should be saved.
    """
    if save_png:
        a = np.array(array, copy=True)*255
        cv2.imwrite(path, a)
        # a = a*255
        # png.from_array(a, 'L8').save(path)
    else:
        print("Save dat file is not implemented")


def create_divided_difference_table(array, depth):
    table = np.zeros(shape=(array.shape[0], depth))
    
    for i in range(array.shape[0]):
        table[i,0] = array[i]

    for j in range(1, depth):
        for i in range(array.shape[0]-j):
            table[i,j] = (table[i+1,j-1] - table[i,j-1])/j

    return table
    


def interpolate(x, array, n_nodes, ddiffs):
    if round(x) == x:
        return(array[int(x)])

    interp = 0
    x_factor = 1

    choices = [int(floor(x)), int(ceil(x))]

    for k in range(n_nodes):
        idx_dn = min(choices)
        idx_up = max(choices)

        # choose index to use
        if idx_dn < 0:
            index = idx_up
            choices.append(index+1)
            D = ddiffs[index-k, k]
        elif idx_up > array.shape[0] - 1:
            index = idx_dn
            choices.append(index-1)
            D = ddiffs[index, k]
        elif abs(ddiffs[idx_dn ,k]) < abs(ddiffs[idx_up-k, k]):
            index = idx_dn
            choices.append(index-1)
            D = ddiffs[index, k]
        else:
            index = idx_up
            choices.append(index+1)
            D = ddiffs[index-k, k]

        # add to interp
        interp = interp + D*x_factor

        # prepare for next iteration
        x_factor = x_factor*(x - index)

    return interp


def eno(array, n_samples, n_nodes):
    """
    Perform ENO on a 1D array at n_out evenly spaced points. Each point
    interpolated will be found using n_nodes from array to perform the 
    ENO. Assumes Von Neumann boundary conditions.
    Args:
        array: A 1D array of values that represent the original data.
        n_out: The number of samples to interpolate. AKA, the length of the returned array.
        n_nodes: The number of original nodes to use in each interpolation.
    Returns:
        A 1D array with length "n_samples".
    """
    x = [float(i)*array.shape[0]/n_samples for i in range(n_samples)]
    new_array = np.zeros(shape=n_samples)
    ddiff_table = create_divided_difference_table(array, n_nodes)

    for i in range(len(x)):
        new_array[i] = interpolate(x[i], array, n_nodes, ddiff_table)

    return new_array


def resample_image(image, n_x, n_y, n_nodes):
    """
    Will perform two iterations of ENO on a 2D array. The first will resample the original array in the "x" direction and the second will resample the array in the "y" direction.
    Args:
        image: A 2D Numpy array representing pixel information of an image.
        n_x: The number of samples in the x direction for the output image.
        n_y: The number of samples in the y direction for the output image.
        n_nodes: The number of original nodes to use in each interpolation.
    Returns:
        A 2D Numpy array that is the resampled image.
    """
    int_image = np.zeros(shape=(image.shape[0], n_x))
    new_image = np.zeros(shape=(n_y, n_x))

    # Resample in x direction (change number of columns)
    for row_idx in range(image.shape[0]):
        row = image[row_idx,:]
        new_row = eno(row, n_x, n_nodes)
        for i in range(new_row.shape[0]):
            try:
                int_image[row_idx,i] = new_row[i]
            except:
                print('i: {}'.format(i))
                print('row_idx: {}'.format(row_idx))
                print('shape matrix: {}'.format(int_image.shape))
                print('shape row: {}'.format(new_row.shape))
                raise(Exception(''))

    # Resample in y direction (change number of rows)
    for col_idx in range(int_image.shape[1]):
        column = int_image[:,col_idx]
        new_column = eno(column, n_y, n_nodes)
        for j in range(new_column.shape[0]):
            try:
                new_image[j,col_idx] = new_column[j]
            except:
                print('j: {}'.format(j))
                print('col_idx: {}'.format(col_idx))
                print('shape matrix: {}'.format(new_image.shape))
                print('shape col: {}'.format(new_col.shape))
                raise(Exception(''))

    return new_image


#######################################################
## Methods to answer questions asked in the project  ##
#######################################################
def partA():
    """
    Perform the following steps:
    1) read the image dolphinorig.dat
    2) resize the image using parameters:
        n_x = 1000
        n_y = 1000
        n_nodes = 4
    3) save the output image to eric-weise-part-1.dat
    """
    print('=== Running Part A ===')
    image = read_image('dolphinorig.dat')
    new_img = resample_image(image, 1000, 1000, 4)
    save_image(new_img, 'results/part-a.png', save_png=True)


def partB():
    """
    Perform the following steps:
    1) read the image dolphinorig.dat
    2) resize the image using parameters:
        n_x = original_image.shape[0]
        n_y = original_image.shape[1]
        n_nodes = 4
    3) save the output image to eric-weise-part-2.dat
    """
    print('=== Running Part B ===')
    image = read_image('dolphinorig.dat')
    new_img = resample_image(image, image.shape[0], image.shape[1], 4)
    save_image(new_img, 'results/part-b.png', save_png=True)


def partC():
    """
    Perform the following steps:
    1) read the image dolphinorig.dat 
    2) resize the image using parameters 
        n_x = 300
        n_y = 100
        n_nodes = 10
    3) save the output image
    """
    print('=== Running Part C ===')
    image = read_image('dolphinorig.dat')
    new_img = resample_image(image, 300, 100, 10)
    save_image(new_img, 'results/part-c.png', save_png=True)


if __name__ == '__main__':
    i = read_image('dolphinorig.dat')
    save_image(i, 'results/part-0.png', save_png=True)
    partA()
    partB()
    partC()
