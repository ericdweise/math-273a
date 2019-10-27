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


def read_image(path):
    """
    Will read in the image and return a numpy array.

    Args:
        path: Path to the image file.

    Returns:
        A Numpy array representation of the file.
    """
    pass


def save_image(array, path):
    """
    Will save a 2D Numpy array to a file.

    Args:
        array: A 2D numpy array.
        path: The path where the image should be saved.
    """


def eno(array, n_samples, n_nodes):
    """
    perform ENO on a 1D array, using n points for each interpolation.
    Assumes Von Neumann boundary conditions.

    Args:
        array: A 1D array of values that represent the original data.
        n_samples: The number of samples to interpolate in the returned array.
        n_nodes: The number of original nodes to use in each interpolation.

    Returns:
        A 1D array with length "n_samples".
    """
    pass


def resample_image(array, n_x, n_y, n_nodes):
    """
    Will perform two iterations of ENO on a 2D array. The first will resample the original array in the "x" direction and the second will resample the array in the "y" direction.

    Args:
        array: A 2D Numpy array representing pixel information of an image.
        n_x: The number of samples in the x direction for the output image.
        n_y: The number of samples in the y direction for the output image.
        n_nodes: The number of original nodes to use in each interpolation.

    Returns:
        A 2D Numpy array that is the resampled image.
    """
    pass


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
    print('Part A has not yet been implemented')


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
    print('Part B has not yet been implemented')


def partC():
    """
    Perform the following steps:
    1) read the image dolphinorig.dat 
    2) resize the image using parameters 
        n_x = 300
        n_y = 100
        n_nodes = 10
    3) save the output image to eric-weise-part-3.dat
    """
    print('Part C has not yet been implemented')


if __name__ == '__main__':
    partA()
    partB()
    partC()
