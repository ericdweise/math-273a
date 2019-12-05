import numpy as np

VERBOSE = True


def conjugate_gradient(A, b, threshold):
    RESIDUAL_CHECK = 20
    basis = standard_basis(A.shape[0])
    conjugate_basis = []
    x = np.zeros(A.shape[0])

    N = A.shape[0]
    n = 0

    for basis_vector in basis:
        n += 1
        p = find_anti_projection(A, basis_vector, conjugate_basis)
        conjugate_basis.append(p)

        x = x + compute_alpha(A, b, p)*p

        if n % RESIDUAL_CHECK == 0:
            residual = find_residual(A, b, x)
            if VERBOSE:
                print(f'\r  * Step {n} of {N} -- residual: {residual}')

            if residual < threshold:
                break

    return x


def compute_alpha(A, b, p):
    return np.matmul(p, b) / inner_product(A, p, p)


def standard_basis(n):
    return [np.eye(n)[i] for i in range(n)]


def inner_product(A, vleft, vright):
    return np.matmul(vleft, np.matmul(A, vright))


def find_residual(A, b, x):
    ORDER=2
    return np.linalg.norm(b - np.matmul(A,x), ORDER)/np.linalg.norm(b, ORDER)


def find_projection(A, v1, v2):
    """Find projection of v1 onto v2"""
    return (inner_product(A, v1, v2)/inner_product(A, v2, v2))*v2


def find_anti_projection(A, vector, on_set):
    for basis_vector in on_set:
        vector = vector - find_projection(A, vector, basis_vector)
    return vector
