import math
import numpy as np
import sys

np.set_printoptions(threshold=sys.maxsize)


def makeMatrixA(n):
    A = np.zeros((n,n))
    for i in range(n):
        A[i,i] = -2
    for i in range(n-1):
        A[i,i+1] = 1
        A[i+1,i] = 1
    return A

def makeMatricesMandN(n):
    M = np.zeros((n,n))
    N = np.zeros((n,n))
    for i in range(n):
        M[i,i] = -2
    for i in range(n-1):
        N[i,i+1] = 1
        N[i+1,i] = 1
    return M,N

def makeMatricesLDU(n):
    L = np.zeros((n,n))
    D = np.zeros((n,n))
    U = np.zeros((n,n))
    for i in range(n):
        D[i,i] = -2
    for i in range(n-1):
        U[i,i+1] = 1
        L[i+1,i] = 1
    return L,D,U


def findIterationMatrixJacobi(n):
    M,N = makeMatricesMandN(n)
    T = np.matmul(np.linalg.inv(M),-1*N)
    return T

def findIterationMatrixGaussSidel(n):
    L,D,U = makeMatricesLDU(n)
    T = np.matmul(np.linalg.inv(D+L),-1*U)
    return T

def findIterationMatrixSOR(n, w):
    L,D,U = makeMatricesLDU(n)
    T = np.matmul(np.linalg.inv(D + w*L), (1-w)*D - w*U)
    return T


def l2norm(m):
    return np.linalg.norm(m,2)

def test_outputs():
    print('===TESTING A creation for n=5 ===')
    A = makeMatrixA(5)
    print('\n--- A ---')
    print(A)

    print('\n\n=== TESTING M and N creation for n=5 ===')
    M,N = makeMatricesMandN(5)
    print('\n--- M ---')
    print(M)
    print('\n--- N ---')
    print(N)
    test_matrix = A==M+N
    assert(test_matrix.all())

    print('\n\n=== TESTING L,D,U CREATION for n = 5 ===')
    L,D,U = makeMatricesLDU(5)
    print('\n--- L ---')
    print(L)
    print('\n--- D ---')
    print(D)
    print('\n--- U ---')
    print(U)
    test_matrix = A==L+D+U
    assert(test_matrix.all())

    print('\n\n=== TESTING JACOBI ITERATION MATRIX ===')
    T = findIterationMatrixJacobi(5)
    print(T)

    print('\n\n=== TESTING Gauss Sidel ITERATION MATRIX ===')
    T = findIterationMatrixGaussSidel(5)
    print(T)

    print('\n\n=== TESTING SOR ITERATION MATRIX ===')
    T = findIterationMatrixSOR(5, 1.5)
    print(T)


def problem3():
    print('\n\n=== Problem 3 ===')
    print('--- Part b ---')
    for n in [100, 200, 400]:
        l2 = l2norm(findIterationMatrixJacobi(n))
        print('  n={} - l2={}'.format(n, l2))
    print('--- Part c ---')
    print('  n = {}'.format(math.log(0.00001,l2norm(findIterationMatrixJacobi(100)))))

def problem4():
    print('\n\n=== Problem 4 ===')
    print('--- Part b ---')
    for n in [100, 200, 400]:
        l2 = l2norm(findIterationMatrixGaussSidel(n))
        print('  n={} - l2={}'.format(n, l2))
    print('--- Part c ---')
    print('  n = {}'.format(math.log(0.00001,l2norm(findIterationMatrixGaussSidel(100)))))

def problem5():
    print('\n\n=== Problem 5 ===')
    print('--- Part b ---')
    for n in [100, 200, 400]:
        l2 = l2norm(findIterationMatrixSOR(n, 1.5))
        print('  n={} - l2={}'.format(n, l2))
    print('--- Part c ---')
    print('  n = {}'.format(math.log(0.00001,l2norm(findIterationMatrixSOR(100, 1.5)))))

def problem6():
    print('\n\n=== Problem 6 ===')
    B =  makeMatrixA(10) -2*np.identity(10)
    I = np.identity(10)
    A0 = np.append(np.append(B, I, 1), np.zeros((10, 80)), 1)
    A1 = np.append(np.append(np.append(I, B, 1), I, 1), np.zeros((10, 70)), 1)
    A2 = np.append(np.append(np.append(np.append(np.zeros((10, 10)), I, 1), B, 1), I, 1), np.zeros((10, 60)), 1)
    A3 = np.append(np.append(np.append(np.append(np.zeros((10, 20)), I, 1), B, 1), I, 1), np.zeros((10, 50)), 1)
    A4 = np.append(np.append(np.append(np.append(np.zeros((10, 30)), I, 1), B, 1), I, 1), np.zeros((10, 40)), 1)
    A5 = np.append(np.append(np.append(np.append(np.zeros((10, 40)), I, 1), B, 1), I, 1), np.zeros((10, 30)), 1)
    A6 = np.append(np.append(np.append(np.append(np.zeros((10, 50)), I, 1), B, 1), I, 1), np.zeros((10, 20)), 1)
    A7 = np.append(np.append(np.append(np.append(np.zeros((10, 60)), I, 1), B, 1), I, 1), np.zeros((10, 10)), 1)
    A8 = np.append(np.append(np.append(np.zeros((10, 70)), I, 1), B, 1), I, 1)
    A9 = np.append(np.append(np.zeros((10, 80)), I, 1), B, 1)
    A = np.append(np.append(np.append(np.append(np.append(np.append(np.append(np.append(np.append(A0,A1,0),A2,0),A3,0),A4,0),A5,0),A6,0),A7,0),A8,0),A9,0)
    with open('problem6-matrix.txt', 'w') as fout:
        fout.write('{}'.format(A))
    M = np.identity(100)*-4
    N = A - M
    T = -np.matmul(np.linalg.inv(M), N)
    print('--- Part b ---')
    for n in [100, 200, 400]:
        l2 = l2norm(T)
        print('  n={} - l2={}'.format(n, l2))
    print('--- Part c ---')
    print('  n = {}'.format(math.log(0.00001,l2norm(T))))

test_outputs()
problem3()
problem4()
problem5()
problem6()
