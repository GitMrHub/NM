import numpy as np

it = 0
file = open("size_and_vector", 'r')
n = int(file.readline())
b = np.array(file.readline().split(), dtype=float)
a = np.array([list(map(float, j.strip().split())) for j in file])
b1 = np.eye(n, dtype=float)
file.close()

def gauss(x):
    x = np.array(x, float)
    print(x[1:] / x[0])
    return x[1:] / x[0]

def gauss_app(C, t):
    C = np.array(C, float)
    t = np.array([[t[i]] for i in range(len(t))], float)
    C[1:, :] = C[1:, :] - t * C[0, :]
    return C

def lu(A):
    LU = np.array(A, float)
    print(LU[0])
    for k in range(LU.shape[0] - 1):
        t = gauss(LU[k:, k])
        LU[k + 1:, k] = t
        LU[k:, k + 1:] = gauss_app(LU[k:, k + 1:], t)
    return LU

def solve_lu(A, b, flag):
    LU = lu(A)
    b = np.array(b, float)
    for i in range(1, len(b)):
        b[i] = b[i] - np.dot(LU[i, :i], b[:i])
        global it
        it += 1
    print()
    for i in range(len(b) - 1, -1, -1):
        b[i] = (b[i] - np.dot(LU[i, i + 1:], b[i + 1:])) / LU[i, i]
        if flag:
            print('Iteration res:', b)
        it += 1
    print('Iterations: ', it)
    return b

root = solve_lu(a, b, True)
print('\nThe root: ', root)
print('\nВектор невязки: ', b - np.dot(a, root), '\n')
res2 = solve_lu(a, b1, False)
print('Determinant: ', np.linalg.det(a))
print('\nReversed matrix\n', res2)
print('\nA * A^-1:\n', np.dot(a, res2))
