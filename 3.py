import numpy as np
from math import sqrt

ITERATION_LIMIT = 1000
it = 0
A = np.array([[6.92, 1.28, 0.79, 1.15, -0.66],
              [0.92, 3.5, 1.3, -1.63, 1.02],
              [1.15, -2.46, 6.1, 2.1, 1.483],
              [1.14, -1.68, -1.217, 9, -3],
              [1.33, 0.16, 2.1, 5.44, -18]])

b = np.array([2.1, 0.72, 3.87, -1.08, 13.8])
x = np.zeros_like(b)

for it_count in range(ITERATION_LIMIT):
    print("\nIterations: ", it)
    print("Current solution:", x)
    error = b - np.dot(A, x)
    print("Вектор невязки:", error)
    print("Норма вектора: ", sqrt(sum(error[i]**2 for i in range(5))))
    x_new = np.zeros_like(x)

    for i in range(A.shape[0]):
        s1 = np.dot(A[i, :i], x[:i])
        s2 = np.dot(A[i, i + 1:], x[i + 1:])
        x_new[i] = (b[i] - s1 - s2) / A[i, i]

    if sqrt(sum(error[i]**2 for i in range(5))) < 1e-5:
        break

    x = x_new
    it += 1

print("\n\nSolution: ", x)
print("Вектор невязки: ",  b - np.dot(A, x))