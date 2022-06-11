import time
import numpy as np
import matplotlib.pyplot as plt

N = 20
A = np.zeros((20, 20))
b = np.zeros(20)
for i in range(20):
    for j in range(20):
        A[i, j] = (np.random.random()-0.5)*20
    b[i] = (np.random.random()-0.5)*20



def kaczmarz(A, b, tol=1e-6):

    m, n = A.shape
    X = np.zeros(n)
    Y = [X]
    dY = [np.max(np.abs(A @ X - b))]
    k = 0

    p = range(m)


    while True:
        i = k % m
        ai = A[i,:]
        X_n = X + (b[i] - ai @ X) / np.linalg.norm(ai)**2 * ai
        err = np.max(np.abs(A @ X_n - b))
        Y.append(X_n)
        dY.append(err)
        if err < tol:
            break
        X = X_n
        k += 1
    return X, Y, dY


t0 = time.time()
X0, Y0, dY0 = kaczmarz(A, b, tol=1e-6)
dt = time.time() - t0
print(f'tempo impiegato: {dt}')
print(X0)

X1 = np.linalg.solve(A, b)
print(X1)

print(np.sqrt(np.sum((X0-X1)**2)))



