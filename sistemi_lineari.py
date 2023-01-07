import time
import numpy as np
import matplotlib.pyplot as plt


def kaczmarz(A, b, tol=1e-6, dense_output=False):
    """
    Aterative algorithm for solving
    linear equation systems

    Parameters
    ----------
    A : 2darray
        matrix of system.
    b : 1darray
        Ordinate or “dependent variable” values.
    tol : float, optional
        required tollerance default 1e-6
    dense_output : bool, optional
        if True all iteration of th solution, error and number
        of iteration are stored and reurned, default is False

    Return
    ------
    x : 1darray
        solution of system
    err : float
        error of solution

    if dense_output :
    s : 2darray
        all iteration of solution
    e : 1darray
        error of all iteration,
    k : int
        number of iteration
    """

    m, n = A.shape
    x = np.zeros(n)
    s = []
    e = []
    k = 0
    p = range(m)

    while True:

        i  = k % m
        ai = A[i,:]
        x = x + (b[i] - ai @ x)/np.sqrt(sum(ai**2))**2 * ai.conj()
        err = np.max(np.abs(A @ x - b))

        if dense_output:
            s.append(x)
            e.append(err)

        if err < tol:
            break
        k += 1

    if not dense_output:
        return x, err
    else:
        return np.array(s), np.array(e), k

if __name__ == "__main__":

    #creo una matrice e un termine noto
    N = 20
    A = np.zeros((N, N))
    b = np.zeros(N)
    for i in range(N):
        for j in range(N):
            A[i, j] = (np.random.random()-0.5)*20
        b[i] = (np.random.random()-0.5)*20


    t0 = time.time()
    sol, err, iter = kaczmarz(A, b, tol=1e-8, dense_output=True)
    dt = time.time() - t0
    x1 = sol[-1]
    print(f'tempo impiegato da kaczmarz: {dt} s')
    print(f'numero di iterazioni: {iter}')
    print('soluzione con kaczmarz\n', x1)

    #soluzione con numpy
    x2 = np.linalg.solve(A, b)
    print('soluzione con numpy\n', x2)

    #distanza delle due soluzioni
    print(f"distanza delle due soluzioni = {np.sqrt(np.sum((x1-x2)**2))}")

    plt.figure(1)
    plt.grid()
    plt.plot(err)
    plt.xlabel('iteration')
    plt.ylabel('error')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()