import numpy as np
from scipy.optimize import root

def newton(f, start, tol, args=()):
    '''
    Generalizzation of newton method for n equations
    
    Parameters
    ----------
    f : callable
        A vector function to find a root of.
    start : ndarray
        Initial guess.
    tol :float, optional, default 1e-8
        required tollerance
    args : tuple, optional
        Extra arguments passed to f
    
    Return
    ------
    x0 : 1darray
        solution of the system
    '''
    
    # initial guess
    x0 = start
    f0 = f(x0, *args)
    # for the computation of jacobian
    nd = len(x0)
    df = np.zeros((nd, nd))
    h  = 1e-8
    s  = np.zeros(nd)
    # count
    #n_iter = 0
    
    while True:
        
        # compute jacobian
        for i in range(nd):
            for j in range(nd):
                s[j] = 1
                xr, xl = x0 + h*s, x0 - h*s
                df[i, j] = (f(xr, *args) - f(xl, *args) )[i]/(2*h)
                s[:] = 0
        
        # update solution
        delta = np.linalg.solve(df, f0)
        x0 = x0 - delta
        f0 = f(x0, *args)
        #n_iter += 1        
        
        # stop condition
        if all(abs(f0) < tol) :
           break
           
    return x0
        
def system(V, a):
    x1, x2 = V
    
    r1 = x1**2 + x2**2 - a
    r2 = x2 - x1**2 + x1/2
    
    R = np.array([r1, r2])
    return R

init = np.array([0.2, 0])
tol = 1e-12

sol = newton(system, init, tol, args=(1, ))
print("Solution with newton: ", *sol)

start = (0.2, 0.0)
sol = root(system, start, method='hybr', args=(1, ))
print("Solution with root:   ", *sol.x)
