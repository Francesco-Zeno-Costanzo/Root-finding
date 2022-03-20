import time
import sympy as sp

x = sp.Symbol('x')
f = sp.tan(x)-x
df = sp.diff(f, x)
t = 1e-13

def tangenti(x0, t):
    iter=1
    while abs(f.subs(x, x0))>=t:
        x0 = x0 - ( f.subs(x,x0) / df.subs(x,x0) )
        iter+=1
        if iter > 10000 and abs(f.subs(x, x0))>500:
            if iter > 10000:
                print('troppe iterazioni')
                break
            if abs(f.subs(x, x0))>500:
                print('la soluzione sta divergendo\nscegliere meglio il punto di partenza')
                break
    return x0, iter


init = 4.4
start_time=time.time()

xs, iter = tangenti(init, t)

a=(time.time() - start_time)
print(iter , " iterazioni necessarie, che in termini di tempo sono:")
print("--- %s seconds ---" %a)

print("xs= %.15f" %xs)
print("c'ha senso?")

if abs(f.subs(x,xs)) <= t:
    print("si per quello che hai chiesto di fare, |f(xs)|= %e" %abs(f.subs(x,xs)))
else:
    print("pare di no:  |f(xs)|= %e" %abs(f.subs(x,xs)))
    print('cambia il valore di init')

