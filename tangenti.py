import sympy as sp

x = sp.Symbol('x')
f = sp.tan(x)-x #funzione di cui trovare gli zeri
df = sp.diff(f, x) #derivata della funzione f
t = 1e-13 #tolleranza

def tangenti(x0, t):
    iter = 1
    while abs(f.subs(x, x0))>=t:
        x0 = x0 - ( f.subs(x,x0) / df.subs(x,x0) )
        iter += 1
        if iter > 10000 or abs(f.subs(x, x0))>500:
            if iter > 10000:
                raise Exception('troppe iterazioni')

            if abs(f.subs(x, x0))>500:
                raise Exception('la soluzione sta divergendo\nscegliere meglio il punto di partenza')

    return x0, iter


#valore iniziale da cui partire
init = 4.4


xs, iter = tangenti(init, t)

print(iter , " iterazioni necessarie")


print("xs= %.15f" %xs)

print("|f(xs)|= %e" %abs(f.subs(x,xs)))


