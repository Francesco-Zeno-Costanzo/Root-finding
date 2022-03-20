import time
import numpy as np
import matplotlib.pyplot as plt


def f(x) :
    return 5.0+4.0*x-np.exp(x)

a=0.0
b=4.0
t=1.0e-15
x=np.linspace(a, b, 1000)
#plot per vedere come scegliere gli estremi
plt.figure(1)
plt.plot(x, f(x))
plt.grid()
plt.show()

##metodo bisezione
start_time=time.time()
fa=f(a)
fb=f(b)
if fa*fb>0:
    print("protrebbero esserci più soluzioni" , fa , fb)

iter=1
while (b-a)>t:
    c=(a+b)/2.0
    fc=f(c)
    if fc*fa>0:
        a=c
    else:
        b=c
    iter+=1

k=(time.time() - start_time)
print(iter , " iterazioni necessarie, che in termini di tempo sono:")
print("--- %s seconds ---" %k)
print("x0 = " ,c)
print("accuracy = " , '{:.2e}' .format(b-a))
print("f (x0)=" ,f(c))

