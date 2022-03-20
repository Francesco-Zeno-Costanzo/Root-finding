import time
import numpy as np
import matplotlib.pyplot as plt


def f(x) :
    return 5.0+4.0*x-np.exp(x)

x1=0.0
x2=4.0
t=1.0e-13
x=np.linspace(x1, x2, 1000)
#plot per vedere come scegliere gli estremi
plt.figure(1)
plt.plot(x, f(x))
plt.grid()
plt.show()

##metodo secanti
start_time=time.time()

iter=1
while abs(x2-x1)>t:
    x3 = x2 - ((x2-x1)/(f(x2)-f(x1)))*f(x2)
    x1, x2 = x2, x3
    iter+=1

a=(time.time() - start_time)
print(iter , " iterazioni necessarie, che in termini di tempo sono:")
print("--- %s seconds ---" %a)
print("x0 = " ,x3)
print("accuracy = %.2e" %(x2-x1))
print("f(x0)=" ,f(x3))
