import numpy as np 
import matplotlib.pyplot as plt
from fvm.resolvedores import TDMA 

L = 1
n = 25
rho = 1
Gamma = 0.1
deltax = L / n
phiA = 1
phiB = 0

D = Gamma / deltax 

# ---------- Caso 2 ----------
u = 2.5
F = rho * u
A = np.zeros((n,n))
b = np.zeros(n)

aW = np.max([F, D + F/2, 0]) 
aE = np.max([-F, D - F/2, 0]) 

Sp1 = - (2 * D + F) 
Su1 = (2 * D + F) * phiA 

Sp5 = - 2 * D
Su5 = 2 * D * phiB 

A[0, 0] = aE - Sp1
b[0] = Su1

for i in range(1, n-1):
    A[i, i] = aW + aE 

for i in range(0, n-1):
    A[i, i+1] = - aE 

for i in range(1, n):
    A[i, i-1] = - aW 

b[n-1] = Su5

A[n-1, n-1] = aW - Sp5

phi2 = TDMA(A, b)

def solucion_analitica(x):
    return 1 + (1 - np.exp(25*x)) / 7.20e10

x = np.linspace(0,1,50)

plt.plot(x, solucion_analitica(x))
plt.plot([deltax / 2 + deltax * i for i in range(0,n)], phi2, 'o')
plt.grid()
plt.show()