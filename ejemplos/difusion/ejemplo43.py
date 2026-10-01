from fvm.resolvedores import TDMA
import numpy as np
import matplotlib.pyplot as plt

n = 10
L = 1
deltax = L / n
TA = 100
Tinf = 20
n2 = 25

aW = 1 / deltax 
aE = 1 / deltax
Sp_interior = - n2 * deltax
aP_interior = aW + aE - Sp_interior
Su_interior = n2 * deltax * Tinf
Su_A = 2 * TA / deltax + n2 * deltax * Tinf
Su_B = n2 * deltax * Tinf
Sp_A = - n2 * deltax - 2 / deltax
Sp_B = - n2 * deltax
aP_A = aE - Sp_A
aP_B = aW - Sp_B

A = np.zeros((n,n))
b = np.zeros(n)

A[0, 0] = aP_A

for i in range(1, n-1):
    A[i, i] = aP_interior
    b[i] = Su_interior

A[n-1, n-1] = aP_B 

for i in range(0, n-1):
    A[i, i+1] = -aE

for i in range(1, n):
    A[i, i-1] = -aW

b[0] = Su_A
b[-1] = Su_B

T = TDMA(A, b)

print(f'T = {T}')

def sol_analitica(x):
    return (TA-Tinf) * np.cosh(np.sqrt(n2)*(L-x)) / np.cosh(np.sqrt(n2)*L) + Tinf

x = np.linspace(0, L, 50)

plt.plot([deltax / 2 + deltax * i for i in range(n)], T, 'o') 
plt.plot(x, sol_analitica(x))
plt.grid()
plt.show()

