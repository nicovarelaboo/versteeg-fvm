from tdma import TDMA
import numpy as np
import matplotlib.pyplot as plt

k = 0.5
n = 5
L = 2e-2
deltax = L / n
A = 1
q = 1000e3
TA = 100
TB = 200

aW = k * A / deltax 
aE = k * A / deltax
aP_interior = aW + aE 
Sp = -2 * k * A / deltax 
Su = q * A * deltax
aP_A = aE - Sp
aP_B = aW - Sp 
Su_A = 2 * k * A / deltax * TA + q * A * deltax
Su_B = 2 * k * A / deltax * TB + q * A * deltax

A = np.zeros((n,n))
b = np.zeros(n)

A[0, 0] = aP_A

for i in range(1, len(A)-1):
    A[i, i] = aP_interior
    b[i] = Su

A[n-1, n-1] = aP_B 

for i in range(0, len(A)-1):
    A[i, i+1] = -aE

for i in range(1, len(A)):
    A[i, i-1] = -aW

b[0] = Su_A

b[-1] = Su_B

T = TDMA(A, b)

print(f'T = {T}')

def sol_analitica(x):
    return ((TB - TA) / L + q / 2 / k * (L - x)) * x + TA

x = np.linspace(0, L, 50)

plt.plot([deltax / 2 + deltax * i for i in range(n)], T, 'o') 
plt.plot(x, sol_analitica(x))
plt.grid()
plt.show()

