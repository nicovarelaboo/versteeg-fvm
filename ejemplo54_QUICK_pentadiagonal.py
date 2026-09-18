import numpy as np 
import matplotlib.pyplot as plt
from tdma import TDMA 

L = 1
n = 5
rho = 1
Gamma = 0.1
deltax = L / n
phiA = 1
phiB = 0

D = Gamma / deltax 

# ---------- CDS ----------
u = 0.2
F = rho * u

A = np.zeros((n,n))
b = np.zeros(n)

aW = D + F / 2
aE = D - F / 2

Sp1 = - (2 * D + F) 
Su1 = (2 * D + F) * phiA 

Sp5 = - (2 * D - F)
Su5 = (2 * D - F) * phiB 

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

phi1 = TDMA(A, b)

plt.figure()

def solucion_analitica(x):
    return 1 + (1 - np.exp(2*x)) / (np.exp(2) - 1)

x = np.linspace(0,1,50)

plt.subplot(1,2,1)
plt.plot(x, solucion_analitica(x))
plt.plot([deltax / 2 + deltax * i for i in range(0,n)], phi1, 'o')
plt.grid()


# ---------- QUICK ----------

u = 0.2
F = rho * u

A = np.zeros((n,n))
b = np.zeros(n)

# Nodo 1

aE1 = D + 1 / 3 * D - 3 / 8 * F 
aW1 = 0 
Sp1 = -( 8 / 3 * D + 2 / 8 * F + F)
Su1 = (8 / 3 * D + 2 / 8 * F + F) * phiA
aP1 = aE1 + aW1 - Sp1 
aWW1 = 0

A[0, 0] = aP1 
A[0, 1] = - aE1 
b[0] = Su1

# Nodo 2

aE2 = - 3 / 8 * F + D
aW2 = 1 / 8 * F + 7 / 8 * F + D 
Sp2 = 1 / 4 * F
Su2 = - 1 / 4 * F * phiA 
aP2 = aE2 + aW2 - Sp2 
aWW2 = 0

A[1,1] = aP2 
A[1,0] = -aW2 
A[1,2] = -aE2 
b[1] = Su2 

# Nodos interiores

aEi = D - 3 / 8 * F 
aWi = D + 6 / 8 * F  + 1 / 8 * F 
aWWi = - 1 / 8 * F
Spi = 0
Sui = 0
aPi = aWWi + aEi + aWi - Spi 

for i in range(2, n-1):
    A[i, i] = aPi  

for i in range(2, n-1):
    A[i, i+1] = - aEi 

for i in range(2, n-1):
    A[i, i-1] = - aWi

for i in range(2, n-1):
    A[i, i-1-1] = -aWWi

# Nodo n

aEn = 0
aWn = D + 1 / 3 * D + 6 / 8 * F 
aWWn = - 1 / 8 * F
Spn = - (8 / 3 * D - F)
Sun = (8 / 3 * D - F) * phiB
aPn = aEn + aWn + aWWn - Spn

A[n-1, n-1] = aPn
A[n-1, n-1-2] = - aWWn
A[n-1, n-1-1] = - aWn
b[n-1] = Sun

phi2 = np.linalg.solve(A, b)

plt.subplot(1,2,2)
plt.plot(x, solucion_analitica(x))
plt.plot([deltax / 2 + deltax * i for i in range(0,n)], phi2, 'o')
plt.grid()
plt.show()

