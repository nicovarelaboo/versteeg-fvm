import numpy as np 
from ejemplos.conveccion.coeficientes_cds import coeficientes_cds
from ejemplos.conveccion.coeficientes_quick import coeficientes_quick
import matplotlib.pyplot as plt

u = 0.2
L = 1
n = 5
rho = 1
Gamma = 0.1
phiA = 1
phiB = 0

deltax = L / n

A_CDS, b_CDS = coeficientes_cds(n, u, L, rho, Gamma, phiA, phiB)

A_quick, b_quick = coeficientes_quick(n, u, L, rho, Gamma, phiA, phiB)

phi_CDS = np.linalg.solve(A_CDS, b_CDS)

phi_quick = np.linalg.solve(A_quick, b_quick)

Pe_dominio = rho * u * L / Gamma 

def solucion_analitica(x):
    return phiA + (phiB - phiA) * (np.exp(Pe_dominio * x / L) - 1) / (np.exp(Pe_dominio) - 1)

xi = [deltax / 2 + deltax * i for i in range(0, n)]
x = np.linspace(0, L, 50)
y = solucion_analitica(x)

plt.figure()

plt.subplot(1, 2, 1)
plt.plot(x, y)
plt.plot(xi, phi_CDS, 'o')
plt.grid()

plt.subplot(1, 2, 2)
plt.plot(x, y)
plt.plot(xi, phi_quick, 'o')
plt.grid()
plt.show()


