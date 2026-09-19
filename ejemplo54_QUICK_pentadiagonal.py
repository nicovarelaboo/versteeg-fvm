import numpy as np 
from coeficientes_quick import coeficientes_quick

u = 0.2
L = 1
n = 5
rho = 1
Gamma = 0.1
phiA = 1
phiB = 0

A, b = coeficientes_quick(n, u, L, rho, Gamma, phiA, phiB)

phi_quick = np.linalg.solve(A,b)

print(phi_quick)
