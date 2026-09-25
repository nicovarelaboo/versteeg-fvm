import numpy as np
from tdma import TDMA

L = 1
n = 4
deltax = L / (n + 1/2)
AA = 0.5
AB = 0.1
pB = 1.0
uA = 1.0
rhoA = 1.0
rhoB = 1.0 
muA = 1.0
Su = 0

i = np.arange(n)

xu = deltax + deltax * i
xP = deltax / 2 + deltax * i

Au = AA - (AA - AB) / L * xu 
AP = AA - (AA - AB) / L * xP

p = np.full(n, pB)
u = np.full(n, uA)

rho = np.full(n, rhoA)
mu = np.full(n, muA)

FW = np.zeros(n)
FE = np.zeros(n)

Sp = np.zeros(n)

b = np.zeros(n)

DW = mu * AP / deltax 
DE = np.zeros(n)
DE[0:n-1] = mu[1:n] * AP[1:n] / deltax
DE[n-1] = 0 #provisionaL

for k in range(n):

    if k == 0:
        FW[k] = rho[k] * AP[k] * (u[k] + uA) / 2 
    else:
        FW[k] = rho[k] * AP[k] * (u[k] + u[k-1]) / 2
        
    if k == n-1:
        r =  (rhoB + rho[k]) / 2 * Au[k] / rhoB / AB
        FB = rhoB * AB * r * u[k] 
        Sp[k] = - FB * r
    else: 
        FE[k] = rho[k+1] * AP[k+1] * (u[k+1] + u[k]) / 2

    if k == n-1:
        b[k] = -(pB - p[k]) * Au[k] + Su *Au[k] * deltax
    elif k == 0:
        b[k] = (1/2 * FW[0] + DW[0]) * uA - (p[1] - p[0]) * Au[0] + Su * Au[0] * deltax
    else:
        b[k] = -(p[k+1] - p[k]) * Au[k] + Su * Au[k] * deltax
    
Sp[0] = - (DW[0] + 1/2 * FW[0])

aW = DW + 1 / 2 * FW
aW[0] = 0
aE = DE - 1 / 2 * FE 
aE[n-1] = 0
aP = aW + aE + (FE - FW) - Sp 

A = np.zeros((n,n))

for k in range(n):
    A[k,k] = aP[k]
    if k-1 >= 0: A[k, k-1] = - aW[k]
    if k+1 < n: A[k, k+1] = - aE[k]

u_star = TDMA(A, b)

