import numpy as np

def coeficientes_cds(n, u, L, rho, Gamma, phiA, phiB):

    deltax = L / n

    F = rho * u
    D = Gamma / deltax
    
    Fe = Fw = FA = FB = F
    De = Dw = D
    DA = DB = 2 * D

    aE = np.full(n, De - Fe / 2) 
    aW = np.full(n, Dw + Fw / 2)
 
    aW[0] = 0

    aE[n-1] = 0

    Sp = np.zeros(n)
    Su = np.zeros(n)

    Sp[0] = -(DA + FA)
    Sp[n-1] = - (DB - FB)  

    Su[0] = (DA + FA) * phiA 
    Su[n-1] = (DB - FB) * phiB 

    aP = aW + aE + (Fe - Fw) - Sp 

    A = np.zeros((n,n))

    for i in range(n):
        
        A[i, i] = aP[i]
        if i-1 >= 0: A[i, i-1] = - aW[i]
        if i+1 < n: A[i, i+1] = - aE[i]

    b = Su

    return A, b

