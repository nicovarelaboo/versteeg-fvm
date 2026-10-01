import numpy as np

def coeficientes_quick(n, u, L, rho, Gamma, phiA, phiB):

    deltax = L / n

    F = rho * u
    D = Gamma / deltax
    
    Fe = Fw = FA = FB = F
    De = Dw = DA_star = DB_star = D

    aE = np.full(n, De - 3 / 8 * Fe) 
    aW = np.full(n, Dw + 6 / 8 * Fw + 1 / 8 * Fe)
    aWW = np.full(n, - 1 / 8 * Fw)

    aE[0] = De + 1 / 3 * DA_star - 3 / 8 * Fe 
    aW[0] = 0
    aWW[0] = 0

    aE[1] = De - 3 / 8 * Fe 
    aW[1] = Dw + 7 / 8 * Fw + 1 / 8 * Fe 
    aWW[1] = 0 

    aE[n-1] = 0
    aW[n-1] = Dw + 1 / 3 * DB_star + 6 / 8 * Fw 
    aWW[n-1] = - 1 / 8 * Fw 

    Sp = np.zeros(n)
    Su = np.zeros(n)

    Sp[0] = - (8 / 3 * DA_star + 2 / 8 * Fe + FA)
    Sp[1] = 1 / 4 * Fw 
    Sp[n-1] = - (8 / 3 * DB_star - FB)

    Su[0] = (8 / 3 * DA_star + 2 / 8 * Fe + FA) * phiA 
    Su[1] = - 1 / 4 * Fw * phiA 
    Su[n-1] = (8 / 3 * DB_star - FB) * phiB 

    aP = aW + aE + aWW + (Fe - Fw) - Sp 

    A = np.zeros((n,n))

    for i in range(n):
        
        A[i, i] = aP[i]
        if i-1 >= 0: A[i, i-1] = - aW[i]
        if i-2 >= 0: A[i, i-2] = - aWW[i]
        if i+1 < n: A[i, i+1] = - aE[i]

    b = Su

    return A, b

