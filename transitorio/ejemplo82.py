import numpy as np 
from tdma import TDMA

def conduccion_implicita_1D(n=5, deltat=2, tfinal=(40,)):

    L = 0.02
    k = 10.0
    rhoc = 10e6
    T0 = 200.0
    TB = 0
    deltax = L / n

    D = k / deltax 

    aW = np.full(n, D)
    aE = np.full(n, D)
    aW[0] = 0
    aE[n-1] = 0
    aP0 = np.full(n, rhoc * deltax / deltat)
    Sp = np.zeros(n)
    Sp[n-1] = - 2 * D
    Su = np.zeros(n)
    Su[n-1] = 2 * D * TB
    
    aP = aP0 + aW + aE - Sp  

    i = np.arange(n)

    A = np.zeros((n,n))
    A[i,i] = aP 
    A[i[1:n],i[0:n-1]] = - aW[1:]
    A[i[0:n-1],i[1:n]] = - aE[:-1]

    resultado = []

    for tiempo in tfinal:

        T = np.full(n, T0)
        index = 0
        max_iter = round(tiempo / deltat)   

        while index < max_iter:

            b = Su + aP0 * T

            T = TDMA(A, b)

            index += 1

        resultado.append(T)

    return resultado 

if __name__ == "__main__":
    resultado = conduccion_implicita_1D()
    print(resultado[0])




