import numpy as np 
from fvm.resolvedores import TDMA

def conduccion_explicita_1D(n=5, deltat=2, tfinal=(10,)):

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
    aP = aP0 = np.full(n, rhoc * deltax / deltat)


    resultado = []

    for tiempo in tfinal:

        T = np.full(n, T0)
        index = 0
        max_iter = round(tiempo / deltat)   

        while index < max_iter:

            Su = np.zeros(n)
            Su[n-1] = 2 * D * TB

            Sp = np.zeros(n)
            Sp[n-1] = - 2 * D 

            i = np.arange(n)

            A = np.zeros((n,n))
            A[i,i] = aP 

            b = np.zeros(n)

            b[0] = aE[0] * T[1] + (aP0[0] - aE[0]) * T[0]

            b[1:n-1] = (aW[1:n-1] * T[0:n-2] + aE[1:n-1] * T[2:n] 
                        + (aP0[1:n-1] - (aW[1:n-1] + aE[1:n-1] - Sp[1:n-1])) * T[1:n-1])

            b[n-1] = aW[n-1] * T[n-2] + (aP0[n-1] - (aW[n-1] - Sp[n-1])) * T[n-1] + Su[n-1]

            T = TDMA(A, b)

            index += 1

        resultado.append(T)

    return resultado 


if __name__ == "__main__":
    resultado = conduccion_explicita_1D()
    print(resultado[0])

