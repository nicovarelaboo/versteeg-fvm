import numpy as np 

def TDMA(A, b):

    n = len(A)

    x = np.zeros(n)

    u = np.zeros(n)
    v = np.zeros(n)
    
    u[0] = b[0] / A[0, 0]
    v[0] = A[0, 1] / A[0, 0]

    u[n-1] = b[n-1] / A[n-1, n-1]
    v[n-1] = A[n-1, n-1-1] / A[n-1, n-1]

    for i in range(1, n-1):
        u[i] = (b[i] - A[i, i-1] * u[i-1]) / (A[i,i] - A[i, i-1] * v[i-1]) 
        v[i] = A[i, i+1] / (A[i, i] - A[i, i-1] * v[i-1])

    x[n-1] = (u[n-1] - v[n-1] * u[n-2]) / (1 - v[n-1] * v[n-2])

    for i in range(n-2, -1, -1):
        x[i] = u[i] - v[i] * x[i + 1]

    return x


def tdma_coeficientes(aW: np.ndarray, aP: np.ndarray, aE: np.ndarray, Su: np.ndarray) -> np.ndarray:

    n = len(aP)

    x = np.zeros(n)

    u = np.zeros(n)
    v = np.zeros(n)
    
    u[0] = Su[0] / aP[0]
    v[0] = - aE[0] / aP[0]

    u[-1] = Su[-1] / aP[-1]
    v[-1] = - aW[-1] / aP[-1]

    for i in range(1, n-1):
        u[i] = (Su[i] + aW[i] * u[i-1]) / (aP[i] + aW[i] * v[i-1]) 
        v[i] = - aE[i] / (aP[i] + aW[i] * v[i-1])

    x[-1] = (u[-1] - v[-1] * u[-2]) / (1 - v[-1] * v[-2])

    for i in range(n-2, -1, -1):
        x[i] = u[i] - v[i] * x[i + 1]

    return x


if __name__ == "__main__":
    aW = np.array([0.0, 100.0, 100.0, 100.0, 100.0])
    aE = np.array([100.0, 100.0, 100.0, 100.0, 0.0])
    aP = np.array([300.0, 200.0, 200.0, 200.0, 300.0])
    Su = np.array([20000.0, 0.0, 0.0, 0.0, 100000.0])

    A = np.diag(aP) - np.diag(aW[1:], -1) - np.diag(aE[:-1], 1)
    print(TDMA(A, Su))
    print(tdma_coeficientes(aW, aP, aE, Su))


