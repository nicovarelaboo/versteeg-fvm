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

