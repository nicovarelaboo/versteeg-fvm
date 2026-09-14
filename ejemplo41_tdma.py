from tdma import TDMA
import numpy as np
import matplotlib.pyplot as plt

k = 1000
deltax = 0.1
A = 10e-3
TA = 100
TB = 500

aW = k * A / deltax 
aE = k * A / deltax
aP_interior = aW + aE 
Sp = -2 * k * A / deltax 
aP_A = aE - Sp
aP_B = aW - Sp 
Su_A = 2 * k * A / deltax * TA
Su_B = 2 * k * A / deltax * TB

A = np.zeros((5,5))
b = np.zeros(5)

A[0, 0] = aP_A

for i in range(1, len(A)-1): A[i, i] = aP_interior

A[4, 4] = aP_B 

for i in range(0, len(A)-1): A[i, i+1] = -aE

for i in range(1, len(A)): A[i, i-1] = -aW

b[0] = Su_A

b[-1] = Su_B

T = TDMA(A, b)

print(f'T = {T}')

# Posiciones de los nodos (centro de cada volumen de control)
x_nodos = np.array([deltax/2 + i*deltax for i in range(5)])

# Solución analítica (recta entre T_A y T_B)
x_analitica = np.array([0, 0.5])
T_analitica = np.array([TA, TB])

plt.plot(x_analitica, T_analitica, '--', label='Solución analítica')
plt.plot(x_nodos, T, 'o', label='TDMA (volúmenes finitos)')
plt.xlabel('x (m)')
plt.ylabel('Temperatura (°C)')
plt.title('Ejemplo 4.1 — Difusión 1D en estado estacionario')
plt.legend()
plt.grid(True)
plt.show()