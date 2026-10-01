import numpy as np 
from fvm.resolvedores import TDMA

np.set_printoptions(precision=3, suppress=True)

k = 1000
q = 500e3
TA = 100
Delta = 0.1
Deltaz = 1e-2
Area = Deltaz * Delta
a = k * Area / Delta
Su123 = q * Area

T2 = np.zeros(4)
T3 = np.zeros(4)
index = 0

while index <= 100:

    # ---------- Línea 1 ----------

    A1 = np.zeros((4,4))
    b1 = np.zeros(4)

    # Nodo 1
    aP1 = 2 * a
    A1[0, 0] = aP1 
    b1[0] = Su123 + a * T2[0]

    # Nodos 2 y 3
    aP23 = 3 * a
    for i in range(1, 4-1):
        A1[i, i] = aP23
        b1[i] = Su123 + a * T2[i]

    # Nodo 4
    Sp4 = - 2 * a
    aP4 = 2 * a - Sp4 
    Su4 = q * Area + 2 * a * TA
    A1[4-1, 4-1] = aP4 
    b1[4-1] = Su4 + a * T2[3]

    for i in range(0, 4-1):
        A1[i, i+1] = -a 

    for i in range(1, 4):
        A1[i, i-1] = -a

    T1 = TDMA(A1, b1)


    # ---------- Línea 2 ----------

    A2 = np.zeros((4,4))
    b2 = np.zeros(4)

    # Nodo 5
    aP5 = 3 * a
    Su5 = 0
    A2[0, 0] = aP5
    b2[0] = Su5 + a * T1[0] + a * T3[0]

    # Nodos 6 y 7
    aP67 = 4 * a
    for i in range(1, 4-1):
        A2[i, i] = aP67
        b2[i] = a * T1[i] + a * T3[i] 

    # Nodo 8
    Sp8 = - 2 * a
    aP8 = 3 * a - Sp8 
    Su8 = 2 * a * TA + a * T1[4-1] + a * T3[4-1]
    A2[4-1, 4-1] = aP8 
    b2[4-1] = Su8

    for i in range(0, 4-1):
        A2[i, i+1] = -a 

    for i in range(1, 4):
        A2[i, i-1] = -a

    T2 = TDMA(A2, b2)

    print(f'T2 = {T2}')


    # ---------- Línea 3 ----------

    A3 = np.zeros((4,4))
    b3 = np.zeros(4)

    # Nodo 9
    aP9 = 2 * a
    Su9 = 0
    A3[0, 0] = aP9
    b3[0] = Su9 + a * T2[0]

    # Nodos 10 y 11
    aP1011 = 3 * a
    for i in range(1, 4-1):
        A3[i, i] = aP1011
        b3[i] = a * T2[i] 

    # Nodo 12
    Sp12 = - 2 * a
    aP12 = 2 * a - Sp12 
    Su12 = 2 * a * TA + a * T2[4-1]
    A3[4-1, 4-1] = aP12 
    b3[4-1] = Su12

    for i in range(0, 4-1):
        A3[i, i+1] = -a 

    for i in range(1, 4):
        A3[i, i-1] = -a

    T3 = TDMA(A3, b3)

    index += 1

print(f'T1 = {T1}')
print(f'T2 = {T2}')
print(f'T3 = {T3}')
print(index)

