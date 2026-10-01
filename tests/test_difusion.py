import numpy as np

from fvm.malla import crear_malla
from fvm.ensamblado import resolver_difusion
from fvm.contornos import Contorno

def prueba_flujo_este():

    m = crear_malla(L=0.5, n=5, area=0.01)

    oeste = Contorno(tipo="dirichlet", valor=200.0)
    este = Contorno(tipo="neumann", valor=100)

    phi = resolver_difusion(malla=m, difusividad=1000, oeste=oeste, este=este)

    esperado = 200.0 + 100.0 * m.nodos

    return np.allclose(phi, esperado)

def prueba_flujo_oeste():

    m = crear_malla(L=0.5, n=5, area=0.01)

    oeste = Contorno(tipo="neumann", valor=-800.0)
    este = Contorno(tipo="dirichlet", valor=500.0)

    phi = resolver_difusion(malla=m, difusividad=1000, oeste=oeste, este=este)

    esperado = 100.0 + 800 * m.nodos

    return np.allclose(phi, esperado)

def prueba_aislamiento_este():

    m = crear_malla(L=0.5, n=5, area=0.01)

    oeste = Contorno(tipo="dirichlet", valor=100.0)
    este = Contorno(tipo="neumann", valor=0.0)

    phi = resolver_difusion(malla=m, difusividad=1000, oeste=oeste, este=este)

    esperado = np.full(m.n, 100)

    return np.allclose(phi, esperado)


if __name__ == "__main__":
    print(f"Prueba flujo este: {prueba_flujo_este()}")
    print(f"Prueba flujo oeste: {prueba_flujo_oeste()}")
    print(f"Prueba aislado este: {prueba_aislamiento_este()}")

