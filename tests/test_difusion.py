import numpy as np

from fvm.malla import crear_malla
from fvm.ensamblado import resolver_difusion
from fvm.contornos import Contorno
from fvm.fuentes import Fuente

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

def prueba_generacion_uniforme():

    m = crear_malla(L=0.02, n=5, area=1.0)

    oeste = Contorno(tipo="dirichlet", valor=100.0)
    este = Contorno(tipo="dirichlet", valor=200.0)

    fuente = Fuente(constante=1000e3, lineal=0.0)

    phi = resolver_difusion(malla=m, difusividad=0.5, oeste=oeste, este=este, fuentes=(fuente,))

    esperado = [150, 218, 254, 258, 230]

    return np.allclose(phi, esperado)

def prueba_aleta():

    m = crear_malla(L=1.0, n=10, area=1.0)

    oeste = Contorno(tipo="dirichlet", valor=100.0)
    este = Contorno(tipo="neumann", valor=0.0)

    fuente = Fuente(constante=25*20.0, lineal=-25.0)

    phi = resolver_difusion(malla=m, difusividad=1.0, oeste=oeste, este=este, fuentes=(fuente,))

    esperado = [80.59909958, 56.94707362, 42.53181608, 33.74951255, 28.40458716,
                25.16080856, 23.2072321, 22.05546367, 21.41756115, 21.13404892]

    return np.allclose(phi, esperado)


if __name__ == "__main__":
    print(f"Prueba flujo este: {prueba_flujo_este()}")
    print(f"Prueba flujo oeste: {prueba_flujo_oeste()}")
    print(f"Prueba aislado este: {prueba_aislamiento_este()}")
    print(f'Prueba generación uniforme: {prueba_generacion_uniforme()}')
    print(f'Prueba aleta: {prueba_aleta()}')

