import numpy as np
from fvm.malla import Malla

def coeficientes_difusion(malla: Malla, difusividad: float | np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    # difusividad en las n+1 caras
    difusividad = np.broadcast_to(difusividad, malla.n + 1)

    # conductancias de las caras interiores (n-1 valores)
    conductancia_difusiva = difusividad[1:-1] * malla.area[1:-1] / malla.delta_nodos

    # aW: n ceros, y rellenar aW[1:] con las conductancias
    aW = np.zeros(malla.n)
    aW[1:] = conductancia_difusiva

    # aE: n ceros, y rellenar aE[:-1] con las conductancias 
    aE = np.zeros(malla.n)
    aE[:-1] = conductancia_difusiva

    return aW, aE 

if __name__ == "__main__":
    from fvm.malla import crear_malla

    m = crear_malla(0.5, 5, 0.01)
    aW, aE = coeficientes_difusion(m, 1000.0)
    print(aW)
    print(aE)

