import numpy as np

from fvm1d.contornos import coeficientes_contorno, Contorno
from fvm1d.difusion import coeficientes_difusion
from tdma import tdma_coeficientes
from fvm1d.malla import Malla 


def resolver_difusion(malla: Malla, difusividad: float | np.ndarray,
                      oeste: Contorno, este: Contorno) -> np.ndarray:
    
    """
    Resuelve la difusión estacionaria 1D de φ con volúmenes finitos.

    Ecuación: d/dx(Γ dφ/dx) = 0, discretizada como
    aP φP = aW φW + aE φE + Su.

    Parameters
    ----------
    malla : Malla
        Malla 1D centrada en celdas.
    difusividad : float o ndarray
        Γ, un valor constante para todo el dominio o uno por cara (n+1 valores).
    oeste, este : Contorno
        Condiciones de contorno en x = 0 y x = L.

    Returns
    -------
    ndarray
        φ en los n nodos.
    """
    
    # 1. aW y aE de las caras interiores (coeficientes_difusion)
    aW, aE = coeficientes_difusion(malla, difusividad)

    # 1. Su y Sp: arrays de n ceros
    Su = np.zeros(malla.n)
    Sp = np.zeros(malla.n)

    # 3. difusividad en las n+1 caras (np.broadcast_to), para usar
    # la de las caras de contorno, índices 0 y -1
    difusividad = np.broadcast_to(difusividad, malla.n + 1)

    # 4. contorno oeste: coeficientes con malla.delta_oeste, y
    _, _, gamma_oeste, kappa_oeste = coeficientes_contorno(oeste, malla.delta_oeste)
    
    Sp[0] += difusividad[0] * malla.area[0] * gamma_oeste
    Su[0] += difusividad[0] * malla.area[0] * kappa_oeste

    # 5. contorno este: igual con malla.delta_este e índices -1
    _, _, gamma_este, kappa_este = coeficientes_contorno(este, malla.delta_este)
    
    Sp[-1] += difusividad[-1] * malla.area[-1] * gamma_este
    Su[-1] += difusividad[-1] * malla.area[-1] * kappa_este

    # 6. aP = aW + aE - Sp
    aP = aW + aE - Sp 

    # 7. resolver con tdma_coeficientes y devolver phi
    phi = tdma_coeficientes(aW=aW, aP=aP, aE=aE, Su=Su)

    return phi


if __name__ == "__main__":
    from fvm1d.malla import crear_malla

    m = crear_malla(0.5, 5, 0.01)
    phi = resolver_difusion(m, 1000.0,
                            Contorno("dirichlet", 100.0),
                            Contorno("dirichlet", 500.0))
    print(phi)

