import numpy as np

from fvm.contornos import coeficientes_contorno, Contorno
from fvm.difusion import coeficientes_difusion
from fvm.resolvedores import tdma_coeficientes
from fvm.malla import Malla 
from fvm.fuentes import Fuente, integrar_fuentes


def resolver_difusion(malla: Malla, difusividad: float | np.ndarray,
                      oeste: Contorno, este: Contorno, 
                      fuentes: tuple[Fuente, ...] | list[Fuente] = ()) -> np.ndarray:
    
    """
    Resuelve la difusión estacionaria 1D de φ con volúmenes finitos.

    Ecuación: d/dx(Γ dφ/dx) + S = 0, discretizada como
    aP φP = aW φW + aE φE + Su.

    Parameters
    ----------
    malla : Malla
        Malla 1D centrada en celdas.
    difusividad : float o ndarray
        Γ, un valor constante para todo el dominio o uno por cara (n+1 valores).
    oeste, este : Contorno
        Condiciones de contorno en x = 0 y x = L.
    fuentes : tuple o list de Fuente (opcional)
        Fuentes linealizadas; por defecto, ninguna.

    Returns
    -------
    ndarray
        φ en los n nodos.
    """
    
    # 1. aW y aE de las caras interiores (coeficientes_difusion)
    aW, aE = coeficientes_difusion(malla, difusividad)

    # 1. Su y Sp: arrays de n ceros
    Su_contorno = np.zeros(malla.n)
    Sp_contorno = np.zeros(malla.n)

    # 3. difusividad en las n+1 caras (np.broadcast_to), para usar
    # la de las caras de contorno, índices 0 y -1
    difusividad = np.broadcast_to(difusividad, malla.n + 1)

    # 4. contorno oeste: coeficientes con malla.delta_oeste, y
    _, _, gamma_oeste, kappa_oeste = coeficientes_contorno(oeste, malla.delta_oeste)
    
    Sp_contorno[0] += difusividad[0] * malla.area[0] * gamma_oeste
    Su_contorno[0] += difusividad[0] * malla.area[0] * kappa_oeste

    # 5. contorno este: igual con malla.delta_este e índices -1
    _, _, gamma_este, kappa_este = coeficientes_contorno(este, malla.delta_este)
    
    Sp_contorno[-1] += difusividad[-1] * malla.area[-1] * gamma_este
    Su_contorno[-1] += difusividad[-1] * malla.area[-1] * kappa_este

    # 6. sumar fuentes
    Su_fuentes, Sp_fuentes = integrar_fuentes(malla=malla, fuentes=fuentes)
    Sp = Sp_contorno + Sp_fuentes
    Su = Su_contorno + Su_fuentes

    # 7. aP = aW + aE - Sp
    aP = aW + aE - Sp 

    # 8. resolver con tdma_coeficientes y devolver phi
    phi = tdma_coeficientes(aW=aW, aP=aP, aE=aE, Su=Su)

    return phi


if __name__ == "__main__":
    from fvm.malla import crear_malla

    m = crear_malla(L=0.5, n=5, area=0.01)
    phi = resolver_difusion(malla=m, difusividad=1000.0,
                            oeste=Contorno("dirichlet", 100.0),
                            este=Contorno("dirichlet", 500.0))
    print(phi)

