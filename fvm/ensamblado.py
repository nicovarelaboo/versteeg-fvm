import numpy as np

from fvm.contornos import Contorno
from fvm.difusion import coeficientes_difusion, contornos_difusion
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
    
    aW, aE = coeficientes_difusion(malla, difusividad)

    Su_dif, Sp_dif = contornos_difusion(malla, difusividad, oeste, este)

    Su_fuentes, Sp_fuentes = integrar_fuentes(malla=malla, fuentes=fuentes)

    Sp = Sp_dif + Sp_fuentes
    Su = Su_dif + Su_fuentes

    aP = aW + aE - Sp 

    phi = tdma_coeficientes(aW=aW, aP=aP, aE=aE, Su=Su)

    return phi


if __name__ == "__main__":
    from fvm.malla import crear_malla

    m = crear_malla(L=0.5, n=5, area=0.01)
    phi = resolver_difusion(malla=m, difusividad=1000.0,
                            oeste=Contorno("dirichlet", 100.0),
                            este=Contorno("dirichlet", 500.0))
    print(phi)

