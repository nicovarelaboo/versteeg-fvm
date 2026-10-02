import numpy as np
from fvm.malla import Malla
from fvm.contornos import Contorno, coeficientes_contorno

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

def contornos_difusion(malla: Malla, difusividad: float | np.ndarray,
                       oeste: Contorno, este: Contorno
                       ) -> tuple[np.ndarray, np.ndarray]:
    """
    Contribución de la difusión a los contornos.

    Parameters
    ----------
    malla : Malla
        Malla 1D centrada en celdas.
    difusividad : float o np.ndarray
        Gamma, escalar o un valor por cara (n + 1).
    oeste, este : Contorno
        Condiciones de contorno de cada extremo.

    Returns
    -------
    Su_dif, Sp_dif : np.ndarray
        Arrays de n valores, nulos salvo en el primer y el último nodo.
    """
    # 1. Arrays de ceros de n valores para Su_dif y Sp_dif.
    Su_dif = np.zeros(malla.n)
    Sp_dif = np.zeros(malla.n)

    # 2. Difusividad en las n + 1 caras (igual que en el ensamblado).
    difusividad = np.broadcast_to(difusividad, malla.n + 1)

    # 3. Pedir gamma y kappa a coeficientes_contorno para el oeste
    #    y sumarlos en el primer nodo.
    _, _, gamma_oeste, kappa_oeste = coeficientes_contorno(oeste, malla.delta_oeste)

    Sp_dif[0] += difusividad[0] * malla.area[0] * gamma_oeste
    Su_dif[0] += difusividad[0] * malla.area[0] * kappa_oeste

    # 4. Lo mismo para el este, en el último nodo.
    _, _, gamma_este, kappa_este = coeficientes_contorno(este, malla.delta_este)
    
    Sp_dif[-1] += difusividad[-1] * malla.area[-1] * gamma_este
    Su_dif[-1] += difusividad[-1] * malla.area[-1] * kappa_este

    # 5. Devolver Su_dif, Sp_dif.

    return Su_dif, Sp_dif

if __name__ == "__main__":
    from fvm.malla import crear_malla

    m = crear_malla(0.5, 5, 0.01)
    aW, aE = coeficientes_difusion(m, 1000.0)
    print(aW)
    print(aE)

