import numpy as np
from fvm.malla import Malla
from fvm.contornos import Contorno, coeficientes_contorno
from fvm.contribuciones import Contribucion


def contribucion_difusion(malla: Malla, difusividad: float | np.ndarray,
                          oeste: Contorno, este: Contorno) -> Contribucion:
    """Contribución de la difusión a la ecuación discreta de cada celda.

    Calcula los coeficientes de las caras interiores y el efecto de los
    contornos sobre el primer y el último nodo.

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
    Contribucion
        aW y aE con las conductancias de las caras interiores (aW[0] y
        aE[-1] nulos); Su y Sp nulos salvo en el primer y el último
        nodo; aP_extra nulo.
    """

    # Difusividad en las n + 1 caras
    difusividad = np.broadcast_to(difusividad, malla.n + 1)

    # ---------- Caras interiores: conductancias, aW y aE ----------
    conductancia_difusiva = difusividad[1:-1] * malla.area[1:-1] / malla.delta_nodos

    # aW: n ceros, y rellenar aW[1:] con las conductancias
    aW_dif = np.zeros(malla.n)
    aW_dif[1:] = conductancia_difusiva

    # aE: n ceros, y rellenar aE[:-1] con las conductancias 
    aE_dif = np.zeros(malla.n)
    aE_dif[:-1] = conductancia_difusiva

    # ---------- Contornos ----------
    # Sp_dif Su_dif: arrays de n ceros
    Sp_dif = np.zeros(malla.n)
    Su_dif = np.zeros(malla.n)

    # Pedir gamma y kappa a coeficientes_contorno para el oeste
    # y sumarlos en el primer nodo.
    _, _, gamma_oeste, kappa_oeste = coeficientes_contorno(oeste, malla.delta_oeste)

    Sp_dif[0] += difusividad[0] * malla.area[0] * gamma_oeste
    Su_dif[0] += difusividad[0] * malla.area[0] * kappa_oeste

    # Lo mismo para el este, en el último nodo.
    _, _, gamma_este, kappa_este = coeficientes_contorno(este, malla.delta_este)
    
    Sp_dif[-1] += difusividad[-1] * malla.area[-1] * gamma_este
    Su_dif[-1] += difusividad[-1] * malla.area[-1] * kappa_este

    # aP_extra: ceros, la difusión no aporta nada a aP aparte de aW + aE
    aP_extra_dif = np.zeros(malla.n)

    # Devolver la Contribucion pasando los campos por nombre (aW=aW, ...)
    return Contribucion(aW=aW_dif, aE=aE_dif, Su=Su_dif, Sp=Sp_dif, aP_extra=aP_extra_dif)

if __name__ == "__main__":
    from fvm.malla import crear_malla

    m = crear_malla(0.5, 5, 0.01)
    oeste = Contorno(tipo='dirichlet', valor=100.0)
    este = Contorno(tipo='dirichlet', valor=500.0)
    c = contribucion_difusion(m, 1000.0, oeste, este)
    print(c.aW)
    print(c.aE)
    print(c.Su)
    print(c.Sp)

