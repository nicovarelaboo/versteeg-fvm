import numpy as np

from fvm.malla import Malla
from fvm.contornos import Contorno, coeficientes_contorno


def coeficientes_upwind(malla: Malla,
                        flujo_masico: float | np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Coeficientes convectivos con el esquema upwind en las caras interiores. 
    Los contornos se tratan en contornos_upwind.

    Parameters
    ----------
    malla : Malla
        Malla 1D centrada en las celdas.
    flujo_masico : float o np.ndarray
        rho * u * A, escalar o un valor por cara (n+1)
    
    Returns
    -------
    aW_conv : np.ndarray
        Contribución convectiva a aW en los n nodos; cero en el primero, 
        cuya cara oeste es de contorno.
    aE_conv : np.ndarray
        Contribución convectiva a aE en los n nodos; cero en el último, 
        cuya cara este es de contorno.
    aP_conv : np.ndarray
        F_e - F_w en cada celda: lo que la convección aporta a aP 
        aparte de aW + aE.
    """
    # flujo másico F en las n+1 caras (np.broadcast_to)
    F = np.broadcast_to(flujo_masico, malla.n+1)

    # F en las caras interiores (n-1 valores)
    F_interiores = F[1:-1]

    # aW_conv: n ceros; aW_conv[1:] = max(F_w, 0)
    aW_conv = np.zeros(malla.n)
    aW_conv[1:] = np.maximum(F_interiores, 0)

    # aE_conv: n ceros; aE_conv[:-1] = max(-F_e, 0)
    aE_conv = np.zeros(malla.n)
    aE_conv[:-1] = np.maximum(-F_interiores, 0)

    # aP_conv: F_e - F_w de cada celda, usando todas las caras
    aP_conv = np.diff(F)
    
    # devolver los tres arrays
    return aW_conv, aE_conv, aP_conv


def contornos_upwind(malla: Malla, flujo_masico: float | np.ndarray,
                     oeste: Contorno, este: Contorno
                     ) -> tuple[np.ndarray, np.ndarray]:
    """
    Contribución convectiva de los contornos con el esquema upwind.

    Solo contribuye el contorno por el que entra el flujo, cuyo valor
    se expresa como φ_b = α φ_P + β. Por el contorno de salida, upwind
    toma φ_P y la condición de contorno no interviene.

    Parameters
    ----------
    malla : Malla
        Malla 1D centrada en las celdas.
    flujo_masico : float o np.ndarray
        rho * u * A, escalar o un valor por cara (n+1).
    oeste, este : Contorno
        Condiciones de contorno de cada extremo.

    Returns
    -------
    Su_conv, Sp_conv : np.ndarray
        Arrays de n valores, nulos salvo en el nodo junto al contorno
        de entrada: Su = F_entrante β y Sp = F_entrante (α - 1).
    """

    # Su_conv y Sp_conv: n ceros.
    Su_conv = np.zeros(malla.n)
    Sp_conv = np.zeros(malla.n)

    # F en las n + 1 caras (np.broadcast_to).
    F = np.broadcast_to(flujo_masico, malla.n+1)

    # alpha y beta del oeste con coeficientes_contorno(oeste, malla.delta_oeste).
    alpha_oeste, beta_oeste, _, _ = coeficientes_contorno(oeste, malla.delta_oeste)

    # Flujo que entra por el oeste: max(F[0], 0).
    # Sumar en el primer nodo según las fórmulas.
    Su_conv[0] += np.maximum(F[0], 0) * beta_oeste
    Sp_conv[0] += np.maximum(F[0], 0) * (alpha_oeste - 1)

    # Lo mismo para el este con F[-1], max(-F[-1], 0) y el último nodo.
    alpha_este, beta_este, _, _ = coeficientes_contorno(este, malla.delta_este)

    Su_conv[-1] += np.maximum(-F[-1], 0) * beta_este
    Sp_conv[-1] += np.maximum(-F[-1], 0) * (alpha_este - 1)

    # Devolver Su_conv, Sp_conv.
    return Su_conv, Sp_conv


if __name__ == "__main__":
    from fvm.malla import crear_malla

    m = crear_malla(L=1.0, n=5, area=1.0)
    F_base = np.array([1.0, 1.1, 1.0, 1.2, 0.9, 1.1])
    oeste = Contorno("dirichlet", 1.0)
    este = Contorno("dirichlet", 2.0)

    for i in ['positivo', 'negativo']:
        k = +1 if i == 'positivo' else -1

        aW_conv, aE_conv, aP_conv = coeficientes_upwind(
            malla=m, flujo_masico=k * F_base)
        print(f'aW con flujo másico {i}: {aW_conv}')
        print(f'aE con flujo másico {i}: {aE_conv}')
        print(f'aP con flujo másico {i}: {aP_conv}')

        Su_conv, Sp_conv = contornos_upwind(
            malla=m, flujo_masico=k * 0.1, oeste=oeste, este=este)
        print(f'Su con flujo másico {i}: {Su_conv}')
        print(f'Sp con flujo másico {i}: {Sp_conv}')


