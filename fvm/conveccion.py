import numpy as np

from fvm.malla import Malla

def coeficientes_conveccion(malla: Malla, flujo_masico: float | np.ndarray
                            ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:

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

if __name__ == "__main__":
    from fvm.malla import crear_malla

    m = crear_malla(L=1.0, n=5, area=1.0)

    for i in ['positivo', 'negativo']:
        k = +1 if i == 'positivo' else -1
        aW_conv, aE_conv, aP_conv = coeficientes_conveccion(malla=m, flujo_masico=k*np.array([1.0, 1.1, 1.0, 1.2, 0.9, 1.1]))
        print(f'aW con flujo másico {i}: {aW_conv}')
        print(f'aE con flujo másico {i}: {aE_conv}')
        print(f'aP con flujo másico {i}: {aP_conv}')


