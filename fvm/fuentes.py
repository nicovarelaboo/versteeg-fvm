from dataclasses import dataclass
import numpy as np

from fvm.malla import Malla

@dataclass(frozen=True)
class Fuente:
    """
    Fuente linealizada por unidad de volumen: S = constante + lineal·φ.
    """

    constante: float | np.ndarray
    lineal: float | np.ndarray  # por unidad de φ y de volumen (debe ser <= 0)

def integrar_fuentes(malla: Malla, fuentes: list[Fuente]) -> tuple[np.ndarray, np.ndarray]:
    """
    Integra las fuentes linealizadas en el volumen de cada celda.

    Parameters
    ----------
    malla : Malla
        Malla 1D centrada en celdas.
    fuentes : list[Fuente]
        Fuentes a sumar; cada campo es un valor uniforme o uno por nodo.

    Returns
    -------
    Su, Sp : ndarray
        Aportaciones a Su y Sp en cada nodo.

    Raises
    ------
    ValueError
        Si algún coeficiente lineal es positivo (regla 3 de Patankar).
    """

    Su = np.zeros(malla.n)
    Sp = np.zeros(malla.n)

    for fuente in fuentes:
        if np.any(np.asarray(fuente.lineal) > 0):
            raise ValueError("El coeficiente lineal de una fuente debe ser <= 0 (regla 3 de Patankar)")
        
        Su += fuente.constante * malla.volumen_celdas

        Sp += fuente.lineal * malla.volumen_celdas

    return Su, Sp

if __name__ == "__main__":
    from fvm.malla import crear_malla

    m = crear_malla(L=0.5, n=5, area=0.01)
    fuentes = [Fuente(constante=1000.0, lineal=-10.0)]
    Su, Sp = integrar_fuentes(malla=m, fuentes=fuentes)
    print(Su)
    print(Sp)

