from dataclasses import dataclass
import numpy as np

from fvm.malla import Malla
from fvm.contribuciones import Contribucion


@dataclass(frozen=True)
class Fuente:
    """
    Fuente linealizada por unidad de volumen: S = constante + lineal·φ.
    """

    constante: float | np.ndarray
    lineal: float | np.ndarray  # por unidad de φ y de volumen (debe ser <= 0)


def contribucion_fuentes(malla: Malla,
                         fuentes: list[Fuente] | tuple[Fuente, ...]
                         ) -> Contribucion:
    """
    Integra las fuentes linealizadas en el volumen de cada celda.

    Parameters
    ----------
    malla : Malla
        Malla 1D centrada en celdas.
    fuentes : list[Fuente] | tuple[Fuente, ...]
        Fuentes a sumar; cada campo es un valor uniforme o uno por nodo.

    Returns
    -------
    Contribucion
        Su y Sp con las fuentes integradas en el volumen de cada celda;
        aW, aE y aP_extra nulos.

    Raises
    ------
    ValueError
        Si algún coeficiente lineal es positivo (regla 3 de Patankar).
    """

    aW = np.zeros(malla.n) 
    aE = np.zeros(malla.n) 
    aP_extra = np.zeros(malla.n)

    Su = np.zeros(malla.n)
    Sp = np.zeros(malla.n)

    for fuente in fuentes:
        if np.any(np.asarray(fuente.lineal) > 0):
            raise ValueError("El coeficiente lineal de una fuente debe ser <= 0 (regla 3 de Patankar)")
        
        Su += fuente.constante * malla.volumen_celdas

        Sp += fuente.lineal * malla.volumen_celdas

    return Contribucion(aW=aW, aE=aE, Su=Su, Sp=Sp, aP_extra=aP_extra)


if __name__ == "__main__":
    from fvm.malla import crear_malla

    m = crear_malla(L=0.5, n=5, area=0.01)
    fuentes = [Fuente(constante=1000.0, lineal=-10.0)]
    c = contribucion_fuentes(malla=m, fuentes=fuentes)
    print(c.Su)
    print(c.Sp)

