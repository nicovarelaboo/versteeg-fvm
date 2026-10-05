"""
    Define la contribución de cada término de la ecuación general
    de transporte (difusión, convección, ...)
"""

from dataclasses import dataclass
import numpy as np 


@dataclass(frozen=True)
class Contribucion:
    """Contribución de un término físico a la ecuación discreta de una celda.

    Cada término de la ecuación de transporte (difusión, convección,
    fuentes, transitorio) devuelve una ``Contribucion``. El ensamblado
    las suma campo a campo y calcula al final

    aP = aW + aE + aP_extra - Sp

    Todos los arrays tienen una componente por celda.

    Attributes
    ----------
    aW : np.ndarray
        Coeficiente del nodo vecino oeste.
    aE : np.ndarray
        Coeficiente del nodo vecino este.
    Su : np.ndarray
        Parte constante de la fuente, con unidades de coeficiente por φ.
        Incluye las condiciones de contorno; las fuentes por unidad de
        volumen van multiplicadas por el volumen de la celda.
    Sp : np.ndarray
        Parte lineal de la fuente, con unidades de coeficiente.
        Incluye las condiciones de contorno; las fuentes por unidad de
        volumen van multiplicadas por el volumen de la celda. Debe ser
        <= 0 para que aP sea positivo y el sistema, estable.
    aP_extra : np.ndarray
        Parte de aP que no proviene de aW + aE. En convección es
        F_e - F_w; en difusión, ceros.
    """

    aW: np.ndarray
    aE: np.ndarray 
    Su: np.ndarray 
    Sp: np.ndarray 
    aP_extra: np.ndarray 





