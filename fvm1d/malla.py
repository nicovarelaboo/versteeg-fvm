from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class Malla:
    L: float
    n: int
    caras: np.ndarray
    nodos: np.ndarray
    dx: np.ndarray
    delta_nodos: np.ndarray
    delta_oeste: float
    delta_este: float
    area: np.ndarray

def crear_malla(L: float, n: int, area: float | np.ndarray) -> Malla:
    caras = np.linspace(0, L, n+1)
    nodos = (caras[:-1] + caras[1:]) / 2
    dx = np.diff(caras)
    delta_nodos = np.diff(nodos)
    delta_este = caras[-1] - nodos[-1]
    delta_oeste = nodos[0] - caras[0]
    area = np.broadcast_to(area, n + 1)
    return Malla(L=L, n=n, caras=caras, nodos=nodos, dx=dx, delta_nodos=delta_nodos,
                 delta_oeste=delta_oeste, delta_este=delta_este, area=area)

if __name__ == "__main__":
    m = crear_malla(0.5, 5, 1.0)
    print(m)

