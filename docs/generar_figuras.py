"""
Genera las figuras del README en docs/img/.

Uso, desde la raíz del repositorio:

    python -m docs.generar_figuras
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from fvm.malla import crear_malla
from fvm.contornos import Contorno
from fvm.fuentes import Fuente
from fvm.ensamblado import resolver_difusion
from ejemplos.conveccion.coeficientes_cds import coeficientes_cds
from ejemplos.conveccion.coeficientes_quick import coeficientes_quick
from ejemplos.SIMPLE.tobera1D import tobera1D, solucion_exacta

CARPETA = Path(__file__).parent / "img"


def coeficientes_upwind(n, u, L, rho, Gamma, phiA, phiB):
    """
    Matriz A y vector b de convección-difusión 1D con upwind (u > 0),
    con φ = phiA en x = 0 y φ = phiB en x = L (ejemplo 5.2 del libro).
    Misma interfaz que coeficientes_cds y coeficientes_quick.
    """
    deltax = L / n
    F = rho * u
    D = Gamma / deltax

    aW = np.full(n, D + F)
    aE = np.full(n, D)
    aW[0] = 0
    aE[n-1] = 0

    Sp = np.zeros(n)
    Su = np.zeros(n)
    Sp[0] = -(2 * D + F)
    Su[0] = (2 * D + F) * phiA
    Sp[n-1] = -2 * D
    Su[n-1] = 2 * D * phiB

    aP = aW + aE - Sp   # F_e - F_w = 0 con u constante

    A = np.diag(aP) - np.diag(aW[1:], -1) - np.diag(aE[:-1], 1)
    return A, Su


def figura_conveccion():
    """Diferencias centradas, upwind y QUICK frente a la solución exacta."""
    L, n, rho, Gamma, u = 1.0, 5, 1.0, 0.1, 2.5
    phiA, phiB = 1.0, 0.0

    Pe = rho * u * L / Gamma           # Péclet del dominio
    Pe_celda = rho * u * (L / n) / Gamma

    x = np.linspace(0, L, 200)
    exacta = phiA + (phiB - phiA) * (np.exp(Pe * x / L) - 1) / (np.exp(Pe) - 1)
    nodos = (np.arange(n) + 0.5) * L / n

    plt.figure()
    plt.plot(x, exacta, "k-", label="Solución exacta")
    for nombre, coeficientes, marca in [("Diferencias centradas", coeficientes_cds, "o"),
                                        ("Upwind", coeficientes_upwind, "s"),
                                        ("QUICK", coeficientes_quick, "^")]:
        A, b = coeficientes(n, u, L, rho, Gamma, phiA, phiB)
        phi = np.linalg.solve(A, b)
        plt.plot(nodos, phi, marca + "--", label=nombre)

    plt.xlabel("x")
    plt.ylabel("φ")
    plt.title(f"Convección-difusión 1D con {n} celdas (Péclet de celda = {Pe_celda:g})")
    plt.legend()
    plt.grid()
    plt.savefig(CARPETA / "conveccion_esquemas.png", dpi=120, bbox_inches="tight")
    plt.close()


def figura_aleta():
    """Aleta con convección al ambiente (ejemplo 4.3), resuelta con el paquete fvm."""
    L, n2, TA, Tinf = 1.0, 25.0, 100.0, 20.0

    x = np.linspace(0, L, 200)
    exacta = (TA - Tinf) * np.cosh(np.sqrt(n2) * (L - x)) / np.cosh(np.sqrt(n2) * L) + Tinf

    plt.figure()
    plt.plot(x, exacta, "k-", label="Solución exacta")
    for n, marca in [(5, "s"), (10, "o")]:
        m = crear_malla(L=L, n=n, area=1.0)
        T = resolver_difusion(malla=m, difusividad=1.0,
                              oeste=Contorno("dirichlet", TA),
                              este=Contorno("neumann", 0.0),
                              fuentes=(Fuente(constante=n2 * Tinf, lineal=-n2),))
        plt.plot(m.nodos, T, marca, label=f"Volúmenes finitos, n = {n}")

    plt.xlabel("x (m)")
    plt.ylabel("T (°C)")
    plt.title("Aleta con convección al ambiente")
    plt.legend()
    plt.grid()
    plt.savefig(CARPETA / "aleta.png", dpi=120, bbox_inches="tight")
    plt.close()


def figura_tobera_velocidad():
    """Velocidad de SIMPLE en la tobera 1D frente a la solución exacta."""
    plt.figure()
    for n, marca in [(10, "s"), (20, "^"), (50, "o")]:
        resultados = tobera1D(n=n)
        plt.plot(resultados["xu"], resultados["u"], marca, label=f"SIMPLE, n = {n}")

    x = np.linspace(0, resultados["L"], 200)
    u, _ = solucion_exacta(x)
    plt.plot(x, u, "k-", label="Solución exacta")

    plt.xlabel("x")
    plt.ylabel("u")
    plt.title("Tobera 1D: velocidad con SIMPLE")
    plt.legend()
    plt.grid()
    plt.savefig(CARPETA / "u_numerica_vs_exacta.png", dpi=120, bbox_inches="tight")
    plt.close()


def figura_tobera_residuos():
    """Residuo de masa de SIMPLE para varios factores de relajación."""
    plt.figure()
    for alphau, alphap in [(0.5, 0.5), (0.5, 0.8), (0.8, 0.2), (1.0, 1.0)]:
        # Sin relajación (1.0, 1.0) el algoritmo diverge: se ignoran los
        # avisos de desbordamiento que produce.
        with np.errstate(all="ignore"):
            resultados = tobera1D(alphau=alphau, alphap=alphap)
        plt.semilogy(resultados["residuos de masa"],
                     label=rf"$\alpha_u$ = {alphau}, $\alpha_p$ = {alphap}")

    plt.xlabel("Iteración")
    plt.ylabel("Residuo de masa normalizado")
    plt.title("Tobera 1D: convergencia de SIMPLE")
    plt.legend()
    plt.grid()
    plt.savefig(CARPETA / "residuos.png", dpi=120, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    CARPETA.mkdir(exist_ok=True)
    figura_conveccion()
    figura_aleta()
    figura_tobera_velocidad()
    figura_tobera_residuos()
    print(f"Figuras guardadas en {CARPETA}")
