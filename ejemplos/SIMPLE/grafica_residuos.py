import matplotlib.pyplot as plt 
from ejemplos.SIMPLE.tobera1D import tobera1D
import os

CARPETA_FIGURAS = os.path.join(os.path.dirname(__file__), "figuras")
os.makedirs(CARPETA_FIGURAS, exist_ok=True)

combinaciones = [(0.5, 0.5), (0.5, 0.8), (0.8, 0.2), (1.0, 1.0)]

plt.figure()

for alphau, alphap in combinaciones:
    resultados = tobera1D(alphau=alphau, alphap=alphap)
    label = f'({alphau}, {alphap})'
    plt.semilogy(resultados["residuos de masa"], label = label)

plt.title(r'Residuos de masa para los pares de valores ($\alpha_u$, $\alpha_p$)')
plt.xlabel("Iteración")
plt.ylabel("Residuo de masa normalizado")
plt.grid()
plt.legend()
plt.savefig(os.path.join(CARPETA_FIGURAS, "residuos.pdf"))
plt.show()

