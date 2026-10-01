import numpy as np
import matplotlib.pyplot as plt 
from ejemplos.SIMPLE.tobera1D import tobera1D, solucion_exacta
import os

CARPETA_FIGURAS = os.path.join(os.path.dirname(__file__), "figuras")
os.makedirs(CARPETA_FIGURAS, exist_ok=True)

valores_n = [10, 20, 50]

fig_u, ax_u = plt.subplots()
fig_p, ax_p = plt.subplots()

for n in valores_n:
    resultados = tobera1D(n=n)
    label = f'n = {n}'
    ax_u.plot(resultados["xu"], resultados["u"], label = label, marker = 'o', linestyle = 'none')
    ax_p.plot(resultados["xP"], resultados["p"], label = label, marker = 'o', linestyle ='none')

x = np.linspace(0, resultados["L"], 100)
u, p = solucion_exacta(x)
ax_u.plot(x, u, label='Solución exacta')
ax_p.plot(x, p, label='Solución exacta')

ax_u.set_title('u numérica vs exacta')
ax_u.set_xlabel("x")
ax_u.set_ylabel("u")
ax_u.grid()
ax_u.legend()
fig_u.savefig(os.path.join(CARPETA_FIGURAS, "u_numerica_vs_exacta.pdf"))


ax_p.set_title('p numérica vs exacta')
ax_p.set_xlabel("x")
ax_p.set_ylabel("p")
ax_p.grid()
ax_p.legend()
fig_p.savefig(os.path.join(CARPETA_FIGURAS, "p_numerica_vs_exacta.pdf"))

plt.show()

