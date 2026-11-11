# versteeg-fvm

Implementación propia en Python del método de volúmenes finitos (FVM),
siguiendo el libro *An Introduction to Computational Fluid Dynamics: The
Finite Volume Method*, de Versteeg y Malalasekera. Es un proyecto personal
para entender el CFD por dentro, desde la discretización hasta el
acoplamiento presión-velocidad.

## Contenido

- `fvm/`: paquete reutilizable para la ecuación de transporte en 1D:
  malla centrada en celdas, condiciones de contorno de Dirichlet y
  Neumann, difusión, fuentes linealizadas (con la regla de Patankar
  del coeficiente lineal negativo) y resolvedor TDMA.
- `ejemplos/`: ejemplos del libro, por tema:
  - `difusion/`: conducción estacionaria con y sin fuentes.
  - `conveccion/`: convección-difusión con diferencias centradas,
    upwind, esquema híbrido y QUICK.
  - `SIMPLE/`: algoritmo SIMPLE en una tobera 1D.
  - `transitorio/`: conducción no estacionaria.
- `tests/`: pruebas frente a soluciones exactas.
- `docs/`: documento de diseño del paquete ([diseño](docs/diseno.md)) y
  script que genera las figuras de este README.

## Resultados

### Convección-difusión: comparación de esquemas

Con un número de Péclet de celda alto, las diferencias centradas oscilan
y QUICK se sale del rango físico, mientras que upwind se mantiene acotado
a costa de difusión numérica.

![Comparación de esquemas convectivos](docs/img/conveccion_esquemas.png)

### Aleta con convección al ambiente

Difusión con una fuente linealizada, resuelta con el paquete `fvm`:

![Aleta](docs/img/aleta.png)

### SIMPLE en una tobera 1D

Velocidad numérica frente a la solución exacta con 10, 20 y 50 celdas:

![Velocidad numérica frente a la exacta](docs/img/u_numerica_vs_exacta.png)

Residuo de masa según los factores de relajación (α_u, α_p): sin
relajación, el algoritmo diverge.

![Residuos de SIMPLE](docs/img/residuos.png)

## Uso

Requiere Python 3 con numpy y matplotlib. Los comandos se ejecutan desde
la raíz del repositorio.

### Resolver un problema con el paquete

Aleta de longitud 1 con el extremo izquierdo a 100 °C, el derecho aislado
y pérdida de calor al ambiente, modelada como una fuente linealizada
S = 500 − 25 T:

```python
from fvm.malla import crear_malla
from fvm.contornos import Contorno
from fvm.fuentes import Fuente
from fvm.ensamblado import resolver_difusion

malla = crear_malla(L=1.0, n=10, area=1.0)

T = resolver_difusion(
    malla=malla,
    difusividad=1.0,
    oeste=Contorno("dirichlet", 100.0),
    este=Contorno("neumann", 0.0),
    fuentes=[Fuente(constante=500.0, lineal=-25.0)],
)
print(T)   # temperatura en los centros de las 10 celdas
```

### Ejecutar los ejemplos, las pruebas y las figuras

    python -m ejemplos.SIMPLE.tobera1D
    python -m tests.test_difusion
    python -m docs.generar_figuras

## Próximos pasos

- Reorganizar el paquete para que cada término (difusión, convección,
  fuentes) aporte su contribución a la ecuación discreta de forma
  independiente (rama `feat/fvm-conveccion`).
- Migrar los ejemplos antiguos al paquete e incorporar a él QUICK y el término
  transitorio.
- Extender el paquete a 2D.
