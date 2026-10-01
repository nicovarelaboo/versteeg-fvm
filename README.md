# Versteeg-Malalasekera FVM

Implementaciones en Python de los distintos ejemplos del método de
volúmenes finitos (FVM), siguiendo el libro "An Introduction to
Computational Fluid Dynamics" de Versteeg y Malalasekera.

## Estructura

- `fvm/`: paquete reutilizable de volúmenes finitos.
- `ejemplos/`: ejemplos del libro, por tema.
- `tests/`: pruebas contra soluciones exactas.
- `docs/`: documentación de diseño ([diseño](docs/diseno.md)).

## Uso

Desde la raíz del repositorio, con las dependencias instaladas
(numpy y matplotlib):

    python -m ejemplos.difusion.ejemplo42
    python -m tests.test_difusion

