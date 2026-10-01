# Diseño de fvm

Paquete de volúmenes finitos (de momento 1D) para la ecuación de transporte de φ,
siguiendo la notación de Versteeg y Malalasekera.

## Estructura del repositorio

- `fvm/`: paquete reutilizable (malla, contornos, difusión, ensamblado
  y resolvedores).
- `ejemplos/`: ejemplos del libro, organizados por tema (`difusion`,
  `conveccion`, `SIMPLE`, `transitorio`). Los más antiguos aún no usan
  el paquete y se irán migrando.
- `tests/`: pruebas contra soluciones exactas.
- `docs/`: documentación de diseño.

Todo se ejecuta como módulo desde la raíz, por ejemplo
`python -m ejemplos.difusion.ejemplo42` o `python -m tests.test_difusion`.

## Decisiones

- Las dataclasses (`frozen`) solo guardan datos; las funciones calculan.
- Todas las distancias de la malla se obtienen restando coordenadas.
- Cada contorno se reduce a cuatro coeficientes lineales en φ_P:
  φ_b = α φ_P + β y (∂φ/∂n)_b = γ φ_P + κ, con n la normal saliente.
  El resto del código no necesita saber de qué tipo es el contorno.
- a_P se calcula al final, después de que todos los términos hayan
  aportado sus contribuciones.
- Notación: mayúsculas para nodos (aP, aW, aE), minúsculas para caras.

## Diagrama de clases (paso 1)

```mermaid
classDiagram
    class Malla {
        <<dataclass frozen>>
        +float L
        +int n
        +ndarray caras
        +ndarray nodos
        +ndarray dx
        +ndarray delta_nodos
        +float delta_oeste
        +float delta_este
        +ndarray area
    }
    class Contorno {
        <<dataclass frozen>>
        +str tipo
        +float valor
    }
    class modulo_malla {
        <<módulo malla.py>>
        +crear_malla(L, n, area) Malla
    }
    class modulo_contornos {
        <<módulo contornos.py>>
        +coeficientes_contorno(contorno, delta) tuple
    }
    class modulo_difusion {
        <<módulo difusion.py>>
        +coeficientes_difusion(malla, difusividad) tuple
    }
    class modulo_ensamblado {
        <<módulo ensamblado.py>>
        +resolver_difusion(malla, difusividad, oeste, este) ndarray
    }
    class modulo_resolvedores {
        <<módulo resolvedores.py>>
        +TDMA(A, b) ndarray
        +tdma_coeficientes(aW, aP, aE, Su) ndarray
    }
    class modulo_pruebas {
        <<módulo tests/test_difusion.py>>
        +prueba_flujo_este() bool
        +prueba_flujo_oeste() bool
        +prueba_aislamiento_este() bool
    }

    modulo_malla ..> Malla : crea
    modulo_contornos ..> Contorno : lee
    modulo_difusion ..> Malla : lee
    modulo_ensamblado ..> modulo_difusion : llama
    modulo_ensamblado ..> modulo_contornos : llama
    modulo_ensamblado ..> modulo_resolvedores : llama
    modulo_pruebas ..> modulo_ensamblado : valida
```

## Diagrama de actividad: resolver_difusion (paso 1)

```mermaid
flowchart TD
    A(["Inicio"]) --> B["crear_malla(L, n, area)"]
    B --> C["Contorno oeste y Contorno este"]
    C --> D["coeficientes_difusion: aW, aE en caras interiores"]
    D --> E["Su = 0, Sp = 0"]
    E --> F["coeficientes_contorno(oeste, delta_oeste): γ, κ"]
    F --> G["Sp[0] += Γ A γ ; Su[0] += Γ A κ"]
    G --> H["coeficientes_contorno(este, delta_este): γ, κ"]
    H --> I["Sp[-1] += Γ A γ ; Su[-1] += Γ A κ"]
    I --> J["aP = aW + aE − Sp"]
    J --> K["tdma_coeficientes(aW, aP, aE, Su)"]
    K --> L(["φ"])
```

## Ampliación prevista (pasos 2 a 5)

Cada término físico aporta su contribución a los mismos arrays y el
ensamblado solo las suma. Las cajas discontinuas aún no existen.

```mermaid
flowchart TD
    M["Malla, contornos y propiedades"] --> T1["Difusión: aW, aE, Su, Sp"]
    M --> T2["Fuentes (paso 2): Su, Sp"]
    M --> T3["Convección (paso 3): aW, aE, Su, Sp"]
    M --> T4["Corrección diferida QUICK (paso 4): Su"]
    M --> T5["Transitorio (paso 5): aP0, Su"]
    T1 --> S["Sumar contribuciones"]
    T2 --> S
    T3 --> S
    T4 --> S
    T5 --> S
    S --> P["aP = aW + aE + ΔF + aP0 − Sp"]
    P --> R["Resolver (TDMA o matriz dispersa)"]
    R --> X(["φ"])

    classDef futuro stroke-dasharray: 5 5
    class T2,T3,T4,T5 futuro
```