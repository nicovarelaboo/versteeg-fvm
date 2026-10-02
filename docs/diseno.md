# Diseño de fvm

Paquete de volúmenes finitos (de momento 1D) para la ecuación de transporte de φ,
siguiendo la notación de Versteeg y Malalasekera.

## Estructura del repositorio

- `fvm/`: paquete reutilizable (malla, contornos, difusión, convección, fuentes,
  ensamblado y resolvedores).
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
- Cada término trata sus propios contornos en su módulo (contornos_difusion,
  contornos_upwind) porque el tratamiento depende del esquema.
- a_P se calcula al final, después de que todos los términos hayan
  aportado sus contribuciones.
- Cada término puede aportar a a_P una parte que no se reduce a aW + aE
  (p. ej. aP_conv = F_e − F_w en convección) y la devuelve por separado;
  el ensamblado la suma junto con las demás.
- Notación: mayúsculas para nodos (aP, aW, aE), minúsculas para caras.
- Fuentes linealizadas por unidad de volumen, S = constante + lineal·φ,
  con lineal <= 0 (regla 3 de Patankar). Se integran con el volumen de
  cada celda, que guarda la malla.
- Cada esquema convectivo tiene sus propias funciones de coeficientes y
  de contornos (coeficientes_upwind, contornos_upwind). Cuando haya un
  segundo esquema, se agruparán en una dataclass o una clase Esquema.

## Diagrama de clases (pasos 1 a 3)

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
        +ndarray volumen_celdas
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
        +contornos_difusion(malla, difusividad, oeste, este) tuple
    }
    class modulo_conveccion {
        <<módulo conveccion.py>>
        +coeficientes_upwind(malla, flujo_masico) tuple
        +contornos_upwind(malla, flujo_masico, oeste, este) tuple
    }
    class modulo_ensamblado {
        <<módulo ensamblado.py>>
        +resolver_difusion(malla, difusividad, oeste, este, fuentes) ndarray
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
        +prueba_generacion_uniforme() bool
        +prueba_aleta() bool
    }
    class Fuente {
        <<dataclass frozen>>
        +float|ndarray constante
        +float|ndarray lineal
    }
    class modulo_fuentes {
        <<módulo fuentes.py>>
        +integrar_fuentes(malla, fuentes) tuple
    }

    modulo_malla ..> Malla : crea
    modulo_contornos ..> Contorno : lee
    modulo_difusion ..> Malla : lee
    modulo_conveccion ..> Malla : lee
    modulo_fuentes ..> Fuente : lee
    modulo_fuentes ..> Malla : lee
    modulo_ensamblado ..> modulo_fuentes : llama
    modulo_ensamblado ..> modulo_difusion : llama
    modulo_difusion ..> modulo_contornos : llama
    modulo_conveccion ..> modulo_contornos : llama
    modulo_ensamblado ..> modulo_resolvedores : llama
    modulo_pruebas ..> modulo_ensamblado : valida
```

## Diagrama de actividad: resolver_difusion (pasos 1 y 2)

```mermaid
flowchart TD
    A(["Inicio"]) --> B["crear_malla(L, n, area)"]
    B --> C["Contorno oeste y Contorno este"]
    C --> D["coeficientes_difusion: aW, aE en caras interiores"]
    D --> E["contornos_difusion: Su_dif, Sp_dif"]
    E --> F["integrar_fuentes: Su_fuentes, Sp_fuentes"]
    F --> G["Su = Su_dif + Su_fuentes ; Sp = Sp_dif + Sp_fuentes"]
    G --> H["aP = aW + aE − Sp"]
    H --> I["tdma_coeficientes(aW, aP, aE, Su)"]
    I --> J(["φ"])
```

## Ampliación prevista (pasos 2 a 5)

Cada término físico aporta su contribución a los mismos arrays y el
ensamblado solo las suma. Las cajas discontinuas aún no están
conectadas al ensamblado.

```mermaid
flowchart TD
    M["Malla, contornos y propiedades"] --> T1["Difusión: aW, aE, Su, Sp"]
    M --> T2["Fuentes (paso 2): Su, Sp"]
    M --> T3["Convección (paso 3): aW, aE, aP_conv, Su, Sp"]
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
    class T3,T4,T5 futuro
```

## Ideas pendientes

- Contornos de segundo orden (gradiente con tres puntos), comparando
  con el primer orden en un estudio de malla.
- Prueba de convergencia del ejemplo 4.2 frente a la solución exacta.
- Constructores alternativos con @classmethod (p. ej. Malla.uniforme).
- Empaquetar las contribuciones de los coeficientes aE, aW, Su, Sp en 
  su propia dataclass
- Agrupar las propiedades del fluido (densidad, difusividad, viscosidad)
  en una dataclass Propiedades; el flujo másico sigue aparte porque en
  SIMPLE cambia en cada iteración.
- Ensamblado con matriz dispersa (scipy.sparse), necesario en 2D.

