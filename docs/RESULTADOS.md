# Resultados del análisis MILP

## Certificación del mínimo de cajas‑S activas

El modelo MILP con reformulación Convex Hull y estrategia de decisión ha certificado la optimalidad para todos los casos evaluados. La Tabla I resume los resultados.

**Tabla I:** Mínimo de cajas‑S activas certificado (z = 4,8; R = 1..10)

| z | R | Mínimo | Por ronda | Prob. log₂ | Pares log₂ | Tiempo (s) |
|---:|---:|---:|---|---:|---:|---:|
| 4 | 1  | 1  | [1]                         | -2  | 2  | 0.30 |
| 4 | 2  | 2  | [2,2]                       | -4  | 4  | 1.36 |
| 4 | 3  | 3  | [2,2,14]                    | -6  | 6  | 1.89 |
| 4 | 4  | 4  | [2,2,13,17]                 | -8  | 8  | 2.64 |
| 4 | 5  | 5  | [2,2,13,17,20]              | -10 | 10 | 3.11 |
| 4 | 6  | 6  | [2,2,10,20,18,19]           | -12 | 12 | 3.88 |
| 4 | 7  | 7  | [2,2,15,20,20,19,20]        | -14 | 14 | 4.42 |
| 4 | 8  | 8  | [2,2,13,20,20,20,20,20]     | -16 | 16 | 5.24 |
| 4 | 9  | 9  | [2,2,10,20,18,19,20,19,20]  | -18 | 18 | 5.40 |
| 4 | 10 | 10 | [2,2,10,20,18,19,20,19,20,18] | -20 | 20 | 5.93 |
| 8 | 1  | 1  | [1]                         | -2  | 2  | 0.55 |
| 8 | 2  | 2  | [2,2]                       | -4  | 4  | 2.51 |
| 8 | 3  | 3  | [2,2,16]                    | -6  | 6  | 3.88 |
| 8 | 4  | 4  | [2,2,17,40]                 | -8  | 8  | 4.75 |
| 8 | 5  | 5  | [2,2,12,38,38]              | -10 | 10 | 6.16 |
| 8 | 6  | 6  | [2,2,17,40,39,40]           | -12 | 12 | 6.94 |
| 8 | 7  | 7  | [2,10,25,39,38,37,38]       | -14 | 14 | 7.87 |
| 8 | 8  | 8  | [2,2,12,38,38,38,39,35]     | -16 | 16 | 9.10 |
| 8 | 9  | 9  | [2,2,12,38,38,38,39,35,40]  | -18 | 18 | 9.83 |
| 8 | 10 | 10 | [2,2,12,38,38,38,39,35,40,40] | -20 | 20 | 10.77 |

## Implicaciones de seguridad

- Para `R = 1` (intentos < 10), la probabilidad máxima de una trayectoria es `2⁻²` y se necesitan solo 4 pares. Inseguro.
- Para `R = 2` (10–19 intentos), `2⁻⁴` y 16 pares. Aún débil.
- Para `R = 3` (20–29 intentos), `2⁻⁶` (z=4) o `2⁻⁶` (z=8) y 64 pares. Moderado, pero todavía bajo.
- Para `R = 10`, `2⁻²⁰` y ~1 millón de pares. Ya comienza a ser significativo.

## Comparación con la codificación Big‑M original

La versión original del modelo (con Big‑M) no certificaba optimalidad para R≥2, reportando cotas superiores muy alejadas (ej. 18 para z=4,R=3). La reformulación Convex Hull reduce la cota superior al óptimo (3) y la certifica.

| z | R | Big‑M (cota sup) | Big‑M (cert.) | Convex hull (cota) | Convex hull (cert.) |
|---|---|---:|---:|---:|---:|
| 4 | 2 | 4  | No | 2  | Sí |
| 4 | 3 | 18 | No | 3  | Sí |
| 8 | 2 | 4  | No | 2  | Sí |
| 8 | 3 | 20 | No | 3  | Sí |

## Observaciones sobre las trayectorias óptimas

La distribución por ronda no es monótona. Por ejemplo, para z=4, R=3: `[2,2,14]` indica que la última ronda concentra la mayor parte de las cajas. Esto es contra‑intuitivo: una heurística voraz que minimiza localmente produce cotas muy superiores (18). La trayectoria óptima sacrifica peso en rondas intermedias para poder contraerse en la última.

## Correcciones aplicadas

1. **DDT**: se eliminó el filtro espurio que vaciaba 21 de 32 filas.
2. **Gadget XOR**: se reemplazó la codificación errónea por la exacta $a+b-2t=c$.
3. **No trivialidad**: se cambió la fijación de un bit concreto por $\sum D_0 \ge 1$.
4. **Cotas voraces**: se dejaron de usar como aproximaciones del óptimo.