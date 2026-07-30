# Resultados detallados

## Entorno

- Python 3, PuLP 3.3.2, HiGHS 1.15.1 (`highspy`), OR-Tools CP-SAT.
- MILP: límite de 120 s por instancia. CP-SAT: 8 hilos de búsqueda.
- Gurobi no utilizado (no requiere licencia).

## Tabla principal

Cinco de los seis casos tienen óptimo exacto certificado por CP-SAT.

| $z$ | $R$ | Mínimo | Por ronda | Estado | Tiempo |
|:---:|:---:|:---:|:---:|---|:---:|
| 4 | 1 | **1** | 1 | `OPTIMAL` | 0,8 s |
| 4 | 2 | **4** | 2,2 | `OPTIMAL` | 2,9 s |
| 4 | 3 | **9** | 2,4,3 | `OPTIMAL` | 68 s |
| 8 | 1 | **1** | 1 | `OPTIMAL` | 2,1 s |
| 8 | 2 | **4** | 2,2 | `OPTIMAL` | 9,2 s |
| 8 | 3 | $[8, 10]$ | 2,4,4 | `FEASIBLE` | 260 s |

Para el caso sin cerrar, ambos extremos están demostrados: existe una trayectoria con 10
cajas activas, y CP-SAT probó que ninguna tiene menos de 8. Corridas más largas con
distinto número de hilos dieron incumbentes de 10 y 12 (CP-SAT no es determinista al
variar los hilos), de modo que 10 es el mejor valor conocido.

## Magnitudes de seguridad

Con $\mathrm{DP}_{\max}(\chi) = 2^{-2}$, una trayectoria con $n$ cajas activas tiene
probabilidad $\le 2^{-2n}$ y requiere $\approx 2^{2n}$ pares.

Un punto metodológico importante: la garantía de seguridad se deriva de la **cota
inferior** del número de cajas activas, no de la superior. Si el mínimo real fuera menor
que el hallado, existiría un ataque mejor. Para los cinco casos certificados la cota
inferior coincide con el mínimo exacto, de modo que la cifra es firme.

| Intentos | $R$ | $n$ | Probabilidad | Pares | Orden decimal |
|:---:|:---:|:---:|:---:|:---:|:---:|
| $< 10$ | 1 | 1 (exacto) | $2^{-2}$ | $2^{2}$ | $\sim 4$ |
| 10–19 | 2 | 4 (exacto) | $2^{-8}$ | $2^{8}$ | $\sim 2.6 \times 10^{2}$ |
| 20–29 | 3 | 9 (exacto, $z{=}4$) | $2^{-18}$ | $2^{18}$ | $\sim 2.6 \times 10^{5}$ |
| 20–29 | 3 | $\ge 8$ ($z{=}8$) | $\le 2^{-16}$ | $\ge 2^{16}$ | $\gtrsim 6.6 \times 10^{4}$ |

## Interpretación de la descomposición por ronda

El desglose es más informativo que el total:

- $R = 2$: `2, 2`.
- $R = 3$, $z=4$: `2, 4, 3` — la trayectoria óptima **vuelve a estrecharse** en la tercera
  ronda en lugar de dispersarse.
- $R = 3$, $z=8$: `2, 4, 4` (mejor trayectoria hallada).

Este patrón es relevante y contradice la intuición inicial de este trabajo. Una búsqueda
voraz que en cada caja elige la salida de menor peso de Hamming produce trayectorias que
se dispersan monótonamente (`2, 2, 14` para $z{=}4$), porque cada decisión local ignora su
efecto en las rondas siguientes. La trayectoria óptima sacrifica peso en la segunda ronda
(4 cajas en lugar de 2) para poder contraerse en la tercera (3 en lugar de 14). Es un
ejemplo claro de por qué la optimización local no basta en criptoanálisis diferencial.

## MILP frente a CP-SAT

| Caso | MILP (HiGHS, 120 s) | CP-SAT | Óptimo |
|---|:---:|:---:|:---:|
| $R{=}1$, $z{=}4$ | 1, certificado | 1, certificado | **1** |
| $R{=}1$, $z{=}8$ | 1, certificado | 1, certificado | **1** |
| $R{=}2$, $z{=}4$ | 9, gap 100 % | 4, certificado | **4** |
| $R{=}2$, $z{=}8$ | 4, gap 100 % | 4, certificado | **4** |
| $R{=}3$, $z{=}4$ | 18, gap 100 % | 9, certificado | **9** |
| $R{=}3$, $z{=}8$ | sin solución factible | $[8, 10]$ | ? |

### Por qué el MILP no cierra

La cota dual permanece en 0 durante toda la ejecución. Dos causas:

1. **Encoding *big-M* de la DDT.** El modelo emplea una variable binaria de selección por
   cada una de las 317 transiciones válidas y por cada caja-S: **12 680 de las 14 140
   variables (≈ 90 %)** para $R{=}2$, $z{=}4$, y 38 040 para $R{=}3$, $z{=}8$. Las
   restricciones *big-M* tienen relajaciones notoriamente débiles.
2. **Paridad sobre $\mathbb{F}_2$.** Al admitir valores fraccionarios, el problema
   continuo se satisface con objetivo cercano a cero mediante asignaciones que no
   corresponden a ninguna diferencia real 0/1.

CP-SAT elimina la primera causa por completo (restricción de tabla, sin variables de
selección ni *big-M*) y maneja la segunda con propagación y aprendizaje de cláusulas sobre
la estructura XOR, que trata de forma nativa.

### Estrategias que no bastaron

Antes de migrar a CP-SAT se probaron dos técnicas habituales sobre el modelo MILP:

| Estrategia | Resultado sobre $R{=}2$, $z{=}4$ |
|---|---|
| Reformulación como problema de decisión ($\sum A \le 3$) | `kTimeLimit` a 180 s |
| Decisión + ruptura de simetría en $z$ | `kTimeLimit` a 240 s |
| Refutar $\sum A \le 1$ (el caso más sencillo) | infactible, pero **92 s** |

Sobre la reformulación de decisión conviene precisar un malentendido frecuente: añadir
$\sum A \le K$ es una cota **superior** del objetivo y no aprieta la relajación por abajo
—el LP sigue satisfaciéndose con objetivo cercano a 0, que cumple trivialmente
$\le K$—. Su beneficio proviene de la propagación de la restricción de cardinalidad, no de
una relajación más ajustada.

La ruptura de simetría es válida (la ronda sin $\iota$ conmuta con la traslación en $z$),
pero la órbita tiene tamaño a lo sumo $z$, de modo que el ahorro máximo es de 4 a 8 veces:
irrelevante frente a un gap del 100 %.

Que refutar el caso más sencillo posible costara 92 s calibra la dificultad del árbol de
búsqueda, y explica por qué ninguna estrategia de búsqueda arregla un modelo cuyo cuello de
botella es el encoding.

## Efecto del reescalado de $\rho$

| $z$ | Offsets distintos (de 25 carriles) | Carriles sin rotar |
|:---:|:---:|:---:|
| 4 | 4 | 7 |
| 8 | 8 | 3 |

En Keccak-f[1600] los 25 desplazamientos son casi todos distintos, lo que evita el
alineamiento vertical de diferencias. Con $z = 4$ esa propiedad se pierde en gran medida.
Es inherente a la reducción de palabra —no un defecto de implementación— y debe tenerse en
cuenta al extrapolar los resultados a $w = 64$.

## Reproducción

```bash
python src/cpsat_keccak.py 120          # certificación (CP-SAT)
python src/experimentos.py --tiempo 120 # modelo MILP
python tests/test_verificaciones.py     # 15 verificaciones
```

Los óptimos de $R = 1$ y $R = 2$ se comprueban automáticamente en
`test_cpsat_certifica_optimos`.
