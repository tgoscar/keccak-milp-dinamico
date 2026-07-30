# Resultados detallados

## Entorno

- Python 3, OR-Tools CP-SAT (motor de certificación), PuLP y HiGHS (modelo MILP
  de comparación).
- Los casos $R \le 3$ se certificaron en la máquina de referencia del proyecto
  ($R{=}3$, $z{=}4$ en 7,2 s; $z{=}8$ en 117,6 s).
- Los casos $R \ge 4$ se midieron con un presupuesto de unos 35 s por caso; sus
  cotas mejoran con más tiempo.

## Resultados por número de rondas

### $z = 4$ (Keccak-f[100], 20 cajas-S por ronda)

| $R$ | Mínimo | Por ronda | Estado | Pares garantizados |
|:---:|:---:|:---:|---|:---:|
| 1 | **1** (exacto) | 1 | certificado | $2^{2}$ |
| 2 | **4** (exacto) | 2,2 | certificado | $2^{8}$ |
| 3 | **9** (exacto) | 2,4,3 | certificado | $2^{18}$ |
| 4 | $[9, 14]$ | - | cota | $2^{18}$ |
| 5 | $[9, 35]$ | - | cota | $2^{18}$ |
| 6 | $[9, 61]$ | - | cota | $2^{18}$ |
| 7 | $[11, 64]$ | - | cota | $2^{22}$ |
| 8 | $[13, 86]$ | - | cota | $2^{26}$ |
| 9 | $[13, 103]$ | - | cota | $2^{26}$ |
| 10 | $[13, 111]$ | - | cota | $2^{26}$ |

### $z = 8$ (Keccak-f[200], 40 cajas-S por ronda)

| $R$ | Mínimo | Por ronda | Estado | Pares garantizados |
|:---:|:---:|:---:|---|:---:|
| 1 | **1** (exacto) | 1 | certificado | $2^{2}$ |
| 2 | **4** (exacto) | 2,2 | certificado | $2^{8}$ |
| 3 | **10** (exacto) | 2,4,4 | certificado | $2^{20}$ |
| 4 | $[10, 60]$ | - | cota | $2^{20}$ |
| 5 | $[10, 62]$ | - | cota | $2^{20}$ |
| 6 | $[10, 106]$ | - | cota | $2^{20}$ |
| 7 | $[10, 156]$ | - | cota | $2^{20}$ |
| 8 | $[10, 238]$ | - | cota | $2^{20}$ |
| 9 | $[10, 239]$ | - | cota | $2^{20}$ |
| 10 | $[10, 282]$ | - | cota | $2^{20}$ |

## Cómo leer estas tablas

**Los seis casos con $R \le 3$ tienen óptimo exacto certificado** (CP-SAT devuelve
`OPTIMAL`, gap cerrado). Para $R \ge 4$ se reportan intervalos
[cota inferior, cota superior]: el extremo inferior está demostrado y el superior
corresponde a la mejor trayectoria hallada dentro del presupuesto.

**La columna de pares garantizados se deriva de la cota INFERIOR**, no de la
superior. Es un punto metodológico importante: si el mínimo real fuera menor que
el hallado, existiría un ataque mejor. Sólo la cota inferior sustenta una
garantía frente a un atacante. Para los casos certificados ambas coinciden.

**Las cotas inferiores se propagan por monotonía.** El mínimo $m(R)$ es no
decreciente en $R$: truncar una trayectoria de $R+1$ rondas produce una de $R$
rondas con un número de cajas activas menor o igual, luego $m(R+1) \ge m(R)$. En
consecuencia

$$m(R) \ge \max\{ LB(r) : r \le R \}$$

Por eso el óptimo certificado $m(3) = 9$ ($z{=}4$) garantiza $m(R) \ge 9$ para todo
$R \ge 3$, y las cotas obtenidas en $R{=}7$ y $R{=}8$ (11 y 13) elevan la garantía
para todas las rondas superiores. La consolidación está implementada en
`consolidar()` y verificada en `test_consolidacion_monotona_de_cotas`.

**Las cotas superiores de $R \ge 4$ son flojas** (hasta 282 para $R{=}10$, $z{=}8$).
No deben interpretarse como estimaciones del óptimo: son simplemente trayectorias
válidas halladas rápido. La experiencia de este proyecto en $R{=}3$ es ilustrativa:
una cota superior de 18 convivía con un óptimo real de 9.

## Interpretación de la descomposición por ronda

- $R = 2$: `2, 2`.
- $R = 3$, $z{=}4$: `2, 4, 3` — la trayectoria óptima **vuelve a estrecharse** en la
  última ronda en lugar de dispersarse.
- $R = 3$, $z{=}8$: `2, 4, 4`.

Este patrón contradice la intuición voraz. Una búsqueda que en cada caja elige la
salida de menor peso de Hamming produce trayectorias de dispersión monótona
(`2, 2, 14`, con 18 cajas en total), porque cada decisión local ignora su efecto
sobre las rondas siguientes. La trayectoria óptima sacrifica peso en la segunda
ronda (4 cajas en lugar de 2) para poder contraerse en la tercera (3 en lugar de
14). Es un ejemplo claro de por qué la optimización local no basta en
criptoanálisis diferencial.

## MILP frente a CP-SAT

| Caso | MILP big-M (HiGHS, 120 s) | MILP envolvente convexa | CP-SAT | Óptimo |
|---|:---:|:---:|:---:|:---:|
| $R{=}1$, $z{=}4$ | 1, certificado | 1 | 1, certificado | **1** |
| $R{=}2$, $z{=}4$ | 9, gap 100 % | 7, sin certificar (200 s) | 4, certificado | **4** |
| $R{=}3$, $z{=}4$ | 18, gap 100 % | no evaluado | 9, certificado | **9** |
| $R{=}3$, $z{=}8$ | sin solución | no evaluado | 10, certificado | **10** |

### Por qué el MILP no cierra

La cota dual permanece en 0 durante toda la ejecución. Dos causas:

1. **Encoding *big-M* de la DDT.** Una variable binaria de selección por cada una
   de las 317 transiciones y por cada caja-S: **12 680 de las 14 140 variables
   (≈ 90 %)** para $R{=}2$, $z{=}4$. Las restricciones *big-M* relajan mal.
2. **Paridad sobre $\mathbb{F}_2$.** Al admitir valores fraccionarios, el problema
   continuo se satisface con objetivo cercano a cero mediante asignaciones que no
   corresponden a ninguna diferencia real 0/1.

La reformulación mediante **envolvente convexa** elimina la primera causa y es
válida si se impone componente a componente (véase `docs/MODELO.md`, §10). Medida:
para $R{=}2$, $z{=}4$ devuelve 7 sin certificar en 200 s. Es decir, mejora sobre el
*big-M* (7 frente a 9) pero sigue lejos de CP-SAT, que certifica 4 en segundos.

### Estrategias que no bastaron

| Estrategia sobre el MILP | Resultado en $R{=}2$, $z{=}4$ |
|---|---|
| Reformulación como problema de decisión ($\sum A \le 3$) | `kTimeLimit` a 180 s |
| Decisión + ruptura de simetría en $z$ | `kTimeLimit` a 240 s |
| Refutar $\sum A \le 1$ (el caso más sencillo) | infactible, pero **92 s** |

Sobre la reformulación de decisión conviene precisar un malentendido frecuente:
añadir $\sum A \le K$ es una cota **superior** del objetivo y no aprieta la
relajación por abajo, pues el LP sigue satisfaciéndose con objetivo cercano a 0.
Su beneficio proviene de la propagación de la cardinalidad. La ruptura de simetría
es válida, pero la órbita tiene tamaño a lo sumo $z$, luego el ahorro máximo es de
4 a 8 veces.

## Paralelismo

`src/cpsat_keccak.py` admite `--workers` (procesos, casos simultáneos) y
`--hilos` (hilos internos de CP-SAT por caso). Sin `--hilos`, los núcleos se
reparten entre los procesos.

```bash
python src/cpsat_keccak.py --rondas 1-10 --tiempo 300 --workers 4
```

Conviene tener presente que CP-SAT **no es determinista** al variar el número de
hilos: en pruebas de $R{=}3$, $z{=}8$ se obtuvo un incumbente de 10 con 8 hilos y
de 12 con 16. Para casos difíciles vale la pena reintentar con configuraciones
distintas antes de concluir que no cierran.

## Efecto del reescalado de $\rho$

| $z$ | Offsets distintos (de 25 carriles) | Carriles sin rotar |
|:---:|:---:|:---:|
| 4 | 4 | 7 |
| 8 | 8 | 3 |

En Keccak-f[1600] los 25 desplazamientos son casi todos distintos, lo que evita el
alineamiento vertical de diferencias. Con $z = 4$ esa propiedad se pierde en gran
medida. Es inherente a la reducción de palabra y debe tenerse en cuenta al
extrapolar a $w = 64$.

## Reproducción

```bash
python src/cpsat_keccak.py --rondas 1-3 --tiempo 300        # certificación
python src/cpsat_keccak.py --rondas 1-10 --workers 4        # exploración
python src/experimentos.py --tiempo 120                     # modelo MILP
python tests/test_verificaciones.py                         # 18 verificaciones
```
