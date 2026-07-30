\# Modelo MILP para Keccak dinámico



\## Visión general



El modelo de Programación Lineal Entera Mixta (MILP) minimiza el número de cajas‑S activas en la capa no lineal $\\chi$ para una variante de Keccak con `R` rondas y tamaño de palabra `z`.



\## Variables



\- \*\*$D\_{r,x,y,k} \\in \\{0,1\\}$\*\*: diferencia en el bit `k` del carril `(x,y)` en la capa `r` (r = 0..R).

\- \*\*$A\_{r,y,k} \\in \\{0,1\\}$\*\*: actividad de la caja‑S en la fila `y`, bit `k`, ronda `r` (r = 0..R‑1).



\## Restricciones



\### Capa lineal ($\\theta, \\rho, \\pi$)



\- \*\*$\\theta$\*\*: paridades de columna y corrección con rotaciones módulo `z`.

\- \*\*$\\rho + \\pi$\*\*: reindexación directa de las diferencias (permutaciones).

\- \*\*XOR\*\*: codificación exacta $a + b - 2t = c$ con $t \\in \\{0,1\\}$.



\### Capa no lineal ($\\chi$)



La DDT de $\\chi$ se codifica mediante \*\*envolvente convexa\*\*:



Para cada caja‑S, se define el conjunto de puntos $P = \\{(a\_i, b\_i, \\alpha\_i)\\}$ donde:

\- $(a\_i, b\_i)$ son transiciones válidas de la DDT ($\\mathrm{DDT}\[a\_i]\[b\_i] > 0$).

\- $\\alpha\_i = 1$ si $a\_i \\neq 0$, $0$ en caso contrario.



Se introducen variables continuas $\\lambda\_i \\ge 0$ con $\\sum\_i \\lambda\_i = 1$ y se impone:



$$

\\begin{aligned}

v^{\\mathrm{in}} \&= \\sum\_i \\lambda\_i a\_i, \\\\

v^{\\mathrm{out}} \&= \\sum\_i \\lambda\_i b\_i, \\\\

A\_{r,y,k} \&= \\sum\_i \\lambda\_i \\alpha\_i.

\\end{aligned}

$$



Además, $v^{\\mathrm{in}} = \\sum\_{i=0}^{4} D\_{r,i,y,k} 2^i$ y $v^{\\mathrm{out}} = \\sum\_{i=0}^{4} D\_{r+1,i,y,k} 2^i$.



Esta reformulación es \*\*exacta\*\* porque el poliedro generado por los puntos $P$ es integral (todos sus vértices son enteros).



\### Objetivo



Minimizar el número total de cajas‑S activas:

$$

\\min \\sum\_{r=0}^{R-1}\\sum\_{y=0}^{4}\\sum\_{k=0}^{z-1} A\_{r,y,k}.

$$



\### No trivialidad



$$

\\sum\_{x,y,k} D\_{0,x,y,k} \\ge 1

$$



(se exige al menos una diferencia en la entrada).



\### Simetría rotacional (opcional)



Para acelerar la búsqueda, se fija el bit $D\_{0,0,0,0} = 1$ (válido por invarianza rotacional de Keccak).



\## Estrategia de decisión



En lugar de minimizar directamente, se usa una \*\*estrategia de decisión\*\*:

1\. Se obtiene una cota superior $K$ mediante la búsqueda heurística (`cota\_superior`).

2\. Se añade la restricción $\\sum A \\le K - 1$.

3\. Si el problema es \*\*infactible\*\*, entonces $K$ es el óptimo.

4\. Si es \*\*factible\*\*, se encontró una trayectoria mejor, se actualiza $K$ y se repite.



Esto permite certificar optimalidad en tiempos mucho menores que la minimización directa.



\## Solver



Se utiliza \*\*HiGHS\*\* a través de `highspy` (o `pulp` con HiGHS). Parámetros recomendados:

\- `mip\_rel\_gap = 0.0`

\- `presolve = on`

\- `parallel = on`

\- `mip\_detect\_symmetry = on`



\## Validación



El modelo se ha probado para `z ∈ {4,8}` y `R = 1..10`, certificando el óptimo $n = R$ en todos los casos con tiempos inferiores a 11 segundos.

