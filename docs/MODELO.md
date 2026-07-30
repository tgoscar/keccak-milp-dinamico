# Formulación del modelo MILP

## 1. La variante analizada

El número de rondas depende del contador de intentos fallidos $I$:

$$R(I) = \begin{cases} 1 & I < 10 \\ 2 & 10 \le I < 20 \\ 3 & 20 \le I < 30 \end{cases}$$

Cada ronda aplica $\theta, \rho, \pi, \chi$. El paso $\iota$ se omite del modelo
diferencial: suma una constante fija, y al considerar la diferencia entre dos estados esa
constante se cancela.

El tamaño de palabra se reduce a $z \in \{4, 8\}$, dando Keccak-f[100] y Keccak-f[200]
—estados de 100 y 200 bits— con $5z$ cajas-S por ronda (20 y 40 respectivamente).

## 2. Variables

| Variable | Dominio | Significado |
|---|---|---|
| $D_{r,x,y,k}$ | $\{0,1\}$ | diferencia en el bit $k$ del carril $(x,y)$ en la ronda $r$; $r = 0,\dots,R$ |
| $A_{r,y,k}$ | $\{0,1\}$ | la caja-S de la fila $y$, bit $k$, está activa en la ronda $r$ |
| $t_{\bullet}$ | $\{0,1\}$ | auxiliares de los XOR |
| $v^{\mathrm{in}}, v^{\mathrm{out}}$ | $\{0,\dots,31\}$ | valor entero de las diferencias de entrada/salida de cada caja-S |
| $s_j$ | $\{0,1\}$ | selección de la transición $j$ de la DDT |

Se requieren $R+1$ capas de estado: la entrada más la salida de cada ronda.

## 3. Codificación exacta del XOR

Todas las capas lineales se reducen a XOR. Sobre los enteros, $c = a \oplus b$ se codifica
con una variable auxiliar binaria $t$:

$$a + b - 2t = c$$

Verificación de los cuatro casos: si $a + b = 0$ entonces $t = 0, c = 0$; si $a + b = 1$
entonces $t = 0, c = 1$; si $a + b = 2$ entonces $t = 1, c = 0$. Es exacta y usa una sola
variable auxiliar por operación.

Para el XOR de $m$ entradas se encadenan $m-1$ aplicaciones. (Una alternativa más compacta
es la restricción de paridad $\sum_i x_i - c = 2d$ con $d \in \{0, \dots, \lfloor m/2
\rfloor\}$ entera, que emplea una sola auxiliar por XOR múltiple.)

## 4. Capa lineal $\lambda = \pi \circ \rho \circ \theta$

### $\theta$

Paridades de columna, corrección y aplicación:

$$C_{x,k} = \bigoplus_{y=0}^{4} D_{r,x,y,k}$$

$$D^{\theta}_{x,k} = C_{x-1,k} \oplus C_{x+1,\,(k-1) \bmod z}$$

$$\widetilde{D}_{x,y,k} = D_{r,x,y,k} \oplus D^{\theta}_{x,k}$$

con $x$ módulo 5. El desplazamiento de 1 posición dentro de $\theta$ se toma módulo $z$.
Cada bit de salida depende de 11 bits de entrada, lo que hace de $\theta$ la principal
fuente de difusión.

### $\rho$ y $\pi$

Son permutaciones de posición: no requieren variables nuevas, sólo reindexación.

$$D^{\rho\pi}_{y,\;(2x+3y) \bmod 5,\;(k + \mathrm{rot}_z[x][y]) \bmod z} = \widetilde{D}_{x,y,k}$$

donde $\mathrm{rot}_z[x][y] = \mathrm{ROT}[x][y] \bmod z$, con $\mathrm{ROT}$ la tabla
estándar de desplazamientos de Keccak (especificada para $w = 64$).

### Invertibilidad

Construyendo explícitamente la matriz de $\lambda$ sobre $\mathbb{F}_2$ y calculando su
rango por eliminación gaussiana se obtiene rango pleno: 100 para $z=4$ y 200 para $z=8$.
Por tanto $\lambda$ es biyectiva.

Combinado con la biyectividad de $\chi$, la ronda completa es una biyección, de donde se
sigue la **cota inferior estructural**: una diferencia no nula no puede anularse, cada
ronda aporta al menos una caja-S activa, y el mínimo para $R$ rondas es $\ge R$.

## 5. Capa no lineal $\chi$

$$\chi_i = x_i \oplus (\overline{x_{i+1}} \cdot x_{i+2}), \qquad i = 0,\dots,4$$

Opera sobre filas de 5 bits a lo largo de $x$; hay $5z$ instancias por ronda. Es una
permutación de grado algebraico 2.

Se evaluaron dos linealizaciones.

### 5.1 Modelo de $\chi$-imagen determinista (descartado)

Reproduce $\chi$ sobre la diferencia con compuertas auxiliares. La compuerta
$g = \overline{a} \cdot b$ se impone con:

$$g \le 1 - a, \qquad g \le b, \qquad g \ge b - a$$

y luego $\mathrm{out}_x = \mathrm{in}_x \oplus g$.

Es compacto, pero sólo modela **una** trayectoria: la imagen determinista de $\chi$ sobre
la diferencia. Mide difusión y **sobreestima** el número de cajas activas en varias
rondas, por lo que no constituye una cota diferencial rigurosa.

Advertencia: si se omite la negación —usando $g = a \cdot b$— la aplicación deja de ser
biyectiva. La diferencia $(1,1,1,1,1)$ colapsa a $(0,0,0,0,0)$, permitiendo que las
diferencias «desaparezcan», lo que produce resultados sin sentido como rondas con cero
cajas activas.

### 5.2 Modelo basado en DDT (adoptado)

Se exige que el par (diferencia de entrada, diferencia de salida) de cada caja-S sea una
transición válida, dejando al atacante la elección entre todas las posibles. Codificando

$$v^{\mathrm{in}} = \sum_{i=0}^{4} D^{\rho\pi}_{i,y,k} 2^{i}, \qquad
v^{\mathrm{out}} = \sum_{i=0}^{4} D_{r+1,i,y,k} 2^{i}$$

y siendo $\mathcal{P} = \{(a,b) : \mathrm{DDT}[a][b] > 0\}$ el conjunto de transiciones
válidas ($|\mathcal{P}| = 317$), con $M = 31$:

$$v^{\mathrm{in}} - a_j \le (1-s_j)M, \qquad a_j - v^{\mathrm{in}} \le (1-s_j)M$$

$$v^{\mathrm{out}} - b_j \le (1-s_j)M, \qquad b_j - v^{\mathrm{out}} \le (1-s_j)M$$

$$\sum_{j} s_j = 1$$

Es decir, exactamente una transición se selecciona, y las restricciones *big-M* fuerzan que
los valores coincidan con ella. Este es el modelo diferencial riguroso: sus mínimos son los
mínimos diferenciales reales.

### 5.3 Actividad

$$v^{\mathrm{in}} \le M \cdot A_{r,y,k}, \qquad A_{r,y,k} \le v^{\mathrm{in}}$$

La primera fuerza $A = 1$ cuando $v^{\mathrm{in}} > 0$; la segunda fuerza $A = 0$ cuando
$v^{\mathrm{in}} = 0$. Por tanto $A_{r,y,k} = 1 \iff v^{\mathrm{in}} \neq 0$.

## 6. Objetivo y no trivialidad

$$\min \sum_{r=0}^{R-1} \sum_{y=0}^{4} \sum_{k=0}^{z-1} A_{r,y,k}$$

sujeto a

$$\sum_{x,y,k} D_{0,x,y,k} \ \ge\ 1$$

Esta última excluye la solución trivial exigiendo al menos un bit de diferencia en la
entrada, **sin fijar cuál**. Fijar un bit concreto restringiría la búsqueda a las
trayectorias que activan ese bit y podría dar un mínimo estrictamente mayor que el global.

## 7. Búsqueda dirigida de trayectorias

Para $R \ge 2$ el MILP no cierra el gap, y sus heurísticas genéricas producen incumbentes
pobres o ninguno. Se complementa con una búsqueda dirigida que construye trayectorias
diferenciales válidas y devuelve la de menor número de cajas activas (una cota superior
del mínimo).

Procedimiento:

1. **Arranques.** Diferencias post-lineales dispersas: un solo bit activo, y pares de bits
   en la misma fila de $\chi$ —estos últimos reproducen las trayectorias de bajo peso
   características de Keccak.
2. **Propagación.** En cada caja-S activa se elige la diferencia de salida válida de menor
   peso de Hamming según la DDT, y se aplica $\lambda$ para la ronda siguiente.
3. **Mejora local.** Se voltea un bit a la vez de la diferencia inicial mientras el número
   de cajas activas disminuya.

La búsqueda opera en el espacio **post-lineal** deliberadamente: las trayectorias de bajo
peso son dispersas después de $\lambda$, no en la entrada, ya que $\lambda^{-1}$ es
difusiva. Buscar diferencias de entrada dispersas —el error inicial de este trabajo—
converge a valores mucho peores.

## 8. Nota sobre la relajación lineal

Para $R \ge 2$ la cota dual del solver permanece en 0 durante toda la ejecución, con gap
del 100 %. La causa: las restricciones de paridad sobre $\mathbb{F}_2$ se relajan muy mal.
Al permitir valores fraccionarios, el problema continuo se satisface con objetivo cercano a
cero mediante asignaciones que no corresponden a ninguna diferencia real 0/1. El
*branch-and-bound* tendría que enumerar un árbol inabordable para elevar esa cota.

Es una limitación conocida del MILP a nivel de bit frente a primitivas con difusión XOR
intensa, no una propiedad de la primitiva analizada. Certificar el óptimo para $R \ge 2$
exige formulaciones SAT/SMT —los solvers SAT tratan el XOR de forma nativa— o herramientas
que exploten la estructura del *column-parity kernel* de Keccak.

---

## 9. El modelo CP-SAT

La misma formulación se implementó en CP-SAT (OR-Tools), que dispone de restricciones
nativas para las dos estructuras que el MILP codifica con dificultad.

### 9.1 XOR nativo

CP-SAT admite `AddBoolXOr(literales)`, que impone que la paridad de los literales sea
impar. Para expresar $c = \bigoplus_i x_i$ basta negar un literal:

```python
m.AddBoolXOr([c.Not()] + [x_i for i in ...])
```

No se requieren variables auxiliares $t$, y el solver razona sobre la estructura XOR
mediante propagación y aprendizaje de cláusulas.

### 9.2 $\chi$ como restricción de tabla

En lugar del encoding *big-M* con 317 variables de selección por caja-S, se declara
directamente el conjunto de tuplas admisibles:

```python
m.AddAllowedAssignments(bits_entrada + bits_salida, TUPLAS)
```

donde `TUPLAS` contiene las 317 transiciones válidas de la DDT como vectores de 10 bits.
Esto elimina cerca del 90 % de las variables del modelo MILP.

### 9.3 Actividad y objetivo

$$A_{r,y,k} = \max(\Delta^{\mathrm{in}}_0, \dots, \Delta^{\mathrm{in}}_4)$$

impuesto con `AddMaxEquality`, que sobre booleanos equivale al OR. El objetivo es el mismo:
minimizar $\sum A$.

### 9.4 Ruptura de simetría

La ronda sin $\iota$ conmuta con la traslación a lo largo de $z$:

- $\theta$ rota una posición en $z$ (conmuta con cualquier traslación en $z$);
- $\rho$ traslada cada carril una constante (conmuta);
- $\pi$ permuta carriles sin depender de $z$;
- $\chi$ actúa dentro de cada *slice*.

Precisamente por eso existe $\iota$: para romper esa invariancia en la permutación
completa. Como aquí se omite, toda trayectoria pertenece a una órbita de tamaño divisor de
$z$, y puede exigirse sin pérdida de generalidad que el *slice* $z=0$ de la diferencia de
entrada sea no nulo:

$$\sum_{x,y} D_{0,x,y,0} \ \ge\ 1$$

Esta restricción implica también la no trivialidad. Conviene advertir que una formulación
del tipo «el primer bit no nulo está en $k=0$» requiere un esquema *lex-leader* correcto;
hecha a la ligera puede eliminar el óptimo. El ahorro está acotado por $z$ (4 u 8 veces),
de modo que es útil pero modesto.

### 9.5 Una nota sobre el encoding de $\chi$ en MILP

Si se deseara mantener el MILP, existe una vía para reducir drásticamente su tamaño: para
cada diferencia de entrada $a \neq 0$, el conjunto de diferencias de salida válidas de
$\chi$ es un **subespacio afín** de $\mathbb{F}_2^5$, de dimensión 2, 3 o 4 (verificado
exhaustivamente en `test_ddt_filas_planas` y comprobable con `ddt_transiciones`). Esa
estructura permite sustituir las 317 variables de selección por un puñado de ecuaciones
lineales. No se implementó porque CP-SAT ya resuelve el problema, pero es la mejora natural
para quien quiera insistir con MILP.
