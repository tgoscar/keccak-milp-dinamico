# Defectos detectados y corregidos

Durante el desarrollo se auditó una implementación previa del modelo. Se encontraron tres
defectos; **dos de ellos hacían el modelo infactible en las seis instancias**, de modo que
ningún resultado publicado a partir de aquella versión era válido.

Cada defecto tiene una verificación automatizada en
[`tests/test_verificaciones.py`](../tests/test_verificaciones.py), de manera que la
corrección sea auditable y no una simple afirmación.

---

## 1. Filtro espurio en el cálculo de la DDT

### Código defectuoso

```python
def precalcular_ddt():
    DDT = {}
    for a_in in range(32):
        posibles = []
        for a_out in range(32):
            count = 0
            for x in range(32):
                if (x ^ keccak_sbox_5bits(x)) == a_in:      # <-- ajeno a la DDT
                    if (keccak_sbox_5bits(x) ^ keccak_sbox_5bits(x ^ a_in)) == a_out:
                        count += 1
            if count > 0:
                posibles.append(a_out)
        DDT[a_in] = posibles
    return DDT
```

### Diagnóstico

La condición `(x ^ S(x)) == a_in` no forma parte de la definición de la Tabla de
Distribución de Diferencias. Restringe el recorrido de `x` a los puntos fijos de la
aplicación $x \mapsto x \oplus \chi(x)$, cuya imagen tiene sólo 11 elementos. En
consecuencia:

- **21 de las 32 filas quedan vacías** (todas aquellas cuyo $a$ no pertenece a esa imagen);
- se omiten **301 de las 317 transiciones válidas**;
- ejemplo concreto: $\mathrm{DDT}[1]$ resulta `[1, 9]` cuando lo correcto es
  `[1, 9, 17, 25]`.

Una fila vacía significa, para el modelo, que ninguna diferencia de salida es admisible
para esa diferencia de entrada; es decir, se declara imposible una transición que sí puede
ocurrir. Un modelo con 21 filas vacías no describe la primitiva.

### Corrección

```python
def ddt_transiciones():
    return {a: sorted({sbox_chi(x) ^ sbox_chi(x ^ a) for x in range(32)})
            for a in range(32)}
```

Se recorren todos los $x$ sin filtro. Resultado verificado: 317 transiciones válidas,
ninguna fila vacía, todas las filas con multiplicidad uniforme.

### Consecuencia adicional: la probabilidad diferencial máxima

Un material previo afirmaba que la probabilidad diferencial máxima de $\chi$ era
$2/4 = 2^{-1}$. La DDT correcta arroja multiplicidad máxima $8/32$ para $a \neq 0$, es
decir:

$$\mathrm{DP}_{\max} = 2^{-2}.$$

La diferencia no es menor: el peso de una trayectoria con $n$ cajas activas es $2^{-2n}$,
no $2^{-n}$, y los pares necesarios $2^{2n}$ en lugar de $2^{n}$. Para $R = 3$ con $z = 4$
($n = 18$) esto supone $2^{36} \approx 6.9 \times 10^{10}$ pares y no
$2^{18} \approx 2.6 \times 10^{5}$: cinco órdenes de magnitud de diferencia en la
estimación de seguridad.

**Verificación:** `test_ddt_transiciones_y_dp_maxima`,
`test_ddt_filtro_espurio_vacia_filas`, `test_ddt_filas_planas`.

---

## 2. *Gadget* XOR que vuelve el modelo infactible

### Código defectuoso

```python
aux = pulp.LpVariable(...)
prob += acum + D[...] - 2*aux <= 0
prob += acum - D[...] - 2*aux <= 0
prob += -acum + D[...] - 2*aux <= 0
prob += acum + D[...] - 2*aux >= 0
```

### Diagnóstico

La primera y la cuarta desigualdad juntas equivalen a la igualdad
$a + b - 2\,\mathrm{aux} = 0$, es decir $a + b = 2\,\mathrm{aux}$. Sobre variables binarias
esto sólo admite $a = b = \mathrm{aux} = 0$ o $a = b = \mathrm{aux} = 1$: **fuerza
$a = b$**.

El caso $a \neq b$ —precisamente aquel en que el XOR vale 1— resulta **infactible**. Como
$\theta$ requiere combinar bits distintos del estado, la infactibilidad se propaga a todo
el modelo. Comprobación exhaustiva de los cuatro casos:

| $a$ | $b$ | Resultado del solver | XOR esperado |
|:---:|:---:|---|:---:|
| 0 | 0 | Optimal, aux = 0 | 0 |
| 0 | 1 | **Infeasible** | 1 |
| 1 | 0 | **Infeasible** | 1 |
| 1 | 1 | Optimal, aux = 1 | 0 |

Ejecutado sobre el modelo completo, $R = 1$, $z = 4$ devuelve `Infeasible` en 0.2 s, de
modo que la función de resolución retornaba `None` en las seis instancias y los
«resultados esperados» documentados nunca se producían.

### Corrección

```python
def xor_lineal(prob, a, b, out):
    t = pulp.LpVariable(..., cat="Binary")
    prob += a + b - 2 * t == out
```

Una única variable auxiliar y una igualdad. Para $a + b \in \{0, 2\}$ se obtiene
$\mathrm{out} = 0$ con $t = a$; para $a + b = 1$ se obtiene $\mathrm{out} = 1$ con $t = 0$.
Exacto en los cuatro casos y más compacto que la versión defectuosa.

**Verificación:** `test_gadget_xor_exacto`, `test_gadget_xor_defectuoso_es_infactible`.

---

## 3. Condición de no trivialidad sesgada

### Código defectuoso

```python
prob += D[(0, 0, 0, 0)] == 1
```

### Diagnóstico

Fijar un bit concreto de la diferencia de entrada no sólo excluye la solución trivial:
restringe la búsqueda a las trayectorias que activan **ese** bit en particular. El mínimo
resultante es el mínimo condicionado a esa elección, que puede ser estrictamente mayor que
el mínimo global.

### Corrección

```python
prob += pulp.lpSum(D[(0, x, y, k)] for x in range(5)
                   for y in range(5) for k in range(z)) >= 1
```

Se exige al menos un bit de diferencia, sin especificar cuál. Cabe señalar que la solución
óptima para $R = 2$, $z = 4$ tiene **38 bits activos** en la diferencia de entrada: las
trayectorias de bajo peso son dispersas *después* de la capa lineal, no en la entrada
(puesto que $\lambda^{-1}$ es difusiva). Fijar un bit concreto habría descartado esa
trayectoria.

---

## 4. Un cuarto riesgo: el estado de terminación informado

No es un defecto del modelo, pero produce conclusiones erróneas con la misma facilidad.

La capa de modelado reporta `Optimal` cuando el solver termina, **incluso si sólo alcanzó
el límite de tiempo**. Un caso observado en este estudio:

```
R=2 z=4 | Optimal | S-boxes=9.0
```

mientras el registro del propio solver indicaba:

```
Status            Time limit reached
Primal bound      9
Dual bound        0
Gap               100% (tolerance: 0.01%)
```

Informar «9, óptimo» es doblemente incorrecto: no es óptimo, y de hecho existe una
trayectoria válida con 4 cajas activas. La incoherencia era además detectable a simple
vista, porque $z = 8$ (estado mayor) arrojaba 4: un estado más grande no puede admitir un
mínimo menor.

**Corrección adoptada:** consultar `highspy` directamente y distinguir `kOptimal` de
`kTimeLimit`, informando además la cota dual y el intervalo
[cota inferior, cota superior]. Véase `resolver_milp` en
[`src/milp_keccak.py`](../src/milp_keccak.py).

---

## Control de correctitud del núcleo

Como referencia independiente, la implementación de SHA-3-256 con los cinco pasos
separados se valida contra `hashlib` en 6 vectores de prueba, incluyendo los casos límite
de relleno (135, 136 y 137 bytes, con *rate* de 136). Si el núcleo no reprodujera FIPS 202,
cualquier análisis diferencial sobre él carecería de sentido.

**Verificación:** `test_sha3_contra_hashlib`.

---

## 5. Cotas por búsqueda voraz tomadas como definitivas

No es un defecto de código, sino un error de interpretación cometido en este trabajo y
corregido después. Se documenta porque afecta a las cifras publicadas en versiones previas
del informe.

Para $R \ge 2$, ante la imposibilidad de que el MILP cerrara el gap, se recurrió a una
búsqueda dirigida de trayectorias: en cada caja-S activa se elegía la diferencia de salida
válida de **menor peso de Hamming**. El resultado se reportó como cota superior del mínimo,
lo cual es correcto, pero se tomó implícitamente como una estimación cercana al óptimo, lo
cual no lo era:

| Caso | Cota voraz | Óptimo certificado (CP-SAT) |
|---|:---:|:---:|
| $R{=}3$, $z{=}4$ | 18 (2,2,14) | **9** (2,4,3) |
| $R{=}3$, $z{=}8$ | 20 (2,2,16) | $[8, 10]$ (2,4,4) |

El error de la cota es de un factor 2. La causa es instructiva: elegir el mínimo peso de
Hamming en cada caja es una decisión **local** que ignora su efecto sobre las rondas
siguientes. La trayectoria óptima hace lo contrario de lo que sugiere la intuición voraz:
sacrifica peso en la segunda ronda (4 cajas en lugar de 2) para poder **contraerse** en la
tercera (3 cajas en lugar de 14).

Consecuencia sobre las cifras de seguridad de $R = 3$, $z = 4$: el peso diferencial pasa de
$2^{-36}$ a $2^{-18}$, y los pares necesarios de $\sim 6.9 \times 10^{10}$ a
$\sim 2.6 \times 10^{5}$. Es decir, **la variante es considerablemente menos segura de lo
que indicaban las versiones previas del informe**. La conclusión cualitativa —la
resistencia crece con los intentos— se mantiene; su magnitud, no.

Lección metodológica: una cota superior no certificada no debe presentarse como una
aproximación del óptimo sin una estimación del error. En criptoanálisis diferencial esto es
especialmente delicado, porque la cifra de seguridad depende de la cota **inferior**, y una
cota superior floja puede dar una falsa sensación de holgura.

**Verificación:** `test_cotas_voraces_no_son_optimas`, `test_cpsat_certifica_optimos`.
