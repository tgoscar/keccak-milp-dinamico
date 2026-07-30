# Análisis diferencial de Keccak con rondas dinámicas mediante MILP y CP-SAT

Modelos de programación matemática para determinar el número mínimo de cajas-S activas
en una variante de la permutación Keccak (SHA-3) cuyo número de rondas se ajusta
dinámicamente según el contador de intentos fallidos de autenticación.

A partir de ese mínimo se acota la probabilidad de las trayectorias diferenciales y el
número de pares de textos que requeriría un ataque, lo que permite cuantificar si la
variante es efectivamente más segura y en qué magnitud.

## Contenido

- **Núcleo SHA-3** validado frente a `hashlib` (FIPS 202).
- **DDT de la caja-S $\chi$** calculada correctamente: 317 transiciones,
  $\mathrm{DP}_{\max} = 2^{-2}$.
- **Modelo CP-SAT** con XOR nativo y $\chi$ como restricción de tabla. Es el motor de
  certificación, con soporte para $R$ arbitrario y paralelismo.
- **Modelo MILP** en dos variantes: *big-M* y envolvente convexa de la DDT.
- **Consolidación de cotas por monotonía**, que propaga las cotas inferiores a rondas
  superiores.
- **18 verificaciones automatizadas** de todas las propiedades afirmadas.
- **Informe en formato IEEE** (LaTeX y PDF).

---

## Resultados

### $z = 4$ (Keccak-f[100], 20 cajas-S por ronda)

| $R$ | Mínimo | Por ronda | Estado | Pares garantizados |
|:---:|:---:|:---:|---|:---:|
| 1 | **1** (exacto) | 1 | certificado | $2^{2}$ |
| 2 | **4** (exacto) | 2,2 | certificado | $2^{8}$ |
| 3 | **9** (exacto) | 2,4,3 | certificado | $2^{18}$ |
| 4 | $[9, 14]$ | — | cota | $2^{18}$ |
| 5 | $[9, 35]$ | — | cota | $2^{18}$ |
| 6 | $[9, 61]$ | — | cota | $2^{18}$ |
| 7 | $[11, 64]$ | — | cota | $2^{22}$ |
| 8 | $[13, 86]$ | — | cota | $2^{26}$ |
| 9 | $[13, 103]$ | — | cota | $2^{26}$ |
| 10 | $[13, 111]$ | — | cota | $2^{26}$ |

### $z = 8$ (Keccak-f[200], 40 cajas-S por ronda)

| $R$ | Mínimo | Por ronda | Estado | Pares garantizados |
|:---:|:---:|:---:|---|:---:|
| 1 | **1** (exacto) | 1 | certificado | $2^{2}$ |
| 2 | **4** (exacto) | 2,2 | certificado | $2^{8}$ |
| 3 | **10** (exacto) | 2,4,4 | certificado | $2^{20}$ |
| 4–10 | $\ge 10$ | — | cota | $2^{20}$ |

Los casos $R \le 3$ tienen **óptimo exacto certificado**. Para $R \ge 4$ se reportan
intervalos: el extremo inferior está demostrado, el superior es la mejor trayectoria
hallada con un presupuesto de unos 35 s por caso, y mejora con más tiempo.

Tres puntos de lectura, detallados en [`docs/RESULTADOS.md`](docs/RESULTADOS.md):

**Los pares garantizados se derivan de la cota inferior**, no de la superior. Si el
mínimo real fuera menor que el hallado, existiría un ataque mejor; sólo la cota inferior
sustenta una garantía.

**Las cotas inferiores se propagan por monotonía.** Como $m(R)$ es no decreciente en $R$,
el óptimo certificado $m(3) = 9$ garantiza $m(R) \ge 9$ para todo $R \ge 3$ sin resolver
esos casos.

**Las cotas superiores de $R \ge 4$ son flojas** y no deben leerse como estimaciones del
óptimo. En $R = 3$ este proyecto llegó a manejar una cota de 18 cuando el óptimo real
era 9.

### MILP frente a CP-SAT

| Caso | MILP *big-M* (120 s) | MILP envolvente convexa | CP-SAT | Óptimo |
|---|:---:|:---:|:---:|:---:|
| $R{=}1$, $z{=}4$ | 1, certificado | 1 | 1, certificado (0,2 s) | **1** |
| $R{=}2$, $z{=}4$ | 9, gap 100 % | 7, sin certificar | 4, certificado (0,7 s) | **4** |
| $R{=}3$, $z{=}4$ | 18, gap 100 % | — | 9, certificado (7,2 s) | **9** |
| $R{=}3$, $z{=}8$ | sin solución | — | 10, certificado (118 s) | **10** |

El MILP no cierra el gap para $R \ge 2$: cerca del 90 % de sus variables son las binarias
de selección *big-M* del encoding de la DDT, y las restricciones de paridad sobre
$\mathbb{F}_2$ se relajan muy mal. La envolvente convexa elimina la primera causa y
mejora la cota, pero tampoco certifica. CP-SAT trata ambas estructuras de forma nativa.

Se conservan los tres modelos: la comparación es en sí misma un resultado metodológico.

---

## Estructura del repositorio

```
.
├── src/
│   ├── sha3_nucleo.py       Núcleo Keccak/SHA-3 (FIPS 202), variante dinámica,
│   │                        caja-S chi y su DDT
│   ├── cpsat_keccak.py      Modelo CP-SAT: certificación, R arbitrario, paralelismo
│   ├── milp_keccak.py       Modelo MILP (big-M y envolvente convexa) y búsqueda
│   │                        dirigida de trayectorias
│   └── experimentos.py      Ejecución de los experimentos con el modelo MILP
├── tests/
│   └── test_verificaciones.py   18 verificaciones
├── notebooks/
│   └── Keccak_MILP_analisis.ipynb   Recorrido explicado del análisis
├── docs/
│   ├── MODELO.md            Formulación de los tres modelos y la monotonía
│   ├── RESULTADOS.md        Resultados detallados e interpretación
│   └── CORRECCIONES.md      Seis defectos detectados y corregidos
├── informe/
│   ├── reporte_ieee.tex     Informe en formato IEEE (español)
│   └── reporte_ieee.pdf
├── requirements.txt
└── LICENSE
```

---

## Instalación

Requiere Python 3.9 o superior.

```bash
git clone https://github.com/tgoscar/keccak-milp-dinamico.git
cd keccak-milp-dinamico
pip install -r requirements.txt
```

> **Incompatibilidad conocida: `highspy` y `ortools` no pueden convivir en un mismo
> proceso de Python.** Ambos empaquetan builds distintos de HiGHS y sus símbolos
> colisionan, en cualquier orden de importación:
>
> ```
> ImportError: libortools.so.9: undefined symbol: _Z19setLocalOptionValue...
> ImportError: highspy/_core...so: undefined symbol: _ZN5Highs13releaseMemoryEv
> ```
>
> Por eso los modelos viven en módulos separados y **deben ejecutarse en procesos
> distintos**: `python src/cpsat_keccak.py` (sólo ortools) y
> `python src/experimentos.py` (sólo highspy/pulp). `src/milp_keccak.py` degrada
> automáticamente a CBC si detecta que `highspy` no puede cargarse. En el cuaderno,
> reinicia el kernel antes de la sección de CP-SAT. Si sólo necesitas uno,
> instala sólo ese: `pip install ortools` o `pip install pulp highspy`.

---

## Uso

### Certificación con CP-SAT

```bash
python src/cpsat_keccak.py                                   # R=1..3, z=4 y 8
python src/cpsat_keccak.py --rondas 1-10 --tiempo 300 --workers 4
python src/cpsat_keccak.py --rondas 3 --z 8 --tiempo 1800 --hilos 8
python src/cpsat_keccak.py --rondas 1-10 --salida resultados.json
```

`--workers` reparte los casos entre procesos; `--hilos` son los hilos internos que
CP-SAT usa en cada caso. Sin `--hilos`, los núcleos disponibles se reparten entre los
procesos. Ten presente que CP-SAT **no es determinista** al variar los hilos: para casos
difíciles conviene reintentar con configuraciones distintas.

### Modelo MILP

```bash
python src/experimentos.py --solo-trayectorias   # cotas por búsqueda dirigida
python src/experimentos.py --tiempo 120          # MILP completo
python src/experimentos.py --z 4 --rondas 3      # un caso concreto
```

### Verificaciones

```bash
python tests/test_verificaciones.py     # o: pytest -q
```

Comprueba los vectores FIPS 202, las propiedades de la DDT, la biyectividad de $\chi$, la
exactitud del *gadget* XOR, la invertibilidad de la capa lineal, la degeneración de los
desplazamientos de $\rho$, la certificación de CP-SAT, la consolidación por monotonía y
las dos formulaciones de la envolvente convexa.

### Validar el núcleo

```bash
python src/sha3_nucleo.py
```

### Usar la variante desde código

```python
from src.sha3_nucleo import keccak_modificado, intentos_a_rondas

intentos_a_rondas(5)    # 1 ronda
intentos_a_rondas(15)   # 2 rondas
intentos_a_rondas(25)   # 3 rondas

digest = keccak_modificado(b"mensaje", intentos=15)
```

---

## Los modelos en breve

El estado de Keccak es $5 \times 5 \times z$ y cada ronda aplica
$\theta, \rho, \pi, \chi$ ($\iota$ se omite: sólo suma una constante y no afecta a las
diferencias). Se asigna una variable binaria a cada bit de diferencia del estado.

- **Capas lineales**: XOR exactos. En MILP, $c = a \oplus b$ se codifica con una auxiliar
  binaria mediante $a + b - 2t = c$; en CP-SAT se usa `AddBoolXOr`. $\rho$ y $\pi$ son
  permutaciones y se imponen por reindexación, con los desplazamientos de $\rho$
  reducidos módulo $z$.
- **Capa no lineal**: el par (diferencia de entrada, diferencia de salida) de cada caja-S
  debe ser una transición válida de la DDT. En CP-SAT mediante `AddAllowedAssignments`;
  en MILP mediante *big-M* o envolvente convexa **impuesta componente a componente**
  (colapsarla en escalares invalida el modelo, véase `docs/CORRECCIONES.md`, defecto 6).
- **Objetivo**: minimizar el total de cajas-S activas, con la diferencia de entrada no
  nula.
- **Ruptura de simetría**: la ronda sin $\iota$ conmuta con la traslación a lo largo de
  $z$, luego puede exigirse que el *slice* $z=0$ de la diferencia de entrada sea no nulo.

Detalle completo en [`docs/MODELO.md`](docs/MODELO.md).

---

## Limitaciones

1. **Sólo $R \le 3$ está certificado.** Para $R \ge 4$ se dispone de cotas; el mínimo
   real está en el intervalo indicado.
2. **Las versiones reducidas difunden peor que el estándar.** Al reducir los
   desplazamientos de $\rho$ módulo $z$, con $z=4$ los 25 offsets colapsan a 4 valores y
   7 carriles quedan sin rotar (con $z=8$: 8 valores, 3 sin rotar). Debe tenerse en
   cuenta al extrapolar a $w = 64$.
3. **La comparación es relativa.** Ninguna configuración se aproxima a las 24 rondas de
   SHA-3. Con $R=3$ y $z=4$, el óptimo de 9 cajas implica $\sim 2^{18}$ pares: suficiente
   para observar el efecto de la variable dinámica, pero lejos de un margen utilizable.
   $R = 1$ es inservible ($\sim 4$ pares).
4. **Canal lateral temporal.** Al depender el número de rondas de un contador observable,
   el tiempo de cómputo filtra información sobre el estado de autenticación. Queda fuera
   del alcance de este análisis, pero es consecuencia directa del diseño.

---

## Defectos detectados y corregidos

| Defecto | Efecto | Corrección |
|---|---|---|
| Filtro espurio en el cálculo de la DDT | 21 de 32 filas vacías, 301 de 317 transiciones omitidas | recorrer todos los $x$ sin filtrar |
| *Gadget* XOR con cuatro desigualdades | forzaba $a = b$: modelo infactible | $a + b - 2t = c$ |
| No trivialidad fijando un bit concreto | sesgaba el mínimo | $\sum D_0 \ge 1$ |
| Estado de terminación informado por la capa de modelado | reportaba «Optimal» con gap del 100 % | consultar el estado real del solver |
| Cotas por búsqueda voraz tomadas como definitivas | 18 frente al óptimo real 9 | certificar con CP-SAT |
| Envolvente convexa colapsada en escalares | restricción de la DDT vacía; devolvía $R$ para todo $R$ | imponerla componente a componente |

Análisis completo en [`docs/CORRECCIONES.md`](docs/CORRECCIONES.md).

---

## Referencias

1. NIST, *SHA-3 Standard: Permutation-Based Hash and Extendable-Output Functions*,
   FIPS PUB 202, 2015.
2. G. Bertoni, J. Daemen, M. Peeters, G. Van Assche, *The Keccak Reference*, v3.0, 2011.
3. N. Mouha, Q. Wang, D. Gu, B. Preneel, «Differential and Linear Cryptanalysis Using
   Mixed-Integer Linear Programming», *Inscrypt*, 2011.
4. J. Daemen, G. Van Assche, «Differential Propagation Analysis of Keccak», *FSE*, 2012.
5. S. Li, X. Dong, Z. Liu, «New MILP Modeling: Improved Conditional Cube Attacks on
   Keccak-based Constructions», IACR ePrint 2017/1030.
6. L. Perron, F. Didier, *CP-SAT*, Google OR-Tools.

---

## Licencia

MIT. Véase [`LICENSE`](LICENSE).
