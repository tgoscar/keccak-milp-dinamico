# Análisis diferencial de una variante de Keccak/SHA-3 con rondas dinámicas

Modelos de **MILP** y **CP-SAT** para determinar el número mínimo de cajas-S activas en
una variante de la permutación Keccak (SHA-3) cuyo número de rondas se ajusta
dinámicamente según el contador de intentos fallidos de autenticación.

A partir de ese mínimo se acota la probabilidad de las trayectorias diferenciales y el
número de pares de textos que requeriría un ataque, lo que permite cuantificar si la
variante es efectivamente más segura y en qué magnitud.



---

## Resultados

Cinco de los seis casos tienen **óptimo certificado** (CP-SAT devuelve `OPTIMAL`, gap
cerrado):

| $z$ | $R$ | Mínimo de cajas-S | Por ronda | Estado | Probabilidad | Pares |
|:---:|:---:|:---:|:---:|---|:---:|:---:|
| 4 | 1 | **1** | 1 | certificado | $2^{-2}$ | $2^{2}$ |
| 4 | 2 | **4** | 2,2 | certificado | $2^{-8}$ | $2^{8}$ |
| 4 | 3 | **9** | 2,4,3 | certificado | $2^{-18}$ | $2^{18}$ |
| 8 | 1 | **1** | 1 | certificado | $2^{-2}$ | $2^{2}$ |
| 8 | 2 | **4** | 2,2 | certificado | $2^{-8}$ | $2^{8}$ |
| 8 | 3 | $[8, 10]$ | 2,4,4 | no cerrado | $\le 2^{-16}$ | $\ge 2^{16}$ |

Siendo $n$ el número de cajas-S activas, la probabilidad de la trayectoria se acota por
$2^{-2n}$ y los pares necesarios por $2^{2n}$, ya que la probabilidad diferencial máxima
de la caja-S $\chi$ es $2^{-2}$.

Las magnitudes de seguridad se derivan de la **cota inferior** del número de cajas
activas: es la que sustenta una garantía frente a un atacante. Para los casos
certificados la cota inferior coincide con el mínimo exacto.

**Lectura.** Al pasar de una a tres rondas, la complejidad de datos del ataque
diferencial asciende de $\sim 4$ pares a $\sim 2.6 \times 10^{5}$: la variable dinámica
endurece la primitiva cuando los intentos fallidos sugieren un ataque, si bien $R = 3$
sigue muy por debajo de un nivel de seguridad utilizable. Las salvedades están en
[Limitaciones](#limitaciones-del-estudio).

### MILP frente a CP-SAT

Ambos modelos describen la misma formulación diferencial, pero su comportamiento difiere
radicalmente:

| Caso | MILP (HiGHS, 120 s) | CP-SAT | Óptimo |
|---|:---:|:---:|:---:|
| $R{=}1$, $z{=}4$ | 1, certificado | 1, certificado (0,8 s) | **1** |
| $R{=}1$, $z{=}8$ | 1, certificado | 1, certificado (2,1 s) | **1** |
| $R{=}2$, $z{=}4$ | 9, gap 100 % | 4, certificado (2,9 s) | **4** |
| $R{=}2$, $z{=}8$ | 4, gap 100 % | 4, certificado (9,2 s) | **4** |
| $R{=}3$, $z{=}4$ | 18, gap 100 % | 9, certificado (68 s) | **9** |
| $R{=}3$, $z{=}8$ | sin solución | $[8, 10]$ (260 s) | ? |

El MILP no cierra el gap para $R \ge 2$: cerca del **90 % de sus variables** son las
binarias de selección *big-M* del encoding de la DDT (12 680 de 14 140 para $R{=}2$,
$z{=}4$), y las restricciones de paridad sobre $\mathbb{F}_2$ se relajan muy mal. CP-SAT
no necesita ese encoding —trata el XOR de forma nativa y $\chi$ como restricción de
tabla— y certifica en segundos lo que el MILP no cierra en minutos.

Se conservan ambos: el MILP es el modelo pedido por el enunciado y CP-SAT la herramienta
de certificación. La comparación es en sí misma un resultado metodológico.

---

## Estructura del repositorio

```
.
├── src/
│   ├── sha3_nucleo.py       Núcleo Keccak/SHA-3 (FIPS 202), variante dinámica,
│   │                        caja-S chi y su DDT
│   ├── milp_keccak.py       Modelo MILP y búsqueda dirigida de trayectorias
│   ├── cpsat_keccak.py      Modelo CP-SAT (certificación de los óptimos)
│   └── experimentos.py      Ejecución de los seis experimentos
├── tests/
│   └── test_verificaciones.py   15 verificaciones de las propiedades del informe
├── notebooks/
│   └── Keccak_MILP_analisis.ipynb   Recorrido explicado del análisis
├── docs/
│   ├── MODELO.md            Formulación de ambos modelos
│   ├── RESULTADOS.md        Resultados detallados y su interpretación
│   └── CORRECCIONES.md      Defectos detectados y corregidos
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

Dependencias: [PuLP](https://github.com/coin-or/pulp) y
[highspy](https://github.com/ERGO-Code/HiGHS) para el modelo MILP;
[OR-Tools](https://developers.google.com/optimization) para CP-SAT. Gurobi no se
requiere.

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
> distintos**:
>
> ```bash
> python src/cpsat_keccak.py     # sólo ortools
> python src/experimentos.py     # sólo highspy/pulp
> ```
>
> `src/milp_keccak.py` degrada automáticamente a CBC si detecta que `highspy` no puede
> cargarse, informando de que el estado de terminación deja de ser fiable. En el cuaderno,
> reinicia el kernel antes de la sección 8. Si sólo necesitas uno de los dos, instala sólo
> ese: `pip install ortools` o `pip install pulp highspy`.

---

## Uso

### Validar el núcleo y la DDT

```bash
python src/sha3_nucleo.py
```

```
SHA-3-256 frente a hashlib: 6/6 vectores coinciden (FIPS 202)
DDT de chi: 317 transiciones válidas, filas vacías: 0
Probabilidad diferencial máxima: 8/32 = 2^-2
DDT[1] = [1, 9, 17, 25]
```

### Certificar los óptimos (CP-SAT, recomendado)

```bash
python src/cpsat_keccak.py            # 120 s por caso
python src/cpsat_keccak.py 300        # más tiempo para R = 3, z = 8
```

### Ejecutar el modelo MILP

```bash
python src/experimentos.py --solo-trayectorias   # cotas por búsqueda dirigida
python src/experimentos.py --tiempo 120          # MILP completo
python src/experimentos.py --z 4 --rondas 3      # un caso concreto
```

### Verificaciones

```bash
python tests/test_verificaciones.py     # o: pytest -q
```

Comprueba los vectores FIPS 202, las propiedades de la DDT, la biyectividad de $\chi$,
la exactitud del *gadget* XOR, la invertibilidad de la capa lineal, la degeneración de
los desplazamientos de $\rho$, la certificación de CP-SAT y las limitaciones conocidas de
la búsqueda voraz.

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

- **Capas lineales** ($\theta, \rho, \pi$): XOR exactos. En MILP, $c = a \oplus b$ se
  codifica con una auxiliar binaria $t$ mediante $a + b - 2t = c$; en CP-SAT se usa
  `AddBoolXOr`. $\rho$ y $\pi$ son permutaciones y se imponen por reindexación, con los
  desplazamientos de $\rho$ reducidos módulo $z$.
- **Capa no lineal** ($\chi$): se exige que el par (diferencia de entrada, diferencia de
  salida) de cada caja-S sea una transición válida de la DDT. En MILP mediante variables
  de selección y *big-M*; en CP-SAT mediante una restricción de tabla con las 317
  transiciones.
- **Objetivo**: minimizar el total de cajas-S activas, con la diferencia de entrada no
  nula.
- **Ruptura de simetría**: la ronda sin $\iota$ conmuta con la traslación a lo largo de
  $z$, luego se puede exigir que el *slice* $z=0$ de la diferencia de entrada sea no
  nulo. El ahorro está acotado por $z$.

Detalle completo en [`docs/MODELO.md`](docs/MODELO.md).

---

## Limitaciones del estudio

1. **Un caso sin cerrar.** Para $R = 3$, $z = 8$ el mínimo está demostrado en
   $[8, 10]$: existe una trayectoria con 10 cajas activas y se ha probado que ninguna
   tiene menos de 8. Los cinco casos restantes tienen óptimo exacto certificado.

2. **Las versiones reducidas difunden peor que el estándar.** Al reducir los
   desplazamientos de $\rho$ módulo $z$, con $z=4$ los 25 offsets colapsan a 4 valores
   distintos y 7 carriles quedan sin rotar (con $z=8$: 8 valores, 3 sin rotar). Es
   inherente a la reducción de palabra y debe tenerse en cuenta al extrapolar a $w = 64$.

3. **La comparación es relativa.** Ninguna de las tres configuraciones se aproxima a las
   24 rondas de SHA-3. Con $R=3$ y $z=4$ el óptimo certificado de 9 cajas implica
   $\sim 2^{18}$ pares: suficiente para observar el efecto de la variable dinámica, pero
   muy lejos de un margen de seguridad utilizable. $R = 1$ es inservible ($\sim 4$ pares).

4. **Canal lateral temporal.** Al depender el número de rondas de un contador observable,
   el tiempo de cómputo filtra información sobre el estado de autenticación. Queda fuera
   del alcance de este análisis, pero es una consecuencia directa del diseño.

---

## Correcciones respecto de versiones previas

| Defecto | Efecto | Corrección |
|---|---|---|
| Filtro espurio en el cálculo de la DDT | 21 de 32 filas vacías, 301 de 317 transiciones omitidas | recorrer todos los $x$ sin filtrar |
| *Gadget* XOR con cuatro desigualdades | forzaba $a = b$: el modelo completo infactible | $a + b - 2t = c$ |
| No trivialidad fijando un bit concreto | sesgaba el mínimo | $\sum D_0 \ge 1$ |
| Estado de terminación informado por la capa de modelado | reportaba «Optimal» con gap del 100 % | consultar el estado real del solver |
| Cotas por búsqueda voraz tomadas como definitivas | 18 frente al óptimo real 9 ($R{=}3$, $z{=}4$) | certificar con CP-SAT |

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
