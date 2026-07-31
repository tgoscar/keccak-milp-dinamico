# Criptoanálisis lineal y diferencial de una variante ligera de Keccak

Segunda parte del trabajo: criptoanálisis lineal (LAT), rediseño de Keccak para
hardware limitado, y verificación formal de los criterios de diseño resultantes.

Parte del análisis diferencial MILP/CP-SAT desarrollado en la rama `main`.

## Resultados principales

**La premisa de partida no se sostiene.** Keccak-f[200] alcanza el 80 % de difusión
en **2 rondas** y Keccak-f[800] en 3. Las 24 rondas del estándar responden a margen
de seguridad, no a difusión.

| Construcción | Estado | Rondas hasta 80 % |
|---|:---:|:---:|
| Keccak-f[200] | 200 bits | **2** |
| Keccak-f[400] | 400 bits | **2** |
| Keccak-f[800] | 800 bits | **3** |
| Variante (Σ por carril) | 200 bits | 3 |

**Cotas diferenciales certificadas** (200 bits, 40 cajas-S por ronda):

| Construcción | R=1 | R=2 | R=3 |
|---|:---:|:---:|:---:|
| Keccak-f[200] | 1 | 4 | **10** |
| Variante (k=3) | 1 | 4 | 7 |

La variante alcanza en 4 rondas la garantía que Keccak da en 3.

**La capa no lineal parcial es insegura.** Admite un subespacio invariante de
2⁸ = 256 diferencias que atraviesan *cualquier* número de rondas con probabilidad 1.

| Capa no lineal | R=1 | R=2 | R=3 | R=6 |
|---|:---:|:---:|:---:|:---:|
| χ completo | 0 | 0 | 0 | 0 |
| χ parcial (1/2) | 2¹⁰⁰ | 2⁸ | 2⁸ | 2⁸ |

## Cuatro criterios de diseño verificados

1. **Número impar de términos por carril.** Con k par, el estado todo-unos se
   cancela consigo mismo y la capa lineal es singular.
2. **La paridad debe cubrir 4 carriles, no 5**, por el mismo argumento.
3. **El núcleo debe ser trivial**: ninguna diferencia debe esquivar la capa no
   lineal. Es un sistema lineal cuyo núcleo se calcula por eliminación gaussiana.
4. **La cota diferencial la fija el número de términos, no las constantes.** Una
   búsqueda sobre 17 125 conjuntos de constantes lo confirma: los mejores según
   dispersión dieron 5 cajas activas, peor que las 7 de partida.

## Estructura

```
src/
  lat.py               Ejercicios 1-3: aproximaciones afines, LAT, DDT
  difusion.py          Matriz de dependencia y medida de difusión
  kernel.py            Cálculo del subespacio que esquiva la no linealidad
  modelo_variante.py   Modelo CP-SAT de la variante (χ completo o parcial)
  variante_k.py        Familia parametrizada con k rotaciones por carril
informe/
  reporte_criptoanalisis.tex / .pdf
```

## Uso

```bash
python src/lat.py           # LAT y DDT de la caja-S de 4 bits
python src/difusion.py      # difusión de Keccak-f[200]
python src/kernel.py        # dimensión del núcleo según la fracción de χ
python src/modelo_variante.py   # cotas diferenciales con CP-SAT
```

Requiere `ortools`. Ejecutar en un proceso separado de los módulos que usan
`highspy` (véase la incompatibilidad documentada en la rama `main`).
