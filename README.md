# Análisis diferencial de Keccak con rondas dinámicas mediante MILP

[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Este repositorio contiene una implementación completa del análisis de seguridad diferencial de una variante de Keccak/SHA-3 con **número de rondas dinámico** (dependiente de intentos fallidos) y tamaños de palabra reducidos (`z = 4, 8`).

## Características

- **Núcleo SHA‑3** validado frente a `hashlib` (FIPS 202).
- **Cálculo correcto de la DDT** de la caja‑S χ (sin filtros espurios).
- **Modelo MILP** con reformulación **Convex Hull** de la DDT (sin Big‑M).
- **Certificación de optimalidad** para **R = 1..10** en tiempos ≤ 11 segundos.
- **Búsqueda heurística** de trayectorias para cotas superiores rápidas.
- **Paralelización** con `ProcessPoolExecutor` (flag `--workers`).
- **Estrategia de decisión** para certificar optimalidad.
- **Informe IEEE** en LaTeX incluido.

## Requisitos

- Python 3.8+
- Dependencias:
  - `pulp >= 2.7`
  - `highspy >= 1.5` (o instalar con `pip install pulp[highs]`)
  - (Opcional) `ortools` para comparación con CP‑SAT.

Instalación:
```bash
pip install -r requirements.txt
