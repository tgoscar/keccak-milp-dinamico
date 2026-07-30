#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Certificación del número mínimo de cajas-S activas en Keccak con rondas
dinámicas, mediante CP-SAT (OR-Tools). Soporta R arbitrario y paralelismo.

Por qué CP-SAT y no MILP
------------------------
La formulación diferencial necesita dos estructuras que el MILP codifica mal:
el XOR sobre GF(2) y la tabla de transiciones válidas de chi. CP-SAT las trata
de forma nativa (AddBoolXOr y AddAllowedAssignments), sin variables auxiliares
ni big-M. En la práctica certifica en segundos casos que el MILP no cierra en
minutos. Véase docs/MODELO.md y docs/RESULTADOS.md.

Consolidación de cotas por monotonía
------------------------------------
El mínimo m(R) es no decreciente en R: truncar una trayectoria de R+1 rondas
produce una de R rondas con un número de cajas activas menor o igual. Por tanto

    m(R) >= max{ LB(r) : r <= R }

y una cota inferior demostrada para pocas rondas se propaga a todas las
mayores. Esto importa porque la garantía de seguridad depende de la cota
INFERIOR: si m(R) >= n, un ataque diferencial necesita del orden de 2^(2n)
pares.

Uso
---
    python src/cpsat_keccak.py                          # R=1..3, z=4 y 8
    python src/cpsat_keccak.py --rondas 1-10 --tiempo 300 --workers 4
    python src/cpsat_keccak.py --rondas 3 --z 8 --tiempo 1800 --hilos 8
    python src/cpsat_keccak.py --rondas 1-10 --salida resultados.json

--workers reparte los casos entre procesos; --hilos son los hilos internos que
CP-SAT usa en cada caso. Sin --hilos, se reparten los núcleos disponibles.
"""

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ortools.sat.python import cp_model

try:
    from sha3_nucleo import ROT, ddt_transiciones
except ImportError:
    from src.sha3_nucleo import ROT, ddt_transiciones

L = 5
DP_MAX_LOG2 = -2          # probabilidad diferencial máxima de chi: 2^-2

# Transiciones válidas de la DDT como tuplas de 10 bits (5 de entrada, 5 de salida)
DDT = ddt_transiciones()
TUPLAS = []
for _a, _bl in DDT.items():
    for _b in _bl:
        TUPLAS.append([(_a >> i) & 1 for i in range(5)] +
                      [(_b >> i) & 1 for i in range(5)])


# ===========================================================================
# Modelo
# ===========================================================================
def construir(R, z, simetria=True):
    """Modelo CP-SAT de propagación diferencial para R rondas y palabra de z bits."""
    m = cp_model.CpModel()
    D = {(r, x, y, k): m.NewBoolVar("D%d_%d_%d_%d" % (r, x, y, k))
         for r in range(R + 1) for x in range(L) for y in range(L) for k in range(z)}
    A = {}

    for r in range(R):
        # --- theta: paridades de columna ---
        C = {(x, k): m.NewBoolVar("C%d_%d_%d" % (r, x, k))
             for x in range(L) for k in range(z)}
        for x in range(L):
            for k in range(z):
                m.AddBoolXOr([C[(x, k)].Not()] + [D[(r, x, y, k)] for y in range(L)])

        # --- theta: corrección ---
        Dth = {(x, k): m.NewBoolVar("Dth%d_%d_%d" % (r, x, k))
               for x in range(L) for k in range(z)}
        for x in range(L):
            for k in range(z):
                m.AddBoolXOr([Dth[(x, k)].Not(), C[((x - 1) % L, k)],
                              C[((x + 1) % L, (k - 1) % z)]])

        # --- theta: estado corregido ---
        Dt = {(x, y, k): m.NewBoolVar("Dt%d_%d_%d_%d" % (r, x, y, k))
              for x in range(L) for y in range(L) for k in range(z)}
        for x in range(L):
            for y in range(L):
                for k in range(z):
                    m.AddBoolXOr([Dt[(x, y, k)].Not(), D[(r, x, y, k)], Dth[(x, k)]])

        # --- rho + pi: alias de variables, sin crear ninguna nueva ---
        Drp = {}
        for x in range(L):
            for y in range(L):
                nx, ny, rot = y, (2 * x + 3 * y) % L, ROT[x][y] % z
                for k in range(z):
                    Drp[(nx, ny, (k + rot) % z)] = Dt[(x, y, k)]

        # --- chi: restricción de tabla (transiciones válidas de la DDT) ---
        for y in range(L):
            for k in range(z):
                entrada = [Drp[(i, y, k)] for i in range(L)]
                salida = [D[(r + 1, i, y, k)] for i in range(L)]
                m.AddAllowedAssignments(entrada + salida, TUPLAS)
                a = m.NewBoolVar("A%d_%d_%d" % (r, y, k))
                m.AddMaxEquality(a, entrada)   # actividad = OR de los bits de entrada
                A[(r, y, k)] = a

    if simetria:
        # La ronda sin iota conmuta con la traslación en z, luego puede exigirse
        # sin pérdida que el slice z=0 de la diferencia de entrada sea no nulo.
        # Implica además la no trivialidad.
        m.AddBoolOr([D[(0, x, y, 0)] for x in range(L) for y in range(L)])
    else:
        m.AddBoolOr([D[(0, x, y, k)] for x in range(L) for y in range(L)
                     for k in range(z)])

    m.Minimize(sum(A.values()))
    return m, A


def resolver(R, z, limite=120, simetria=True, hilos=8):
    """Resuelve un caso y devuelve cotas y estado REAL de terminación."""
    m, A = construir(R, z, simetria)
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = float(limite)
    s.parameters.num_search_workers = max(1, int(hilos))

    t0 = time.time()
    estado = s.Solve(m)
    segundos = time.time() - t0

    hay_sol = estado in (cp_model.OPTIMAL, cp_model.FEASIBLE)
    certificado = (estado == cp_model.OPTIMAL)
    ub = int(s.ObjectiveValue()) if hay_sol else None
    lb = int(s.BestObjectiveBound()) if hay_sol else R

    por_ronda = None
    if hay_sol:
        por_ronda = [sum(s.Value(A[(r, y, k)]) for y in range(L) for k in range(z))
                     for r in range(R)]

    return {"R": R, "z": z, "cajas_por_ronda": 5 * z,
            "estado": s.StatusName(estado), "certificado": certificado,
            "cota_inferior": max(lb, R), "cota_superior": ub,
            "por_ronda": por_ronda, "segundos": round(segundos, 1),
            "limite": limite, "hilos": hilos}


# ===========================================================================
# Consolidación por monotonía
# ===========================================================================
def consolidar(resultados):
    """Propaga las cotas inferiores hacia R mayores usando m(R) >= m(R-1).

    Añade a cada resultado:
      cota_inferior_consolidada : max{LB(r) : r <= R, mismo z}
      exacto                    : True si la cota consolidada iguala a la superior
      origen_cota               : R del que proviene la cota inferior empleada
    """
    por_z = {}
    for d in resultados:
        por_z.setdefault(d["z"], []).append(d)

    for z, lista in por_z.items():
        lista.sort(key=lambda d: d["R"])
        mejor, origen = 0, None
        for d in lista:
            if d["cota_inferior"] > mejor:
                mejor, origen = d["cota_inferior"], d["R"]
            d["cota_inferior_consolidada"] = mejor
            d["origen_cota"] = origen
            ub = d["cota_superior"]
            d["exacto"] = bool(ub is not None and ub == mejor)
    return resultados


def magnitudes(d):
    """Probabilidad y pares garantizados a partir de la cota INFERIOR."""
    n = d.get("cota_inferior_consolidada", d["cota_inferior"])
    peso = -DP_MAX_LOG2 * n           # 2n
    return {"prob_log2_garantizada": -peso, "pares_log2_garantizados": peso}


# ===========================================================================
# Ejecución (con paralelismo entre casos)
# ===========================================================================
def _tarea(args):
    R, z, limite, hilos = args
    return resolver(R, z, limite=limite, hilos=hilos)


def ejecutar(rondas, zetas, limite, workers, hilos=None):
    nucleos = os.cpu_count() or 4
    if hilos is None:
        hilos = max(1, nucleos // max(1, workers))
    casos = [(R, z, limite, hilos) for z in zetas for R in rondas]

    if workers <= 1:
        salida = [_tarea(c) for c in casos]
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            salida = list(ex.map(_tarea, casos))
    return consolidar(salida)


def tabla(resultados):
    cab = (" z |  R | mínimo        | por ronda            | estado    "
           "| pares (garant.) | tiempo")
    print(cab); print("-" * len(cab))
    for d in sorted(resultados, key=lambda a: (a["z"], a["R"])):
        mg = magnitudes(d)
        lb, ub = d["cota_inferior_consolidada"], d["cota_superior"]
        if d["exacto"]:
            minimo = "%d (exacto)" % ub
        elif ub is None:
            minimo = ">= %d" % lb
        else:
            minimo = "[%d, %d]" % (lb, ub)
        pr = ",".join(map(str, d["por_ronda"])) if d["por_ronda"] else "-"
        if len(pr) > 20:
            pr = pr[:17] + "..."
        print(" %d | %2d | %-13s | %-20s | %-9s | 2^%-13d | %5.1fs"
              % (d["z"], d["R"], minimo, pr, d["estado"],
                 mg["pares_log2_garantizados"], d["segundos"]))
    print("\nLa columna de pares se deriva de la cota INFERIOR consolidada: es la")
    print("que sustenta una garantía frente a un atacante. Las cotas inferiores se")
    print("propagan a R mayores porque m(R) es no decreciente en R.")


def rango(texto):
    """Interpreta '3', '1-10' o '1,3,5' como lista de rondas."""
    if "-" in texto:
        a, b = texto.split("-", 1)
        return list(range(int(a), int(b) + 1))
    return [int(v) for v in texto.split(",")]


def main():
    ap = argparse.ArgumentParser(
        description="Certificación CP-SAT del mínimo de cajas-S activas")
    ap.add_argument("--rondas", default="1-3",
                    help="rondas a analizar: '3', '1-10' o '1,3,5' (por omisión 1-3)")
    ap.add_argument("--z", default="4,8", help="tamaños de palabra (por omisión 4,8)")
    ap.add_argument("--tiempo", type=int, default=120,
                    help="límite de tiempo por caso, en segundos")
    ap.add_argument("--workers", type=int, default=1,
                    help="procesos en paralelo (casos simultáneos)")
    ap.add_argument("--hilos", type=int, default=None,
                    help="hilos internos de CP-SAT por caso")
    ap.add_argument("--salida", default=None, help="guardar resultados en JSON")
    a = ap.parse_args()

    rondas = rango(a.rondas)
    zetas = [int(v) for v in a.z.split(",")]

    print("CP-SAT | rondas %s | z %s | %d s por caso | %d proceso(s)\n"
          % (a.rondas, a.z, a.tiempo, a.workers))
    t0 = time.time()
    res = ejecutar(rondas, zetas, a.tiempo, a.workers, a.hilos)
    print()
    tabla(res)
    print("\nTiempo total: %.1f s" % (time.time() - t0))

    if a.salida:
        for d in res:
            d.update(magnitudes(d))
        with open(a.salida, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=2, ensure_ascii=False)
        print("Resultados guardados en %s" % a.salida)


if __name__ == "__main__":
    main()
