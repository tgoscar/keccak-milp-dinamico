#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Certificación del número mínimo de cajas-S activas mediante CP-SAT (OR-Tools).

Misma formulación diferencial que el modelo MILP de `milp_keccak.py`, pero con
las restricciones que CP-SAT trata de forma nativa:

  * XOR nativo (`AddBoolXOr`) en lugar de variables auxiliares con a + b - 2t = c.
  * chi mediante restricción de TABLA (`AddAllowedAssignments`) con las 317
    transiciones válidas de la DDT: sin big-M y sin variables de selección.

El modelo MILP dedica cerca del 90 % de sus variables al encoding big-M de la
DDT (12 680 de 14 140 para R=2, z=4), cuya relajación lineal es muy débil. CP-SAT
no necesita ese encoding y, además, aprende cláusulas sobre la estructura XOR;
en la práctica certifica en segundos lo que el MILP no cierra en minutos.

Ruptura de simetría
-------------------
La ronda sin iota conmuta con la traslación a lo largo de z: theta rota una
posición en z, rho traslada cada carril una constante, pi permuta carriles y chi
actúa dentro del slice. (iota existe precisamente para romper esa invariancia.)
Por tanto toda trayectoria pertenece a una órbita de tamaño divisor de z, y se
puede exigir sin pérdida que el slice z=0 de la diferencia de entrada sea no
nulo. El ahorro está acotado por z, de modo que es útil pero modesto.

Uso:
    python src/cpsat_keccak.py            # los 6 casos, 120 s cada uno
    python src/cpsat_keccak.py 300        # 300 s por caso
"""

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ortools.sat.python import cp_model

try:
    from sha3_nucleo import ROT, ddt_transiciones
except ImportError:
    from src.sha3_nucleo import ROT, ddt_transiciones

L = 5
DP_MAX_LOG2 = -2

# Transiciones válidas de la DDT como tuplas de 10 bits (5 de entrada, 5 de salida)
DDT = ddt_transiciones()
TUPLAS = []
for _a, _bl in DDT.items():
    for _b in _bl:
        TUPLAS.append([(_a >> i) & 1 for i in range(5)] +
                      [(_b >> i) & 1 for i in range(5)])


def construir(R, z, simetria=True):
    """Construye el modelo CP-SAT. Devuelve (modelo, dict de actividades)."""
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
                # C = XOR de los 5 bits  <=>  XOR(not C, bits...) impar
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

        # --- chi: restricción de tabla + actividad ---
        for y in range(L):
            for k in range(z):
                entrada = [Drp[(i, y, k)] for i in range(L)]
                salida = [D[(r + 1, i, y, k)] for i in range(L)]
                m.AddAllowedAssignments(entrada + salida, TUPLAS)
                a = m.NewBoolVar("A%d_%d_%d" % (r, y, k))
                m.AddMaxEquality(a, entrada)      # a = OR de los bits de entrada
                A[(r, y, k)] = a

    if simetria:
        m.AddBoolOr([D[(0, x, y, 0)] for x in range(L) for y in range(L)])
    else:
        m.AddBoolOr([D[(0, x, y, k)] for x in range(L) for y in range(L)
                     for k in range(z)])

    m.Minimize(sum(A.values()))
    return m, A


def resolver(R, z, limite=120, simetria=True, hilos=8):
    """Resuelve y devuelve un dict con el mínimo, las cotas y el estado.

    `certificado` es True sólo si CP-SAT devuelve OPTIMAL. En caso contrario se
    informa el intervalo [cota_inferior, cota_superior] realmente demostrado.
    """
    m, A = construir(R, z, simetria)
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = float(limite)
    s.parameters.num_search_workers = hilos

    t0 = time.time()
    estado = s.Solve(m)
    segundos = time.time() - t0

    nombre = s.StatusName(estado)
    hay_sol = estado in (cp_model.OPTIMAL, cp_model.FEASIBLE)
    ub = int(s.ObjectiveValue()) if hay_sol else None
    lb = int(s.BestObjectiveBound()) if hay_sol else R
    certificado = (estado == cp_model.OPTIMAL)

    por_ronda = None
    if hay_sol:
        por_ronda = [sum(s.Value(A[(r, y, k)]) for y in range(L) for k in range(z))
                     for r in range(R)]

    # El peso diferencial se calcula con la cota INFERIOR del número de cajas:
    # es la que sustenta una garantía de seguridad.
    peso_garantizado = -DP_MAX_LOG2 * max(lb, R)
    return {
        "R": R, "z": z, "estado": nombre, "certificado": certificado,
        "cota_inferior": max(lb, R), "cota_superior": ub, "por_ronda": por_ronda,
        "segundos": segundos, "pares_log2_garantizado": peso_garantizado,
    }


if __name__ == "__main__":
    limite = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    print("Certificación con CP-SAT — límite %d s por caso\n" % limite)
    cab = " z | R | mínimo    | por ronda | estado    | pares    | tiempo"
    print(cab)
    print("-" * len(cab))
    for z in (4, 8):
        for R in (1, 2, 3):
            d = resolver(R, z, limite=limite)
            minimo = ("%d (exacto)" % d["cota_superior"] if d["certificado"]
                      else "[%s, %s]" % (d["cota_inferior"], d["cota_superior"]))
            print(" %d | %d | %-9s | %-9s | %-9s | 2^%-6d | %5.1fs"
                  % (z, R, minimo,
                     ",".join(map(str, d["por_ronda"])) if d["por_ronda"] else "-",
                     d["estado"], d["pares_log2_garantizado"], d["segundos"]))
    print("\nLos pares se derivan de la cota INFERIOR de cajas activas (2^2n), que es")
    print("la que sustenta una garantía de seguridad frente a un atacante.")
