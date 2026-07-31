#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modelo CP-SAT de la variante ligera propuesta.

Estado: 5 carriles de w bits (5w bits), organizado en rodajas (bit-slicing).
Ronda:
  1) Capa lineal LAMBDA, solo XOR y rotaciones, dentro de cada carril:
         L_i <- L_i (+) rot(L_i, a_i) (+) rot(L_i, b_i)
  2) Capa no lineal CHI sobre columnas de 5 bits (una por indice z).
     En la version PARCIAL solo se aplica a un subconjunto de columnas; el resto
     pasa por la identidad. Esto reduce el numero de compuertas AND por ronda.

La difusion entre carriles proviene unicamente de chi; la difusion a lo largo
del carril, unicamente de las rotaciones. Es la estructura de Ascon.
"""
import sys, time
from ortools.sat.python import cp_model

L = 5
DP_MAX_LOG2 = -2

def sbox_chi(x):
    b = [(x >> i) & 1 for i in range(5)]
    o = [b[i] ^ ((~b[(i+1) % 5]) & b[(i+2) % 5]) & 1 for i in range(5)]
    return sum(o[i] << i for i in range(5))

DDT = {a: sorted({sbox_chi(x) ^ sbox_chi(x ^ a) for x in range(32)}) for a in range(32)}
TUPLAS = []
for a, bl in DDT.items():
    for b in bl:
        TUPLAS.append([(a >> i) & 1 for i in range(5)] + [(b >> i) & 1 for i in range(5)])


def columnas_activas(r, w, fraccion):
    """Que columnas llevan chi en la ronda r. fraccion=1 -> todas."""
    if fraccion == 1:
        return set(range(w))
    paso = int(round(1 / fraccion))
    return {z for z in range(w) if (z + r) % paso == 0}


def construir(R, w, rots, fraccion=1.0, simetria=True):
    m = cp_model.CpModel()
    D = {(r, i, z): m.NewBoolVar("D%d_%d_%d" % (r, i, z))
         for r in range(R + 1) for i in range(L) for z in range(w)}
    A = {}

    for r in range(R):
        # --- capa lineal: T[i][z] = D[i][z] (+) D[i][z-a] (+) D[i][z-b] ---
        T = {}
        for i in range(L):
            a, b = rots[i]
            for z in range(w):
                t = m.NewBoolVar("T%d_%d_%d" % (r, i, z))
                m.AddBoolXOr([t.Not(), D[(r, i, z)],
                              D[(r, i, (z - a) % w)], D[(r, i, (z - b) % w)]])
                T[(i, z)] = t

        # --- capa no lineal (posiblemente parcial) ---
        activas = columnas_activas(r, w, fraccion)
        for z in range(w):
            entrada = [T[(i, z)] for i in range(L)]
            salida = [D[(r + 1, i, z)] for i in range(L)]
            if z in activas:
                m.AddAllowedAssignments(entrada + salida, TUPLAS)
                a = m.NewBoolVar("A%d_%d" % (r, z))
                m.AddMaxEquality(a, entrada)
                A[(r, z)] = a
            else:
                for i in range(L):           # identidad
                    m.Add(salida[i] == entrada[i])

    if simetria:
        m.AddBoolOr([D[(0, i, 0)] for i in range(L)])
    else:
        m.AddBoolOr([D[(0, i, z)] for i in range(L) for z in range(w)])

    m.Minimize(sum(A.values()))
    return m, A


def resolver(R, w, rots, fraccion=1.0, limite=120, hilos=8):
    m, A = construir(R, w, rots, fraccion)
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = float(limite)
    s.parameters.num_search_workers = hilos
    t0 = time.time(); est = s.Solve(m); dt = time.time() - t0
    hay = est in (cp_model.OPTIMAL, cp_model.FEASIBLE)
    ub = int(s.ObjectiveValue()) if hay else None
    lb = int(s.BestObjectiveBound()) if hay else 0
    pr = None
    if hay:
        pr = [sum(s.Value(v) for (rr, z), v in A.items() if rr == r) for r in range(R)]
    return {"R": R, "w": w, "fraccion": fraccion, "estado": s.StatusName(est),
            "certificado": est == cp_model.OPTIMAL, "lb": max(lb, 0), "ub": ub,
            "por_ronda": pr, "segundos": round(dt, 1)}


if __name__ == "__main__":
    w = 40
    rots = [(19, 28), (21, 39), (1, 6), (10, 17), (7, 31)]
    print("Variante: 5 carriles x %d bits = %d bits\n" % (w, 5 * w))
    print(" chi        | R | minimo      | por ronda | estado   | tiempo")
    print("------------+---+-------------+-----------+----------+--------")
    for frac, nom in ((1.0, "completo"), (0.5, "parcial 1/2")):
        for R in (1, 2):
            d = resolver(R, w, rots, frac, limite=60)
            mm = "%d (exacto)" % d["ub"] if d["certificado"] else "[%s, %s]" % (d["lb"], d["ub"])
            print(" %-10s | %d | %-11s | %-9s | %-8s | %5.1fs"
                  % (nom, R, mm, ",".join(map(str, d["por_ronda"])) if d["por_ronda"] else "-",
                     d["estado"], d["segundos"]))
