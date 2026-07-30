#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modelo MILP mejorado para Keccak dinámico con:
- Reformulación Convex Hull de la DDT (sin Big‑M).
- Estrategia de decisión (preguntar por cota < cota_superior).
- Rompimiento de simetría rotacional.
- Parámetros avanzados de HiGHS.
"""

import math
import os
import tempfile
import time

import pulp

try:
    from sha3_nucleo import ROT, ddt_transiciones
except ImportError:
    from src.sha3_nucleo import ROT, ddt_transiciones

L = 5
DP_MAX_LOG2 = -2
DDT = ddt_transiciones()
DDT_LIST = [(a, b) for a in range(32) for b in DDT[a]]   # todos los pares válidos

# ===========================================================================
# 1. Capas lineales (igual que antes)
# ===========================================================================
def rotl(x, n, w):
    n %= w
    return ((x << n) | (x >> (w - n))) & ((1 << w) - 1) if n else x

def theta_z(A, w):
    C = [A[x][0] ^ A[x][1] ^ A[x][2] ^ A[x][3] ^ A[x][4] for x in range(L)]
    D = [C[(x - 1) % L] ^ rotl(C[(x + 1) % L], 1, w) for x in range(L)]
    return [[A[x][y] ^ D[x] for y in range(L)] for x in range(L)]

def rho_pi_z(A, w):
    R_ = [[rotl(A[x][y], ROT[x][y] % w, w) for y in range(L)] for x in range(L)]
    B = [[0] * L for _ in range(L)]
    for x in range(L):
        for y in range(L):
            B[y][(2 * x + 3 * y) % L] = R_[x][y]
    return B

def capa_lineal(A, w):
    return rho_pi_z(theta_z(A, w), w)

# ===========================================================================
# 2. Búsqueda de trayectorias (cota superior) - sin cambios
# ===========================================================================
DDT_MIN_PESO = {a: min([b for b in DDT[a] if b], key=lambda b: bin(b).count("1"))
                for a in range(1, 32)}

def propagar(CI, R, w):
    total, por_ronda = 0, []
    for r in range(R):
        activas = 0
        siguiente = [[0] * L for _ in range(L)]
        for y in range(L):
            for k in range(w):
                a = sum(((CI[x][y] >> k) & 1) << x for x in range(L))
                if a:
                    activas += 1
                    b = DDT_MIN_PESO[a]
                    for x in range(L):
                        siguiente[x][y] |= ((b >> x) & 1) << k
        por_ronda.append(activas)
        total += activas
        if r < R - 1:
            CI = capa_lineal(siguiente, w)
    return total, por_ronda

def cota_superior(R, w, pasadas=4):
    cero = lambda: [[0] * L for _ in range(L)]
    arranques = []
    for x in range(L):
        for y in range(L):
            for k in range(w):
                c = cero(); c[x][y] = 1 << k
                arranques.append(c)
    for y in range(L):
        for k in range(w):
            for x1 in range(L):
                for x2 in range(x1 + 1, L):
                    c = cero(); c[x1][y] = 1 << k; c[x2][y] = 1 << k
                    arranques.append(c)

    mejor, mejorC, mejorPR = None, None, None
    for c in arranques:
        t, pr = propagar(c, R, w)
        if mejor is None or t < mejor:
            mejor, mejorC, mejorPR = t, [f[:] for f in c], pr

    for _ in range(pasadas):
        mejoro = False
        for x in range(L):
            for y in range(L):
                for k in range(w):
                    cand = [f[:] for f in mejorC]
                    cand[x][y] ^= (1 << k)
                    if not any(any(f) for f in cand):
                        continue
                    t, pr = propagar(cand, R, w)
                    if t < mejor:
                        mejor, mejorC, mejorPR, mejoro = t, cand, pr, True
        if not mejoro:
            break
    return mejor, mejorPR

# ===========================================================================
# 3. Modelo MILP mejorado (Convex Hull + Simetría)
# ===========================================================================
_contador = [0]

def xor_lineal(prob, a, b, out):
    _contador[0] += 1
    t = pulp.LpVariable("xt_%d" % _contador[0], cat="Binary")
    prob += a + b - 2 * t == out

def construir_milp(R, z, limite_superior=None):
    """
    Construye el modelo MILP con reformulación Convex Hull para la DDT.
    Si se proporciona limite_superior, se añade la restricción de decisión:
        total_cajas < limite_superior
    """
    _contador[0] = 0
    prob = pulp.LpProblem("Keccak_R%d_z%d" % (R, z), pulp.LpMinimize)

    # Variables de diferencia
    D = pulp.LpVariable.dicts(
        "D", ((r, x, y, k) for r in range(R + 1) for x in range(L)
              for y in range(L) for k in range(z)), cat="Binary")
    # Variables de actividad
    A = pulp.LpVariable.dicts(
        "A", ((r, y, k) for r in range(R) for y in range(L) for k in range(z)),
        cat="Binary")

    # ROMBO DE SIMETRÍA: fijamos el bit (0,0,0) de la entrada a 1
    prob += D[(0, 0, 0, 0)] == 1

    # Preparamos la lista de pares de la DDT (global)
    pares = DDT_LIST

    for r in range(R):
        # --- theta ---
        C = pulp.LpVariable.dicts("C%d" % r, ((x, k) for x in range(L) for k in range(z)), cat="Binary")
        for x in range(L):
            for k in range(z):
                acc = D[(r, x, 0, k)]
                for y in range(1, L):
                    nxt = pulp.LpVariable("cc_%d_%d_%d_%d" % (r, x, k, y), cat="Binary")
                    xor_lineal(prob, acc, D[(r, x, y, k)], nxt)
                    acc = nxt
                prob += C[(x, k)] == acc

        Dth = pulp.LpVariable.dicts("Dth%d" % r, ((x, k) for x in range(L) for k in range(z)), cat="Binary")
        for x in range(L):
            for k in range(z):
                xor_lineal(prob, C[((x - 1) % L, k)], C[((x + 1) % L, (k - 1) % z)], Dth[(x, k)])

        Dt = pulp.LpVariable.dicts("Dt%d" % r, ((x, y, k) for x in range(L) for y in range(L) for k in range(z)), cat="Binary")
        for x in range(L):
            for y in range(L):
                for k in range(z):
                    xor_lineal(prob, D[(r, x, y, k)], Dth[(x, k)], Dt[(x, y, k)])

        # --- rho + pi (reindexación) ---
        Drp = pulp.LpVariable.dicts("Drp%d" % r, ((x, y, k) for x in range(L) for y in range(L) for k in range(z)), cat="Binary")
        for x in range(L):
            for y in range(L):
                nx, ny, rot = y, (2 * x + 3 * y) % L, ROT[x][y] % z
                for k in range(z):
                    prob += Drp[(nx, ny, (k + rot) % z)] == Dt[(x, y, k)]

        # --- chi: reformulación Convex Hull ---
        for y in range(L):
            for k in range(z):
                # Variables lambda (continuas, no negativas, suma = 1)
                lambdas = pulp.LpVariable.dicts(
                    f"lam_{r}_{y}_{k}",
                    range(len(pares)),
                    lowBound=0,
                    cat="Continuous"
                )
                prob += pulp.lpSum(lambdas) == 1

                # vin y vout como combinación convexa de los pares
                vin_expr = pulp.lpSum(lambdas[i] * a for i, (a, _) in enumerate(pares))
                vout_expr = pulp.lpSum(lambdas[i] * b for i, (_, b) in enumerate(pares))
                # Actividad: 1 si a != 0
                act_expr = pulp.lpSum(lambdas[i] * (1 if a != 0 else 0) for i, (a, _) in enumerate(pares))

                # Conectar con las variables D de entrada y salida
                # vin = Σ D[r,i,y,k] * 2^i
                prob += vin_expr == pulp.lpSum(D[(r, i, y, k)] * (1 << i) for i in range(L))
                prob += vout_expr == pulp.lpSum(D[(r + 1, i, y, k)] * (1 << i) for i in range(L))
                # Actividad
                prob += A[(r, y, k)] == act_expr

    # Objetivo: minimizar total de cajas activas
    total_cajas = pulp.lpSum(A[(r, y, k)] for r in range(R) for y in range(L) for k in range(z))
    prob += total_cajas

    # Si se da una cota superior, añadimos restricción de decisión
    if limite_superior is not None:
        prob += total_cajas <= limite_superior - 1

    return prob, A

# ===========================================================================
# 4. Resolución con HiGHS (con parámetros mejorados)
# ===========================================================================
def resolver_milp(R, z, limite_tiempo=60, cota_superior_ub=None):
    """
    Resuelve el MILP con HiGHS.
    Si cota_superior_ub no es None, se usa la estrategia de decisión.
    Devuelve (certificado, incumbente, cota_dual, estado_str)
    """
    import highspy

    prob, _ = construir_milp(R, z, limite_superior=cota_superior_ub)
    ruta = tempfile.mktemp(suffix=".lp")
    prob.writeLP(ruta)

    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("time_limit", float(limite_tiempo))
    h.setOptionValue("mip_rel_gap", 0.0)
    # Parámetros adicionales para fortalecer la búsqueda
    h.setOptionValue("presolve", "on")
    h.setOptionValue("mip_max_nodes", 1e9)
    h.setOptionValue("parallel", "on")
    h.setOptionValue("mip_detect_symmetry", "on")
    h.setOptionValue("mip_pool_cleanup", "on")  # limpia soluciones dominadas

    h.readModel(ruta)
    h.run()

    estado = h.getModelStatus()
    certificado = (estado == highspy.HighsModelStatus.kOptimal)
    obj = h.getObjectiveValue()
    dual = h.getInfo().mip_dual_bound
    os.remove(ruta)

    incumbente = int(round(obj)) if obj not in (float("inf"), float("-inf")) else None
    return certificado, incumbente, dual, str(estado).split(".")[-1]

# ===========================================================================
# 5. Análisis con estrategia de decisión
# ===========================================================================
def analizar_con_decision(R, z, limite_tiempo=60):
    """
    Combina búsqueda de trayectorias (cota superior) y MILP con decisión
    para intentar certificar el óptimo.
    """
    # 1. Obtener cota superior mediante búsqueda heurística
    ub_tray, por_ronda = cota_superior(R, z)
    print(f"[*] Cota superior inicial para R={R}, z={z}: {ub_tray} (por ronda: {por_ronda})")

    # 2. Resolver MILP con restricción de decisión: total <= ub_tray - 1
    t0 = time.time()
    certificado, incumbente, dual, nombre = resolver_milp(
        R, z, limite_tiempo, cota_superior_ub=ub_tray
    )
    segundos = time.time() - t0

    # 3. Interpretar resultado
    if certificado and incumbente is not None:
        # Si el problema de decisión es factible, hemos encontrado una trayectoria mejor
        # que la heurística. Actualizamos la cota superior.
        ub = incumbente
        # Pero como era <= ub_tray -1, entonces ub < ub_tray
        print(f"[*] ¡Encontrada trayectoria mejor! Nueva cota superior: {ub}")
    else:
        # Si es infactible o tiempo límite, mantenemos la cota superior
        ub = ub_tray

    # 4. Cota inferior: máximo entre R y la cota dual (si está definida)
    if dual in (float("inf"), float("-inf")) or dual != dual:  # nan
        lb = R
    else:
        lb = max(R, int(math.ceil(dual - 1e-6)))

    # 5. Si el MILP encontró una solución mejor, actualizamos por_ronda (opcional)
    # (No tenemos la distribución por ronda de la solución MILP, la dejamos de la heurística)

    peso = -DP_MAX_LOG2 * ub
    return {
        "R": R, "z": z, "cajas_por_ronda": 5 * z,
        "cota_inferior": lb, "cota_superior": ub,
        "por_ronda": por_ronda,
        "certificado": certificado and (incumbente == ub),
        "incumbente_milp": incumbente,
        "cota_dual": dual,
        "estado_solver": nombre,
        "segundos": segundos,
        "prob_log2": -peso,
        "pares_log2": peso,
    }

# (Opcional) Mantener función analizar antigua por compatibilidad
def analizar(R, z, limite_tiempo=60):
    return analizar_con_decision(R, z, limite_tiempo)


if __name__ == "__main__":
    print("Prueba rápida: R=1, z=4")
    d = analizar_con_decision(1, 4, limite_tiempo=10)
    print(d)