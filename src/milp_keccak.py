#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modelo MILP para el número mínimo de cajas-S activas en Keccak con rondas
dinámicas, y búsqueda dirigida de trayectorias diferenciales válidas.

Contenido:
  * `cota_superior(R, z)`  -> mejor trayectoria diferencial válida hallada
                              (cota SUPERIOR del mínimo).
  * `construir_milp(R, z)` -> modelo MILP completo (PuLP).
  * `resolver_milp(R, z)`  -> resuelve con HiGHS y devuelve el estado REAL de
                              terminación, el incumbente y la cota dual.
  * `analizar(R, z)`       -> combina ambos y devuelve el intervalo
                              [cota inferior, cota superior].

Sobre el estado del solver
--------------------------
La capa de modelado (PuLP) puede reportar "Optimal" aunque HiGHS solo haya
alcanzado el límite de tiempo con gap del 100 %. Aquí se consulta highspy
directamente (kOptimal frente a kTimeLimit) para no informar optimalidad
donde no la hay.

Sobre las cotas
---------------
  * Cota inferior: R. La capa lineal es biyectiva sobre GF(2) y chi es
    biyectiva, luego la ronda es una biyección y una diferencia no nula no
    puede anularse: cada ronda aporta al menos una caja-S activa.
  * Cota superior: para R = 1 el MILP certifica el óptimo. Para R >= 2 la
    relajación lineal es demasiado débil (la cota dual se queda en 0) y la
    optimalidad no se certifica; se reporta la mejor trayectoria válida.
"""

import math
import os
import tempfile
import time

import pulp

try:
    from sha3_nucleo import ROT, ddt_transiciones
except ImportError:  # ejecución desde la raíz del repositorio
    from src.sha3_nucleo import ROT, ddt_transiciones

L = 5                       # lado del estado (5 x 5 carriles)
DP_MAX_LOG2 = -2            # probabilidad diferencial máxima de chi: 2^-2

DDT = ddt_transiciones()
PARES = [(a, b) for a, bl in DDT.items() for b in bl]
# Para cada diferencia de entrada, la salida válida de menor peso de Hamming.
DDT_MIN_PESO = {a: min([b for b in DDT[a] if b], key=lambda b: bin(b).count("1"))
                for a in range(1, 32)}


# ===========================================================================
# 1. Capas lineales aplicadas directamente sobre la diferencia
# ===========================================================================
def rotl(x, n, w):
    """Rotación circular a la izquierda sobre w bits."""
    n %= w
    return ((x << n) | (x >> (w - n))) & ((1 << w) - 1) if n else x


def theta_z(A, w):
    """theta sobre carriles de w bits (A[x][y] es un entero de w bits)."""
    C = [A[x][0] ^ A[x][1] ^ A[x][2] ^ A[x][3] ^ A[x][4] for x in range(L)]
    D = [C[(x - 1) % L] ^ rotl(C[(x + 1) % L], 1, w) for x in range(L)]
    return [[A[x][y] ^ D[x] for y in range(L)] for x in range(L)]


def rho_pi_z(A, w):
    """rho (offsets reducidos mod w) seguido de pi."""
    R_ = [[rotl(A[x][y], ROT[x][y] % w, w) for y in range(L)] for x in range(L)]
    B = [[0] * L for _ in range(L)]
    for x in range(L):
        for y in range(L):
            B[y][(2 * x + 3 * y) % L] = R_[x][y]
    return B


def capa_lineal(A, w):
    """lambda = pi o rho o theta."""
    return rho_pi_z(theta_z(A, w), w)


# ===========================================================================
# 2. Búsqueda dirigida de trayectorias diferenciales válidas (cota superior)
# ===========================================================================
def propagar(CI, R, w):
    """Propaga una diferencia post-lineal eligiendo salidas válidas de chi.

    En cada caja-S activa se toma la diferencia de salida válida de menor peso
    de Hamming (según la DDT) y se aplica la capa lineal para la ronda
    siguiente. Devuelve (total de cajas activas, lista por ronda).
    """
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
    """Mejor trayectoria válida hallada -> cota superior del mínimo.

    Arranques: (a) un solo bit activo; (b) dos bits en la misma fila de chi,
    que reproducen las trayectorias de bajo peso características de Keccak.
    Después, mejora local volteando un bit a la vez.
    """
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
# 3. Modelo MILP
# ===========================================================================
_contador = [0]


def xor_lineal(prob, a, b, out):
    """Impone out = a XOR b de forma exacta: a + b - 2t = out, con t binaria.

    Verificación: (0,0)->0, (0,1)->1, (1,0)->1, (1,1)->0. Véase
    tests/test_verificaciones.py.
    """
    _contador[0] += 1
    t = pulp.LpVariable("xt_%d" % _contador[0], cat="Binary")
    prob += a + b - 2 * t == out


def construir_milp(R, z):
    """Construye el modelo MILP para R rondas y tamaño de palabra z."""
    _contador[0] = 0
    prob = pulp.LpProblem("Keccak_R%d_z%d" % (R, z), pulp.LpMinimize)

    # Diferencias del estado: R+1 capas
    D = pulp.LpVariable.dicts(
        "D", ((r, x, y, k) for r in range(R + 1) for x in range(L)
              for y in range(L) for k in range(z)), cat="Binary")
    # Actividad de cada caja-S
    A = pulp.LpVariable.dicts(
        "A", ((r, y, k) for r in range(R) for y in range(L) for k in range(z)),
        cat="Binary")

    for r in range(R):
        # --- theta: paridades de columna ---
        C = pulp.LpVariable.dicts("C%d" % r, ((x, k) for x in range(L)
                                              for k in range(z)), cat="Binary")
        for x in range(L):
            for k in range(z):
                acc = D[(r, x, 0, k)]
                for y in range(1, L):
                    nxt = pulp.LpVariable("cc_%d_%d_%d_%d" % (r, x, k, y), cat="Binary")
                    xor_lineal(prob, acc, D[(r, x, y, k)], nxt)
                    acc = nxt
                prob += C[(x, k)] == acc

        # --- theta: corrección D_theta[x][k] = C[x-1][k] XOR ROTL(C[x+1])[k] ---
        Dth = pulp.LpVariable.dicts("Dth%d" % r, ((x, k) for x in range(L)
                                                  for k in range(z)), cat="Binary")
        for x in range(L):
            for k in range(z):
                xor_lineal(prob, C[((x - 1) % L, k)], C[((x + 1) % L, (k - 1) % z)],
                           Dth[(x, k)])

        # --- theta: estado corregido ---
        Dt = pulp.LpVariable.dicts("Dt%d" % r, ((x, y, k) for x in range(L)
                                                for y in range(L) for k in range(z)),
                                   cat="Binary")
        for x in range(L):
            for y in range(L):
                for k in range(z):
                    xor_lineal(prob, D[(r, x, y, k)], Dth[(x, k)], Dt[(x, y, k)])

        # --- rho + pi: permutación pura, por reindexación ---
        Drp = pulp.LpVariable.dicts("Drp%d" % r, ((x, y, k) for x in range(L)
                                                  for y in range(L) for k in range(z)),
                                    cat="Binary")
        for x in range(L):
            for y in range(L):
                nx, ny, rot = y, (2 * x + 3 * y) % L, ROT[x][y] % z
                for k in range(z):
                    prob += Drp[(nx, ny, (k + rot) % z)] == Dt[(x, y, k)]

        # --- chi: transición válida de la DDT por caja-S ---
        for y in range(L):
            for k in range(z):
                vin = pulp.LpVariable("vin_%d_%d_%d" % (r, y, k), 0, 31, cat="Integer")
                prob += vin == pulp.lpSum(Drp[(i, y, k)] * (1 << i) for i in range(L))
                vout = pulp.LpVariable("vout_%d_%d_%d" % (r, y, k), 0, 31, cat="Integer")
                prob += vout == pulp.lpSum(D[(r + 1, i, y, k)] * (1 << i) for i in range(L))
                sel = []
                for idx, (a, b) in enumerate(PARES):
                    s = pulp.LpVariable("sel_%d_%d_%d_%d" % (r, y, k, idx), cat="Binary")
                    sel.append(s)
                    prob += vin - a <= (1 - s) * 31
                    prob += a - vin <= (1 - s) * 31
                    prob += vout - b <= (1 - s) * 31
                    prob += b - vout <= (1 - s) * 31
                prob += pulp.lpSum(sel) == 1
                # A = 1 si y solo si vin != 0
                prob += vin <= 31 * A[(r, y, k)]
                prob += A[(r, y, k)] <= vin

    # Objetivo: minimizar el total de cajas-S activas
    prob += pulp.lpSum(A[(r, y, k)] for r in range(R) for y in range(L) for k in range(z))
    # No trivialidad: al menos un bit de diferencia en la entrada.
    # (Fijar un bit concreto sesgaría el mínimo; véase docs/CORRECCIONES.md.)
    prob += pulp.lpSum(D[(0, x, y, k)] for x in range(L) for y in range(L)
                       for k in range(z)) >= 1
    return prob, A


def resolver_milp(R, z, limite_tiempo=60):
    """Resuelve el MILP con HiGHS y devuelve el estado REAL de terminación.

    Returns
    -------
    (certificado, incumbente, cota_dual, nombre_estado)
        certificado : True solo si HiGHS reporta kOptimal (gap cerrado).
        incumbente  : mejor solución entera hallada, o None si no halló ninguna.
    """
    try:
        import highspy
    except ImportError as e:
        # highspy y ortools empaquetan builds distintos de HiGHS y sus simbolos
        # colisionan: no pueden convivir en un mismo proceso. Si CP-SAT ya se
        # cargo, se degrada a CBC informando que el estado no es fiable.
        print("  [aviso] highspy no disponible en este proceso (%s)." % type(e).__name__)
        print("          Se usa CBC; su estado de terminacion no es fiable.")
        prob, _ = construir_milp(R, z)
        prob.solve(pulp.PULP_CBC_CMD(msg=0, timeLimit=limite_tiempo))
        obj = pulp.value(prob.objective)
        inc = int(round(obj)) if obj is not None else None
        # CBC puede devolver valores por debajo de la cota estructural R, que son
        # imposibles (la ronda es biyectiva). Se descartan en lugar de propagarlos.
        if inc is not None and inc < R:
            print("          CBC devolvio %d < R=%d (imposible); se descarta." % (inc, R))
            inc = None
        return False, inc, float("nan"), "CBC(sin certificar)"

    prob, _ = construir_milp(R, z)
    ruta = tempfile.mktemp(suffix=".lp")
    prob.writeLP(ruta)

    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("time_limit", float(limite_tiempo))
    h.setOptionValue("mip_rel_gap", 0.0)
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
# 4. Análisis combinado
# ===========================================================================
def analizar(R, z, limite_tiempo=60):
    """Devuelve un dict con cotas, estado y magnitudes de seguridad."""
    ub_tray, por_ronda = cota_superior(R, z)
    t0 = time.time()
    certificado, incumbente, dual, nombre = resolver_milp(R, z, limite_tiempo)
    segundos = time.time() - t0

    ub = ub_tray if incumbente is None else min(ub_tray, incumbente)
    if dual in (float("inf"), float("-inf")) or dual != dual:  # nan
        lb = R
    else:
        lb = max(R, int(math.ceil(dual - 1e-6)))

    peso = -DP_MAX_LOG2 * ub          # exponente: 2n
    return {
        "R": R, "z": z, "cajas_por_ronda": 5 * z,
        "cota_inferior": lb, "cota_superior": ub, "por_ronda": por_ronda,
        "certificado": certificado, "incumbente_milp": incumbente,
        "cota_dual": dual, "estado_solver": nombre, "segundos": segundos,
        "prob_log2": -peso, "pares_log2": peso,
    }


if __name__ == "__main__":
    print("Cajas-S activas en Keccak con rondas dinámicas")
    print("(probabilidad de trayectoria <= 2^-2n ; pares necesarios ~ 2^2n)\n")
    cab = " z | R | cota inf | cota sup | por ronda | estado del solver  | mínimo"
    print(cab)
    print("-" * len(cab))
    for z in (4, 8):
        for R in (1, 2, 3):
            d = analizar(R, z, limite_tiempo=60)
            estado = "óptimo certificado" if d["certificado"] else "límite de tiempo"
            minimo = ("%d (exacto)" % d["cota_superior"] if d["certificado"]
                      else "[%d, %d]" % (d["cota_inferior"], d["cota_superior"]))
            print(" %d | %d |   %3d    |   %3d    | %-9s | %-18s | %s"
                  % (z, R, d["cota_inferior"], d["cota_superior"],
                     ",".join(map(str, d["por_ronda"])), estado, minimo))


# ===========================================================================
# 5. Variante: envolvente convexa (Convex Hull) de la DDT, sin big-M
# ===========================================================================
# Motivación: el encoding big-M de chi consume cerca del 90 % de las variables
# del modelo y relaja mal. Una alternativa es describir el conjunto de
# transiciones válidas como la envolvente convexa de sus 317 puntos en
# {0,1}^10 (5 bits de entrada, 5 de salida), mediante multiplicadores lambda.
#
# CORRECCIÓN IMPORTANTE. La envolvente debe imponerse COMPONENTE A COMPONENTE:
# una ecuación por cada uno de los 10 bits. Colapsar esas 10 ecuaciones en dos
# ecuaciones escalares con potencias de dos,
#
#     sum_j lambda_j * a_j  ==  sum_i D_in[i] * 2^i        <-- INCORRECTO
#
# invalida el modelo: una combinación convexa de los enteros a_j alcanza valores
# que no corresponden a ninguna transición válida, de modo que la restricción de
# la DDT queda vacía. En pruebas, esa versión escalar admitía trayectorias con
# transiciones imposibles (por ejemplo, diferencia de entrada no nula con
# diferencia de salida nula, imposible por ser chi biyectiva) y devolvía valores
# por debajo del óptimo certificado. Véase docs/CORRECCIONES.md, defecto 6.
#
# La formulación por componentes sí es exacta: para un conjunto S de puntos 0/1
# se cumple conv(S) ∩ {0,1}^n = S, porque todo punto 0/1 de la envolvente es
# vértice del cubo unidad y por tanto pertenece a S.
#
# Resultado medido: es válida, pero NO competitiva. Para R=2, z=4 devuelve una
# cota de 7 sin certificar en 200 s, mientras CP-SAT certifica el óptimo (4) en
# segundos. Se conserva por completitud metodológica.

def construir_milp_hull(R, z, cota_decision=None):
    """Modelo MILP con la DDT descrita por su envolvente convexa (por componentes).

    cota_decision: si se indica K, añade sum(A) <= K y convierte el problema en
    uno de decisión (útil para intentar refutar la existencia de trayectorias
    con K o menos cajas activas).
    """
    _contador[0] = 0
    prob = pulp.LpProblem("Keccak_hull_R%d_z%d" % (R, z), pulp.LpMinimize)

    D = pulp.LpVariable.dicts(
        "D", ((r, x, y, k) for r in range(R + 1) for x in range(L)
              for y in range(L) for k in range(z)), cat="Binary")
    A = pulp.LpVariable.dicts(
        "A", ((r, y, k) for r in range(R) for y in range(L) for k in range(z)),
        cat="Binary")

    for r in range(R):
        C = pulp.LpVariable.dicts("C%d" % r, ((x, k) for x in range(L)
                                              for k in range(z)), cat="Binary")
        for x in range(L):
            for k in range(z):
                acc = D[(r, x, 0, k)]
                for y in range(1, L):
                    nxt = pulp.LpVariable("hc_%d_%d_%d_%d" % (r, x, k, y), cat="Binary")
                    xor_lineal(prob, acc, D[(r, x, y, k)], nxt)
                    acc = nxt
                prob += C[(x, k)] == acc

        Dth = pulp.LpVariable.dicts("hDth%d" % r, ((x, k) for x in range(L)
                                                   for k in range(z)), cat="Binary")
        for x in range(L):
            for k in range(z):
                xor_lineal(prob, C[((x - 1) % L, k)],
                           C[((x + 1) % L, (k - 1) % z)], Dth[(x, k)])

        Dt = pulp.LpVariable.dicts("hDt%d" % r, ((x, y, k) for x in range(L)
                                                 for y in range(L) for k in range(z)),
                                   cat="Binary")
        for x in range(L):
            for y in range(L):
                for k in range(z):
                    xor_lineal(prob, D[(r, x, y, k)], Dth[(x, k)], Dt[(x, y, k)])

        Drp = {}
        for x in range(L):
            for y in range(L):
                nx, ny, rot = y, (2 * x + 3 * y) % L, ROT[x][y] % z
                for k in range(z):
                    Drp[(nx, ny, (k + rot) % z)] = Dt[(x, y, k)]

        # --- chi mediante envolvente convexa, una ecuación por bit ---
        for y in range(L):
            for k in range(z):
                lam = pulp.LpVariable.dicts("lam_%d_%d_%d" % (r, y, k),
                                            range(len(PARES)), lowBound=0,
                                            cat="Continuous")
                prob += pulp.lpSum(lam.values()) == 1
                for i in range(L):
                    prob += pulp.lpSum(lam[j] * ((PARES[j][0] >> i) & 1)
                                       for j in range(len(PARES))) == Drp[(i, y, k)]
                    prob += pulp.lpSum(lam[j] * ((PARES[j][1] >> i) & 1)
                                       for j in range(len(PARES))) == D[(r + 1, i, y, k)]
                # Actividad derivada de los bits REALES de entrada, no de los lambda
                for i in range(L):
                    prob += A[(r, y, k)] >= Drp[(i, y, k)]
                prob += A[(r, y, k)] <= pulp.lpSum(Drp[(i, y, k)] for i in range(L))

    total = pulp.lpSum(A[(r, y, k)] for r in range(R) for y in range(L)
                       for k in range(z))
    prob += total
    prob += pulp.lpSum(D[(0, x, y, k)] for x in range(L) for y in range(L)
                       for k in range(z)) >= 1
    if cota_decision is not None:
        prob += total <= cota_decision
    return prob, A, D


def trayectoria_es_valida(D_valores, R, z, sup=None):
    """Comprueba que una trayectoria respeta la DDT en todas las cajas activas.

    D_valores: función (r, x, y, k) -> 0/1. Devuelve (numero_de_violaciones,
    lista de las primeras violaciones como (r, y, k, a, b)).
    """
    violaciones = []
    for r in range(R):
        estado = [[sum(D_valores(r, x, y, k) << k for k in range(z))
                   for y in range(L)] for x in range(L)]
        CI = capa_lineal(estado, z)
        for y in range(L):
            for k in range(z):
                a = sum(((CI[x][y] >> k) & 1) << x for x in range(L))
                b = sum(D_valores(r + 1, x, y, k) << x for x in range(L))
                if (a or b) and b not in DDT[a]:
                    violaciones.append((r, y, k, a, b))
    return len(violaciones), violaciones[:5]
