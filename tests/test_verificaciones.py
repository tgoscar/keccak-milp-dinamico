#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verificaciones de las propiedades afirmadas en el informe.

Ejecutar:  python tests/test_verificaciones.py     (o bien: pytest -q)

Cada función comprueba una afirmación concreta del informe, de modo que los
resultados publicados sean auditables sin necesidad de rehacer los experimentos.
"""

import hashlib
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from sha3_nucleo import (ROT, sbox_chi, ddt_completa, ddt_transiciones,
                         dp_maxima, sha3_256, intentos_a_rondas)

L = 5


# ---------------------------------------------------------------------------
def test_sha3_contra_hashlib():
    """El núcleo por pasos separados reproduce SHA-3-256 (FIPS 202)."""
    vectores = [b"", b"abc", b"a" * 135, b"a" * 136, b"a" * 137,
                b"Plataforma del Estado - modulo de autenticacion"]
    for m in vectores:
        assert sha3_256(m).hex() == hashlib.sha3_256(m).hexdigest()
    print("  [ok] SHA-3-256 coincide con hashlib en %d vectores" % len(vectores))


# ---------------------------------------------------------------------------
def test_chi_es_biyectiva():
    """chi es una permutación de 5 bits."""
    imagenes = {sbox_chi(x) for x in range(32)}
    assert len(imagenes) == 32
    print("  [ok] chi es biyectiva (32 imágenes distintas)")


def test_chi_sin_negacion_no_es_biyectiva():
    """Omitir el NOT de chi destruye la biyectividad.

    Contraejemplo: la diferencia (1,1,1,1,1) colapsa a (0,0,0,0,0), lo que
    permitiría que una diferencia no nula 'desaparezca'.
    """
    def chi_sin_not(x):
        b = [(x >> i) & 1 for i in range(5)]
        o = [b[i] ^ (b[(i + 1) % 5] & b[(i + 2) % 5]) for i in range(5)]
        return sum(o[i] << i for i in range(5))

    assert chi_sin_not(0b11111) == 0
    assert len({chi_sin_not(x) for x in range(32)}) < 32
    print("  [ok] sin el NOT, chi no es biyectiva (11111 -> 00000)")


# ---------------------------------------------------------------------------
def test_ddt_transiciones_y_dp_maxima():
    """La DDT tiene 317 transiciones válidas, ninguna fila vacía, DP_max = 2^-2."""
    trans = ddt_transiciones()
    assert sum(1 for a in range(32) if not trans[a]) == 0, "hay filas vacías"
    assert sum(len(v) for v in trans.values()) == 317
    assert trans[1] == [1, 9, 17, 25]
    mult, expo = dp_maxima()
    assert (mult, expo) == (8, -2)
    print("  [ok] DDT: 317 transiciones, 0 filas vacías, DP_max = 8/32 = 2^-2")


def test_ddt_filtro_espurio_vacia_filas():
    """Reproduce el defecto corregido: el filtro espurio vacía 21 de 32 filas."""
    def ddt_con_bug():
        d = {}
        for a in range(32):
            posibles = []
            for b in range(32):
                for x in range(32):
                    if (x ^ sbox_chi(x)) == a:              # filtro ajeno a la DDT
                        if (sbox_chi(x) ^ sbox_chi(x ^ a)) == b:
                            posibles.append(b)
                            break
            d[a] = posibles
        return d

    con_bug = ddt_con_bug()
    correcta = ddt_transiciones()
    vacias = [a for a in range(32) if not con_bug[a]]
    omitidas = sum(len(set(correcta[a]) - set(con_bug[a])) for a in range(32))
    assert len(vacias) == 21
    assert omitidas == 301
    print("  [ok] el filtro espurio deja 21 filas vacías y omite 301 transiciones")


def test_ddt_filas_planas():
    """Cada fila no nula de la DDT tiene multiplicidad uniforme."""
    ddt = ddt_completa()
    for a in range(1, 32):
        valores = {ddt[a][b] for b in range(32) if ddt[a][b] > 0}
        assert len(valores) == 1, "fila %d con multiplicidades distintas" % a
    print("  [ok] todas las filas de la DDT tienen multiplicidad uniforme")


# ---------------------------------------------------------------------------
def test_gadget_xor_exacto():
    """a + b - 2t = c reproduce el XOR en los cuatro casos."""
    import pulp
    for va in (0, 1):
        for vb in (0, 1):
            prob = pulp.LpProblem("t", pulp.LpMinimize)
            a = pulp.LpVariable("a", cat="Binary")
            b = pulp.LpVariable("b", cat="Binary")
            c = pulp.LpVariable("c", cat="Binary")
            t = pulp.LpVariable("t", cat="Binary")
            prob += a == va
            prob += b == vb
            prob += a + b - 2 * t == c
            prob += 0
            prob.solve(pulp.PULP_CBC_CMD(msg=0))
            assert pulp.LpStatus[prob.status] == "Optimal"
            assert int(round(pulp.value(c))) == (va ^ vb)
    print("  [ok] el gadget XOR exacto resuelve los 4 casos")


def test_gadget_xor_defectuoso_es_infactible():
    """Reproduce el defecto corregido: el gadget previo fuerza a = b."""
    import pulp
    infactibles = 0
    for va, vb in [(0, 1), (1, 0)]:
        prob = pulp.LpProblem("t", pulp.LpMinimize)
        a = pulp.LpVariable("a", cat="Binary")
        b = pulp.LpVariable("b", cat="Binary")
        aux = pulp.LpVariable("aux", cat="Binary")
        prob += a == va
        prob += b == vb
        prob += a + b - 2 * aux <= 0
        prob += a - b - 2 * aux <= 0
        prob += -a + b - 2 * aux <= 0
        prob += a + b - 2 * aux >= 0
        prob += 0
        prob.solve(pulp.PULP_CBC_CMD(msg=0))
        if pulp.LpStatus[prob.status] == "Infeasible":
            infactibles += 1
    assert infactibles == 2
    print("  [ok] el gadget defectuoso es infactible cuando a != b")


# ---------------------------------------------------------------------------
def _rango_gf2(filas, n):
    filas = filas[:]
    r = 0
    for col in range(n):
        piv = next((i for i in range(r, len(filas)) if (filas[i] >> col) & 1), None)
        if piv is None:
            continue
        filas[r], filas[piv] = filas[piv], filas[r]
        for i in range(len(filas)):
            if i != r and (filas[i] >> col) & 1:
                filas[i] ^= filas[r]
        r += 1
    return r


def _matriz_capa_lineal(w):
    """Matriz de lambda = pi o rho o theta sobre GF(2), filas como bitmasks."""
    idx = lambda x, y, z: x * (L * w) + y * w + z
    n = L * L * w
    M = [0] * n
    for x in range(L):
        for y in range(L):
            a = (x + 3 * y) % L
            sh = ROT[a][x] % w
            for z in range(w):
                c = (z - sh) % w
                m = 1 << idx(a, x, c)
                for yy in range(L):
                    m ^= 1 << idx((a - 1) % L, yy, c)
                for yy in range(L):
                    m ^= 1 << idx((a + 1) % L, yy, (c - 1) % w)
                M[idx(x, y, z)] = m
    return M, n


def test_capa_lineal_invertible():
    """lambda tiene rango pleno sobre GF(2) para z = 4 y z = 8."""
    for w in (4, 8):
        M, n = _matriz_capa_lineal(w)
        assert _rango_gf2(M, n) == n, "lambda singular para z=%d" % w
        print("  [ok] z=%d: lambda invertible (rango %d de %d)" % (w, n, n))


def test_degeneracion_offsets_rho():
    """El reescalado de rho mod z degrada la difusión.

    z=4: 4 valores distintos y 7 carriles sin rotar.
    z=8: 8 valores distintos y 3 carriles sin rotar.
    """
    esperado = {4: (4, 7), 8: (8, 3)}
    for w in (4, 8):
        planos = [ROT[x][y] % w for x in range(L) for y in range(L)]
        distintos, ceros = len(set(planos)), planos.count(0)
        assert (distintos, ceros) == esperado[w]
        print("  [ok] z=%d: %d offsets distintos, %d carriles sin rotar"
              % (w, distintos, ceros))


# ---------------------------------------------------------------------------
def test_variable_dinamica():
    """El contador de intentos determina las rondas según lo especificado."""
    assert intentos_a_rondas(0) == 1
    assert intentos_a_rondas(9) == 1
    assert intentos_a_rondas(10) == 2
    assert intentos_a_rondas(19) == 2
    assert intentos_a_rondas(20) == 3
    assert intentos_a_rondas(29) == 3
    print("  [ok] intentos -> rondas: 1 (<10), 2 (10-19), 3 (>=20)")


def test_cotas_superiores_reproducibles():
    """La búsqueda dirigida reproduce las cotas publicadas."""
    from milp_keccak import cota_superior
    esperado = {(1, 4): 1, (2, 4): 4, (3, 4): 18,
                (1, 8): 1, (2, 8): 4, (3, 8): 20}
    for (R, z), valor in sorted(esperado.items()):
        n, _ = cota_superior(R, z)
        assert n == valor, "R=%d z=%d: se obtuvo %d, se esperaba %d" % (R, z, n, valor)
    print("  [ok] cotas superiores reproducidas: 1/4/18 (z=4) y 1/4/20 (z=8)")


def test_cota_inferior_estructural():
    """La cota inferior R nunca supera a la cota superior hallada."""
    from milp_keccak import cota_superior
    for z in (4, 8):
        for R in (1, 2, 3):
            n, por_ronda = cota_superior(R, z)
            assert n >= R
            assert all(c >= 1 for c in por_ronda), "una ronda sin cajas activas"
    print("  [ok] cada ronda aporta >= 1 caja activa (cota inferior R)")


def test_cpsat_certifica_optimos():
    """CP-SAT certifica los óptimos de R = 1 y R = 2 (casos rápidos).

    Se omite si OR-Tools no está instalado. Los casos R = 3 no se incluyen aquí
    por tiempo de ejecución: véase docs/RESULTADOS.md.
    """
    try:
        from cpsat_keccak import resolver
    except ImportError as e:
        if "undefined symbol" in str(e):
            print("  [omitido] conflicto highspy/ortools en este proceso: "
                  "ejecutar esta prueba por separado")
        else:
            print("  [omitido] OR-Tools no instalado (pip install ortools)")
        return
    esperado = {(1, 4): 1, (1, 8): 1, (2, 4): 4, (2, 8): 4}
    for (R, z), valor in sorted(esperado.items()):
        d = resolver(R, z, limite=60)
        assert d["certificado"], "R=%d z=%d no certificado" % (R, z)
        assert d["cota_superior"] == valor, (
            "R=%d z=%d: se obtuvo %s, se esperaba %d" % (R, z, d["cota_superior"], valor))
    print("  [ok] CP-SAT certifica: R=1 -> 1 y R=2 -> 4 (z = 4 y 8)")


def test_cotas_voraces_no_son_optimas():
    """La búsqueda voraz da cotas correctas pero no óptimas para R = 3.

    Documenta una limitación real del método: elegir en cada caja la salida de
    menor peso de Hamming es una decisión local y no alcanza el óptimo global
    (18 frente al óptimo certificado 9 para R = 3, z = 4).
    """
    from milp_keccak import cota_superior
    n, _ = cota_superior(3, 4)
    assert n == 18, "la cota voraz cambió: %d" % n
    assert n > 9, "la cota voraz debería ser peor que el óptimo certificado"
    print("  [ok] cota voraz R=3 z=4 = 18 > 9 (óptimo certificado por CP-SAT)")

# ---------------------------------------------------------------------------
if __name__ == "__main__":
    pruebas = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    print("Ejecutando %d verificaciones\n" % len(pruebas))
    fallos = 0
    for fn in pruebas:
        print("%s:" % fn.__name__)
        try:
            fn()
        except AssertionError as e:
            fallos += 1
            print("  [FALLO] %s" % e)
    print("\n%d/%d verificaciones superadas" % (len(pruebas) - fallos, len(pruebas)))
    sys.exit(1 if fallos else 0)
