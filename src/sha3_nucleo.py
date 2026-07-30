#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Núcleo Keccak / SHA-3 (FIPS 202) y variante de rondas dinámicas.

Contiene:
  * Los cinco pasos de Keccak como funciones puras (theta, rho, pi, chi, iota).
  * La permutación estándar de 24 rondas y la esponja SHA-3-256, validadas
    contra `hashlib`.
  * La variante dinámica: el número de rondas depende del contador de intentos.
  * La caja-S chi de 5 bits y su DDT (Tabla de Distribución de Diferencias)
    calculada correctamente.

Referencia: NIST FIPS PUB 202, "SHA-3 Standard", 2015.
"""

import hashlib

MASK64 = 0xFFFFFFFFFFFFFFFF

# Desplazamientos de rho, indexados ROT[x][y] (especificados para w = 64).
ROT = [
    [0, 36, 3, 41, 18],
    [1, 44, 10, 45, 2],
    [62, 6, 43, 15, 61],
    [28, 55, 25, 21, 56],
    [27, 20, 39, 8, 14],
]

# Constantes de ronda iota para Keccak-f (24 rondas).
RC = [
    0x0000000000000001, 0x0000000000008082, 0x800000000000808A, 0x8000000080008000,
    0x000000000000808B, 0x0000000080000001, 0x8000000080008081, 0x8000000000008009,
    0x000000000000008A, 0x0000000000000088, 0x0000000080008009, 0x000000008000000A,
    0x000000008000808B, 0x800000000000008B, 0x8000000000008089, 0x8000000000008003,
    0x8000000000008002, 0x8000000000000080, 0x000000000000800A, 0x800000008000000A,
    0x8000000080008081, 0x8000000000008080, 0x0000000080000001, 0x8000000080008008,
]


# ---------------------------------------------------------------------------
# Los cinco pasos de Keccak (w = 64)
# ---------------------------------------------------------------------------
def ROTL64(x, n):
    """Rotación circular a la izquierda sobre 64 bits."""
    n %= 64
    if n == 0:
        return x & MASK64
    return ((x << n) & MASK64) | ((x & MASK64) >> (64 - n))


def theta(A):
    """Mezcla de columnas: cada bit se combina con dos paridades de columna."""
    C = [A[x][0] ^ A[x][1] ^ A[x][2] ^ A[x][3] ^ A[x][4] for x in range(5)]
    D = [C[(x - 1) % 5] ^ ROTL64(C[(x + 1) % 5], 1) for x in range(5)]
    return [[A[x][y] ^ D[x] for y in range(5)] for x in range(5)]


def rho(A):
    """Rotación de cada carril según su desplazamiento."""
    return [[ROTL64(A[x][y], ROT[x][y]) for y in range(5)] for x in range(5)]


def pi_(A):
    """Permutación de posiciones: (x, y) -> (y, 2x + 3y mod 5)."""
    B = [[0] * 5 for _ in range(5)]
    for x in range(5):
        for y in range(5):
            B[y][(2 * x + 3 * y) % 5] = A[x][y]
    return B


def chi(A):
    """Única capa no lineal: opera sobre filas de 5 bits a lo largo de x."""
    return [[(A[x][y] ^ ((~A[(x + 1) % 5][y]) & A[(x + 2) % 5][y])) & MASK64
             for y in range(5)] for x in range(5)]


def iota(A, rc):
    """Suma de la constante de ronda al carril (0, 0)."""
    B = [fila[:] for fila in A]
    B[0][0] ^= rc
    return B


def ronda(estado, rc):
    """Una ronda completa de Keccak-f."""
    return iota(chi(pi_(rho(theta(estado)))), rc)


def keccak_f(estado, rondas=24):
    """Permutación Keccak-f[1600]. Con `rondas` < 24 se obtiene la variante."""
    for i in range(rondas):
        estado = ronda(estado, RC[i])
    return estado


# ---------------------------------------------------------------------------
# Construcción de esponja SHA-3-256
# ---------------------------------------------------------------------------
def esponja_256(mensaje: bytes, rondas=24) -> bytes:
    """Esponja SHA-3-256 (rate = 1088 bits = 136 bytes)."""
    rate = 136
    estado = [[0] * 5 for _ in range(5)]

    # Relleno pad10*1 con el sufijo de dominio de SHA-3 (0x06 ... 0x80)
    p = bytearray(mensaje)
    p.append(0x06)
    while len(p) % rate != 0:
        p.append(0x00)
    p[-1] |= 0x80

    # Absorción
    for off in range(0, len(p), rate):
        bloque = p[off:off + rate]
        for i in range(rate // 8):
            estado[i % 5][i // 5] ^= int.from_bytes(bloque[i * 8:(i + 1) * 8], "little")
        estado = keccak_f(estado, rondas)

    # Exprimido: 32 bytes = 256 bits
    salida = bytearray()
    for i in range(4):
        salida += estado[i % 5][i // 5].to_bytes(8, "little")
    return bytes(salida)


def sha3_256(mensaje: bytes) -> bytes:
    """SHA-3-256 estándar (24 rondas)."""
    return esponja_256(mensaje, rondas=24)


# ---------------------------------------------------------------------------
# Variante dinámica: las rondas dependen del contador de intentos
# ---------------------------------------------------------------------------
def intentos_a_rondas(intentos: int) -> int:
    """Número de rondas en función de los intentos fallidos.

    intentos < 10 -> 1 ronda ; < 20 -> 2 rondas ; < 30 -> 3 rondas.
    A partir de 30 se satura en 3, que es el máximo analizado en el estudio.
    """
    if intentos < 10:
        return 1
    if intentos < 20:
        return 2
    return 3


def keccak_modificado(mensaje: bytes, intentos: int) -> bytes:
    """Variante dinámica: SHA-3-256 con un número de rondas según `intentos`."""
    return esponja_256(mensaje, rondas=intentos_a_rondas(intentos))


# ---------------------------------------------------------------------------
# Caja-S chi de 5 bits y su DDT
# ---------------------------------------------------------------------------
def sbox_chi(x: int) -> int:
    """chi sobre una fila de 5 bits: chi_i = x_i XOR ((NOT x_{i+1}) AND x_{i+2})."""
    b = [(x >> i) & 1 for i in range(5)]
    o = [b[i] ^ ((~b[(i + 1) % 5]) & b[(i + 2) % 5]) & 1 for i in range(5)]
    return sum(o[i] << i for i in range(5))


def ddt_completa():
    """DDT con multiplicidades: ddt[a][b] = #{x : S(x) XOR S(x^a) = b}."""
    ddt = [[0] * 32 for _ in range(32)]
    for a in range(32):
        for x in range(32):
            ddt[a][sbox_chi(x) ^ sbox_chi(x ^ a)] += 1
    return ddt


def ddt_transiciones():
    """DDT como diccionario a -> lista de diferencias de salida alcanzables.

    NOTA: no debe filtrarse el recorrido de x. Una implementación previa
    incluía la condición `if (x ^ S(x)) == a`, ajena a la definición de DDT,
    que dejaba 21 de 32 filas vacías (véase docs/CORRECCIONES.md).
    """
    return {a: sorted({sbox_chi(x) ^ sbox_chi(x ^ a) for x in range(32)})
            for a in range(32)}


def dp_maxima():
    """Probabilidad diferencial máxima de la caja-S para diferencias no nulas.

    Devuelve (multiplicidad_maxima, exponente) con prob = mult/32 = 2^exponente.
    Para chi resulta 8/32 = 2^-2.
    """
    ddt = ddt_completa()
    m = max(ddt[a][b] for a in range(1, 32) for b in range(32))
    from math import log2
    return m, int(log2(m / 32))


if __name__ == "__main__":
    # Validación del núcleo contra la biblioteca estándar
    vectores = [b"", b"abc", b"a" * 135, b"a" * 136, b"a" * 137,
                b"Plataforma del Estado - modulo de autenticacion"]
    ok = sum(sha3_256(m).hex() == hashlib.sha3_256(m).hexdigest() for m in vectores)
    print("SHA-3-256 frente a hashlib: %d/%d vectores coinciden (FIPS 202)"
          % (ok, len(vectores)))

    trans = ddt_transiciones()
    mult, expo = dp_maxima()
    print("DDT de chi: %d transiciones válidas, filas vacías: %d"
          % (sum(len(v) for v in trans.values()),
             sum(1 for a in range(32) if not trans[a])))
    print("Probabilidad diferencial máxima: %d/32 = 2^%d" % (mult, expo))
    print("DDT[1] =", trans[1])

    for i in (0, 5, 12, 25):
        print("intentos=%2d -> %d ronda(s)" % (i, intentos_a_rondas(i)))
