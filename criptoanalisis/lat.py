#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Criptoanalisis lineal: aproximaciones afines, tabla LAT y DDT.

Ejercicios 1 a 3 de la actividad.
"""
S = [0xc,0x2,0x6,0xf,0x9,0x0,0x4,0xd,0x3,0xe,0xb,0x8,0xa,0x7,0x1,0x5]
par = lambda v: bin(v).count("1") & 1


def ejercicio1():
    """Aproximar x OR y con las funciones afines de dos variables."""
    entradas = [(0,0),(0,1),(1,0),(1,1)]
    OR = {(x,y): x|y for x,y in entradas}
    funcs = [("1", lambda x,y: 1), ("x", lambda x,y: x),
             ("y", lambda x,y: y), ("x+y", lambda x,y: x^y)]
    aciertos = [sum(1 for x,y in entradas if f(x,y)==OR[(x,y)]) for _,f in funcs]
    return [(n, a, a/4, abs(a/4-0.5)) for (n,_),a in zip(funcs, aciertos)]


def ejercicio2():
    """Aproximacion y3 = x2 (+) x3 sobre la S-box de 4 bits."""
    filas, ok = [], 0
    for x in range(16):
        x2, x3 = (x>>2)&1, (x>>1)&1
        y3 = (S[x]>>1)&1
        c = (y3 == (x2 ^ x3)); ok += c
        filas.append((x, S[x], x2^x3, y3, c))
    return filas, ok


def lat():
    """LAT[a][b] = #{x : a.x = b.S(x)} - 8."""
    return [[sum(1 for x in range(16) if par(a&x)==par(b&S[x])) - 8
             for b in range(16)] for a in range(16)]


def ddt():
    D = [[0]*16 for _ in range(16)]
    for da in range(16):
        for x in range(16):
            D[da][S[x]^S[x^da]] += 1
    return D


if __name__ == "__main__":
    print("=== 1. Aproximaciones afines de x OR y ===")
    for n,a,p,s in ejercicio1():
        print("   %-4s : %d/4 = %.2f   sesgo %.2f" % (n,a,p,s))
    print("   Todas aciertan 3/4: la no linealidad de OR es 1.\n")

    print("=== 2. y3 = x2 (+) x3 ===")
    filas, ok = ejercicio2()
    print("   coincidencias: %d/16 = %.4f   sesgo %.4f   LAT = %+d"
          % (ok, ok/16, abs(ok/16-0.5), ok-8))
    print("   textos necesarios ~ 1/sesgo^2 = %d\n" % round(1/(abs(ok/16-0.5)**2)))

    T = lat()
    print("=== 3. LAT ===")
    h = "0123456789abcdef"
    print("   a\\b " + "".join("%4s"%h[b] for b in range(16)))
    for a in range(16):
        print("    %s  "%h[a] + "".join("%4d"%T[a][b] for b in range(16)))
    mx = max(abs(T[a][b]) for a in range(16) for b in range(16) if (a,b)!=(0,0))
    D = ddt(); dmx = max(D[a][b] for a in range(1,16) for b in range(16))
    print("\n   |LAT| maximo = %d -> p = %.4f, sesgo = %.4f" % (mx,(8+mx)/16,mx/16))
    print("   DDT maximo   = %d" % dmx)
    print("   (una S-box de 4 bits optima alcanza 4 en ambas metricas)")
