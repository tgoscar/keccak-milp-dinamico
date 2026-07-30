#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ejecuta los seis experimentos del estudio: z en {4, 8} x R en {1, 2, 3}.

Uso:
    python src/experimentos.py                    # los 6 casos, 60 s de MILP cada uno
    python src/experimentos.py --tiempo 300       # más tiempo por caso
    python src/experimentos.py --solo-trayectorias  # sin MILP (segundos)
    python src/experimentos.py --z 4 --rondas 2   # un caso concreto

Salida: tabla por consola y `resultados.json` con todos los campos.

Nota: para R >= 2 el MILP no cierra el gap (la relajación lineal es demasiado
débil), de modo que aumentar el tiempo rara vez cambia el resultado. La opción
--solo-trayectorias reproduce las cotas superiores publicadas en segundos.
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from milp_keccak import analizar, cota_superior, DP_MAX_LOG2


def tabla(filas):
    cab = (" z | R | S-boxes/ronda | cota inf | cota sup | por ronda | "
           "prob.     | pares    | estado")
    print(cab)
    print("-" * len(cab))
    for d in filas:
        estado = ("óptimo certificado" if d.get("certificado")
                  else d.get("estado_solver", "sólo trayectoria"))
        print(" %d | %d |      %3d      |   %3d    |   %3d    | %-9s | 2^%-6d | 2^%-6d | %s"
              % (d["z"], d["R"], d["cajas_por_ronda"], d["cota_inferior"],
                 d["cota_superior"], ",".join(map(str, d["por_ronda"])),
                 d["prob_log2"], d["pares_log2"], estado))


def main():
    ap = argparse.ArgumentParser(description="Experimentos MILP de Keccak dinámico")
    ap.add_argument("--tiempo", type=int, default=60,
                    help="límite de tiempo del MILP por caso, en segundos")
    ap.add_argument("--solo-trayectorias", action="store_true",
                    help="omitir el MILP y calcular sólo las cotas superiores")
    ap.add_argument("--z", type=int, choices=[4, 8], help="ejecutar un solo z")
    ap.add_argument("--rondas", type=int, choices=[1, 2, 3], help="ejecutar un solo R")
    ap.add_argument("--salida", default="resultados.json")
    args = ap.parse_args()

    zs = [args.z] if args.z else [4, 8]
    rs = [args.rondas] if args.rondas else [1, 2, 3]

    filas = []
    for z in zs:
        for R in rs:
            if args.solo_trayectorias:
                ub, por_ronda = cota_superior(R, z)
                peso = -DP_MAX_LOG2 * ub
                d = {"R": R, "z": z, "cajas_por_ronda": 5 * z,
                     "cota_inferior": R, "cota_superior": ub,
                     "por_ronda": por_ronda, "certificado": False,
                     "estado_solver": "sólo trayectoria",
                     "prob_log2": -peso, "pares_log2": peso}
            else:
                d = analizar(R, z, limite_tiempo=args.tiempo)
            filas.append(d)

    print()
    tabla(filas)
    print("\nLa probabilidad de la trayectoria se acota por 2^-2n y el número de "
          "pares\nnecesarios por 2^2n, siendo n el número de cajas-S activas "
          "(DP_max de chi = 2^-2).")

    with open(args.salida, "w", encoding="utf-8") as f:
        json.dump(filas, f, indent=2, ensure_ascii=False)
    print("\nResultados guardados en %s" % args.salida)


if __name__ == "__main__":
    main()
