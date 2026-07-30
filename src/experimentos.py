#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ejecuta los experimentos de Keccak dinámico con MILP, ahora con soporte para
R >= 4, paralelización y estrategia de decisión para certificar optimalidad.

Uso:
    python src/experimentos.py                    # todos los casos (1..5 rondas, z=4,8)
    python src/experimentos.py --tiempo 300       # más tiempo por caso
    python src/experimentos.py --solo-trayectorias  # sin MILP (solo cotas superiores)
    python src/experimentos.py --z 4 --rondas 5   # un caso concreto
    python src/experimentos.py --workers 4        # usar 4 procesos en paralelo
"""

import argparse
import json
import os
import sys
import concurrent.futures
import multiprocessing

# Forzamos spawn para evitar problemas con highspy en Linux/Windows
if __name__ == "__main__":
    multiprocessing.set_start_method('spawn', force=True)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from milp_keccak import analizar_con_decision, cota_superior, DP_MAX_LOG2


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


def ejecutar_caso(z, R, solo_trayectorias, tiempo):
    """Función independiente para ejecutar un caso (z, R) en un proceso separado."""
    if solo_trayectorias:
        ub, por_ronda = cota_superior(R, z)
        peso = -DP_MAX_LOG2 * ub
        d = {
            "R": R, "z": z, "cajas_por_ronda": 5 * z,
            "cota_inferior": R, "cota_superior": ub,
            "por_ronda": por_ronda, "certificado": False,
            "estado_solver": "sólo trayectoria",
            "prob_log2": -peso, "pares_log2": peso,
        }
    else:
        # Usamos la nueva función con estrategia de decisión
        d = analizar_con_decision(R, z, limite_tiempo=tiempo)
    return d


def main():
    ap = argparse.ArgumentParser(description="Experimentos MILP de Keccak dinámico (paralelo)")
    ap.add_argument("--tiempo", type=int, default=60,
                    help="límite de tiempo del MILP por caso, en segundos")
    ap.add_argument("--solo-trayectorias", action="store_true",
                    help="omitir el MILP y calcular sólo las cotas superiores")
    ap.add_argument("--z", type=int, choices=[4, 8], help="ejecutar un solo z")
    ap.add_argument("--rondas", type=int, help="ejecutar un solo R (cualquier entero >=1)")
    ap.add_argument("--salida", default="resultados.json")
    ap.add_argument("--workers", type=int, default=None,
                    help="número de procesos en paralelo (por defecto: número de núcleos)")
    args = ap.parse_args()

    zs = [args.z] if args.z else [4, 8]
    # Si no se especifica --rondas, ejecutamos de 1 a 5 (puedes ampliar)
    if args.rondas:
        rs = [args.rondas]
    else:
        rs = list(range(1, 11))  # 1,2,3,4,5

    casos = [(z, R) for z in zs for R in rs]
    print(f"Ejecutando {len(casos)} casos en paralelo con {args.workers or 'auto'} workers...")

    filas = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(ejecutar_caso, z, R, args.solo_trayectorias, args.tiempo): (z, R)
            for z, R in casos
        }
        for future in concurrent.futures.as_completed(futures):
            z, R = futures[future]
            try:
                resultado = future.result()
                filas.append(resultado)
                print(f"✅ Completado: z={z}, R={R}")
            except Exception as e:
                print(f"❌ Error en z={z}, R={R}: {e}")
                filas.append({"error": str(e), "z": z, "R": R})

    filas.sort(key=lambda d: (d.get("z", 0), d.get("R", 0)))
    print()
    tabla(filas)
    print("\nLa probabilidad de la trayectoria se acota por 2^-2n y el número de pares")
    print("necesarios por 2^2n, siendo n el número de cajas-S activas (DP_max de chi = 2^-2).")

    with open(args.salida, "w", encoding="utf-8") as f:
        json.dump(filas, f, indent=2, ensure_ascii=False)
    print("\nResultados guardados en %s" % args.salida)


if __name__ == "__main__":
    main()