#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fusiona los resultados de varias ejecuciones individuales de `cpsat_keccak.py`.

Cada ejecución guarda su propio JSON (una por caso, para no perder trabajo si algo
se interrumpe). Esta utilidad los une, se queda con lo mejor de cada caso —la cota
inferior más alta y la superior más baja entre todas las corridas—, aplica la
consolidación por monotonía y muestra la tabla final.

Uso:
    python src/fusionar.py resultados/*.json
    python src/fusionar.py resultados/*.json --salida final.json

En PowerShell los comodines no se expanden solos, de modo que conviene indicar la
carpeta:
    py src\\fusionar.py resultados --salida final.json
"""

import argparse
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cpsat_keccak import consolidar, magnitudes


def cargar(rutas):
    """Acepta ficheros, comodines o carpetas."""
    ficheros = []
    for r in rutas:
        if os.path.isdir(r):
            ficheros.extend(sorted(glob.glob(os.path.join(r, "*.json"))))
        else:
            expandido = sorted(glob.glob(r))
            ficheros.extend(expandido if expandido else [r])

    registros = []
    for f in ficheros:
        try:
            with open(f, encoding="utf-8") as fh:
                d = json.load(fh)
            registros.extend(d if isinstance(d, list) else [d])
        except (OSError, json.JSONDecodeError) as e:
            print("  [aviso] no se pudo leer %s (%s)" % (f, type(e).__name__))
    return ficheros, registros


def mejor_por_caso(registros):
    """Combina varias corridas del mismo caso quedándose con las mejores cotas."""
    mejor = {}
    for d in registros:
        clave = (d["z"], d["R"])
        if clave not in mejor:
            mejor[clave] = dict(d)
            continue
        m = mejor[clave]
        # la cota inferior más alta demostrada
        if d.get("cota_inferior", 0) > m.get("cota_inferior", 0):
            m["cota_inferior"] = d["cota_inferior"]
        # la cota superior más baja hallada
        ub_d, ub_m = d.get("cota_superior"), m.get("cota_superior")
        if ub_d is not None and (ub_m is None or ub_d < ub_m):
            m["cota_superior"] = ub_d
            m["por_ronda"] = d.get("por_ronda")
        # una certificación previa no se pierde
        if d.get("certificado"):
            m["certificado"] = True
            m["estado"] = d.get("estado", m.get("estado"))
        m["segundos"] = max(m.get("segundos", 0) or 0, d.get("segundos", 0) or 0)
    return list(mejor.values())


def main():
    ap = argparse.ArgumentParser(description="Fusiona resultados de cpsat_keccak.py")
    ap.add_argument("rutas", nargs="+", help="ficheros JSON, comodines o una carpeta")
    ap.add_argument("--salida", default=None, help="guardar el resultado fusionado")
    a = ap.parse_args()

    ficheros, registros = cargar(a.rutas)
    if not registros:
        print("No se encontraron resultados en:", ", ".join(a.rutas))
        return 1
    print("Leídos %d fichero(s), %d registro(s)\n" % (len(ficheros), len(registros)))

    res = consolidar(mejor_por_caso(registros))
    for d in res:
        d.update(magnitudes(d))

    cab = " z |  R | mínimo        | por ronda            | estado    | pares (garant.)"
    print(cab); print("-" * len(cab))
    certificados = 0
    for d in sorted(res, key=lambda x: (x["z"], x["R"])):
        lb = d["cota_inferior_consolidada"]
        ub = d["cota_superior_consolidada"]
        if d["exacto"]:
            minimo = "%d (exacto)" % ub; certificados += 1
        elif ub is None:
            minimo = ">= %d" % lb
        else:
            minimo = "[%d, %d]" % (lb, ub)
        pr = ",".join(map(str, d["por_ronda"])) if d.get("por_ronda") else "-"
        if len(pr) > 20:
            pr = pr[:17] + "..."
        print(" %d | %2d | %-13s | %-20s | %-9s | 2^%d"
              % (d["z"], d["R"], minimo, pr, d.get("estado", "-"),
                 d["pares_log2_garantizados"]))

    print("\n%d de %d casos con óptimo exacto." % (certificados, len(res)))
    print("Los pares se derivan de la cota INFERIOR consolidada.")

    if a.salida:
        with open(a.salida, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=2, ensure_ascii=False)
        print("Guardado en %s" % a.salida)
    return 0


if __name__ == "__main__":
    sys.exit(main())
