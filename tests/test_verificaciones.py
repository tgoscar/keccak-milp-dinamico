#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pruebas unitarias para la variante de Keccak y el modelo MILP.
Verifica que los resultados óptimos coincidan con los esperados.
"""

import sys
import os
import unittest
import json

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.milp_keccak import analizar_con_decision, cota_superior
from src.sha3_nucleo import ddt_transiciones

class TestKeccakMILP(unittest.TestCase):
    
    def test_ddt_correcta(self):
        """Verifica que la DDT tenga todas las filas no vacías y la multiplicidad máxima."""
        DDT = ddt_transiciones()
        vacias = [a for a in range(32) if not DDT[a]]
        self.assertEqual(vacias, [], f"DDT tiene filas vacías: {vacias}")
        
        # Verificar multiplicidad máxima = 8 (para a != 0)
        from src.sha3_nucleo import ddt_completa
        cnt = [[0]*32 for _ in range(32)]
        for a in range(32):
            for x in range(32):
                from src.sha3_nucleo import sbox_chi
                cnt[a][sbox_chi(x) ^ sbox_chi(x ^ a)] += 1
        maxmult = max(cnt[a][b] for a in range(1, 32) for b in range(32))
        self.assertEqual(maxmult, 8, f"Multiplicidad máxima es {maxmult}, debería ser 8")
    
    def test_cota_inferior_estructural(self):
        """La cota inferior R siempre se cumple (cada ronda aporta ≥1 caja activa)."""
        for z in (4, 8):
            for R in range(1, 11):
                # La búsqueda de trayectorias siempre debe devolver al menos R
                ub, _ = cota_superior(R, z)
                self.assertGreaterEqual(ub, R, f"R={R}, z={z}: cota superior {ub} < R")
    
    def test_optimalidad_R1(self):
        """R=1 siempre debe ser 1 (óptimo certificado)."""
        for z in (4, 8):
            d = analizar_con_decision(1, z, limite_tiempo=10)
            self.assertEqual(d["cota_superior"], 1)
            self.assertEqual(d["cota_inferior"], 1)
            self.assertTrue(d["certificado"] or d["estado_solver"] == "kInfeasible")
    
    def test_optimalidad_R2_R3(self):
        """Para R=2 y R=3 el óptimo certificado debe ser R (2 y 3 respectivamente)."""
        expected = {2: 2, 3: 3}
        for z in (4, 8):
            for R in (2, 3):
                d = analizar_con_decision(R, z, limite_tiempo=60)
                self.assertEqual(d["cota_superior"], expected[R])
                self.assertEqual(d["cota_inferior"], expected[R])
                self.assertTrue(d["certificado"])
    
    def test_optimalidad_hasta_R10(self):
        """Para R=4..10 el óptimo certificado debe ser R."""
        for z in (4, 8):
            for R in range(4, 11):
                d = analizar_con_decision(R, z, limite_tiempo=60)
                self.assertEqual(d["cota_superior"], R)
                self.assertEqual(d["cota_inferior"], R)
                self.assertTrue(d["certificado"])
    
    def test_consistencia_resultados_json(self):
        """Comprueba que el JSON generado por experimentos tenga datos coherentes."""
        json_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "resultados.json")
        if not os.path.exists(json_path):
            self.skipTest("resultados.json no encontrado, saltando prueba")
        with open(json_path, "r") as f:
            datos = json.load(f)
        
        for d in datos:
            z = d["z"]
            R = d["R"]
            self.assertIn(z, (4, 8))
            self.assertIn(R, range(1, 11))
            self.assertEqual(d["cota_inferior"], R)
            self.assertEqual(d["cota_superior"], R)
            self.assertTrue(d["certificado"] or (R == 1 and d["estado_solver"] == "kInfeasible"))
            self.assertEqual(d["prob_log2"], -2 * R)
            self.assertEqual(d["pares_log2"], 2 * R)


if __name__ == "__main__":
    unittest.main(verbosity=2)