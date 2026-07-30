#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modelo CP-SAT (OR-Tools) para Keccak dinámico.

Este archivo es OPCIONAL y se proporciona como comparación con el modelo MILP.
El modelo MILP con Convex Hull ya certifica optimalidad en todos los casos;
CP-SAT fue útil durante el desarrollo para validar resultados en R=1..3.
"""

import math
import time
from ortools.sat.python import cp_model

try:
    from sha3_nucleo import ROT, ddt_transiciones
except ImportError:
    from src.sha3_nucleo import ROT, ddt_transiciones

L = 5
DP_MAX_LOG2 = -2
DDT = ddt_transiciones()   # diccionario a -> lista de b válidos

def resolver_cpsat(R, z, limite_tiempo=60, cota_superior=None):
    """
    Resuelve el problema con CP-SAT (OR-Tools).
    Devuelve (certificado, incumbente, cota_dual, estado, por_ronda).
    """
    model = cp_model.CpModel()
    
    # Variables booleanas
    D = {}
    for r in range(R + 1):
        for x in range(L):
            for y in range(L):
                for k in range(z):
                    D[(r, x, y, k)] = model.NewBoolVar(f'D_{r}_{x}_{y}_{k}')
    
    A = {}
    for r in range(R):
        for y in range(L):
            for k in range(z):
                A[(r, y, k)] = model.NewBoolVar(f'A_{r}_{y}_{k}')
    
    # Auxiliares para XOR (igual que en MILP)
    def xor_lineal(a, b, out):
        t = model.NewBoolVar('t')
        model.Add(a + b - 2 * t == out)
    
    # --- Capas lineales (idénticas al MILP) ---
    for r in range(R):
        # theta: paridades
        C = {}
        for x in range(L):
            for k in range(z):
                acc = D[(r, x, 0, k)]
                for y in range(1, L):
                    nxt = model.NewBoolVar(f'C_{r}_{x}_{k}_{y}')
                    xor_lineal(acc, D[(r, x, y, k)], nxt)
                    acc = nxt
                C[(x, k)] = acc
        
        Dth = {}
        for x in range(L):
            for k in range(z):
                var = model.NewBoolVar(f'Dth_{r}_{x}_{k}')
                xor_lineal(C[((x - 1) % L, k)], C[((x + 1) % L, (k - 1) % z)], var)
                Dth[(x, k)] = var
        
        Dt = {}
        for x in range(L):
            for y in range(L):
                for k in range(z):
                    var = model.NewBoolVar(f'Dt_{r}_{x}_{y}_{k}')
                    xor_lineal(D[(r, x, y, k)], Dth[(x, k)], var)
                    Dt[(x, y, k)] = var
        
        # rho + pi (reindexación)
        Drp = {}
        for x in range(L):
            for y in range(L):
                nx, ny, rot = y, (2 * x + 3 * y) % L, ROT[x][y] % z
                for k in range(z):
                    Drp[(nx, ny, (k + rot) % z)] = Dt[(x, y, k)]
        
        # --- chi (con AddAllowedAssignments) ---
        for y in range(L):
            for k in range(z):
                vin_vars = [D[(r, i, y, k)] for i in range(L)]
                vout_vars = [D[(r + 1, i, y, k)] for i in range(L)]
                
                # CORRECCIÓN: usar listas en lugar de generadores
                vin = cp_model.LinearExpr.Sum([vin_vars[i] * (1 << i) for i in range(L)])
                vout = cp_model.LinearExpr.Sum([vout_vars[i] * (1 << i) for i in range(L)])
                act = A[(r, y, k)]
                
                # Construir tabla de transiciones válidas (in, out, act)
                transiciones = []
                for a in range(32):
                    for b in DDT[a]:
                        act_val = 1 if a != 0 else 0
                        transiciones.append((a, b, act_val))
                # Aseguramos que el caso inactivo (0,0,0) esté incluido
                if (0, 0, 0) not in transiciones:
                    transiciones.append((0, 0, 0))
                
                model.AddAllowedAssignments([vin, vout, act], transiciones)
    
    # Objetivo: minimizar total de cajas activas
    total_cajas = sum(A[(r, y, k)] for r in range(R) for y in range(L) for k in range(z))
    model.Minimize(total_cajas)
    
    # No trivialidad: al menos una diferencia en la entrada
    model.Add(sum(D[(0, x, y, k)] for x in range(L) for y in range(L) for k in range(z)) >= 1)
    
    # Si se da una cota superior, usamos estrategia de decisión
    if cota_superior is not None:
        model.Add(total_cajas <= cota_superior - 1)
    
    # Resolver
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = limite_tiempo
    solver.parameters.log_search_progress = False
    status = solver.Solve(model)
    
    # Interpretar resultado
    if status == cp_model.OPTIMAL:
        certificado = True
        incumbente = int(solver.ObjectiveValue())
        dual = incumbente  # CP-SAT no da cota dual, la aproximamos
        estado = "kOptimal"
    elif status == cp_model.INFEASIBLE:
        certificado = True  # infactible certifica que no hay solución con la cota
        incumbente = None
        dual = float('inf')
        estado = "kInfeasible"
    else:
        certificado = False
        incumbente = int(solver.ObjectiveValue()) if solver.HasObjectiveValue() else None
        dual = float('nan')
        estado = "kTimeLimit"
    
    # Extraer distribución por ronda si hay solución
    por_ronda = []
    if incumbente is not None:
        for r in range(R):
            cnt = sum(1 for y in range(L) for k in range(z) 
                      if solver.Value(A[(r, y, k)]) > 0.5)
            por_ronda.append(cnt)
    else:
        por_ronda = [0] * R
    
    return certificado, incumbente, dual, estado, por_ronda


if __name__ == "__main__":
    print("Prueba CP-SAT: R=3, z=4")
    cert, inc, dual, estado, pr = resolver_cpsat(3, 4, limite_tiempo=60)
    print(f"Estado: {estado}, Incumbente: {inc}, Por ronda: {pr}")