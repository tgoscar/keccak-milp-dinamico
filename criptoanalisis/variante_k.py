#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Familia parametrizada: k rotaciones por carril.

    L_i <- L_i (+) rot(L_i, r_{i,1}) (+) ... (+) rot(L_i, r_{i,k-1})

Hipotesis: la cota diferencial la fija k (cuantas columnas alcanza un bit
activo), no la eleccion fina de las constantes. Con k terminos, un bit activo
llega a k columnas, luego la trayectoria minima crece con k.
Coste: (k-1) palabras XOR por carril, 5(k-1) por ronda.
"""
import time
from ortools.sat.python import cp_model
L = 5
def sbox_chi(x):
    b=[(x>>i)&1 for i in range(5)]
    return sum((b[i]^((~b[(i+1)%5])&b[(i+2)%5])&1)<<i for i in range(5))
DDT={a:sorted({sbox_chi(x)^sbox_chi(x^a) for x in range(32)}) for a in range(32)}
TUPLAS=[[(a>>i)&1 for i in range(5)]+[(b>>i)&1 for i in range(5)]
        for a,bl in DDT.items() for b in bl]

def construir(R,w,rots):
    """rots[i] = tupla de desplazamientos (sin incluir el 0, que va implicito)."""
    m=cp_model.CpModel()
    D={(r,i,z):m.NewBoolVar("D%d_%d_%d"%(r,i,z))
       for r in range(R+1) for i in range(L) for z in range(w)}
    A={}
    for r in range(R):
        T={}
        for i in range(L):
            for z in range(w):
                t=m.NewBoolVar("T%d_%d_%d"%(r,i,z))
                lits=[t.Not(), D[(r,i,z)]]+[D[(r,i,(z-s)%w)] for s in rots[i]]
                m.AddBoolXOr(lits); T[(i,z)]=t
        for z in range(w):
            ent=[T[(i,z)] for i in range(L)]; sal=[D[(r+1,i,z)] for i in range(L)]
            m.AddAllowedAssignments(ent+sal,TUPLAS)
            a=m.NewBoolVar("A%d_%d"%(r,z)); m.AddMaxEquality(a,ent); A[(r,z)]=a
    m.AddBoolOr([D[(0,i,0)] for i in range(L)])
    m.Minimize(sum(A.values()))
    return m,A

def resolver(R,w,rots,limite=120,hilos=8):
    m,A=construir(R,w,rots)
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=float(limite)
    s.parameters.num_search_workers=hilos
    t0=time.time(); est=s.Solve(m); dt=time.time()-t0
    hay=est in (cp_model.OPTIMAL,cp_model.FEASIBLE)
    pr=[sum(s.Value(v) for (rr,z),v in A.items() if rr==r) for r in range(R)] if hay else None
    return {"estado":s.StatusName(est),"cert":est==cp_model.OPTIMAL,
            "lb":int(s.BestObjectiveBound()) if hay else 0,
            "ub":int(s.ObjectiveValue()) if hay else None,"pr":pr,"seg":round(dt,1)}
