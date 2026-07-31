"""
Criterio de diseno: subespacio de diferencias que NO activan ninguna caja-S.

Si una diferencia d cumple que, en cada ronda r, su imagen tras la capa lineal
es nula en todas las columnas que llevan chi, entonces chi actua sobre cero y la
ronda entera es LINEAL para esa diferencia. La trayectoria tiene probabilidad 1.

La condicion es un sistema lineal sobre GF(2):
        P_r ( Lambda^(r+1) (d) ) = 0     para r = 0 .. R-1
donde P_r proyecta sobre las columnas activas de la ronda r. El diseno es seguro
frente a este fallo si y solo si el nucleo de ese sistema es trivial.
"""
L = 5

def matriz_lambda(w, rots):
    """Lambda como lista de filas (bitmask): fila (i,z) = bits de los que depende."""
    n = 5*w
    idx = lambda i,z: i*w + z
    M = [0]*n
    for i in range(L):
        a,b = rots[i]
        for z in range(w):
            M[idx(i,z)] = (1<<idx(i,z)) ^ (1<<idx(i,(z-a)%w)) ^ (1<<idx(i,(z-b)%w))
    return M, n

def aplicar(M, v, n):
    """y = M v sobre GF(2), v como bitmask."""
    y = 0
    for i in range(n):
        if bin(M[i] & v).count("1") & 1:
            y |= 1 << i
    return y

def columnas_activas(r, w, fraccion):
    if fraccion >= 1: return set(range(w))
    paso = int(round(1/fraccion))
    return {z for z in range(w) if (z+r) % paso == 0}

def dim_nucleo(R, w, rots, fraccion):
    """Dimension del subespacio de diferencias que esquivan todas las cajas-S."""
    M, n = matriz_lambda(w, rots)
    idx = lambda i,z: i*w + z
    # ecuaciones: para cada ronda r y cada (i,z) con z activa: fila de Lambda^(r+1)
    filas = []
    pot = [1<<j for j in range(n)]          # Lambda^0 aplicado a cada base = identidad
    # construimos Lambda^(r+1) por columnas: col j = Lambda^(r+1) e_j
    cols = [1<<j for j in range(n)]
    for r in range(R):
        cols = [aplicar(M, c, n) for c in cols]      # ahora cols = Lambda^(r+1) e_j
        act = columnas_activas(r, w, fraccion)
        for i in range(L):
            for z in act:
                fila = 0
                p = idx(i,z)
                for j in range(n):
                    if (cols[j] >> p) & 1:
                        fila |= 1 << j
                filas.append(fila)
    # rango sobre GF(2)
    rango = 0
    for c in range(n):
        piv = next((k for k in range(rango, len(filas)) if (filas[k]>>c)&1), None)
        if piv is None: continue
        filas[rango], filas[piv] = filas[piv], filas[rango]
        for k in range(len(filas)):
            if k != rango and (filas[k]>>c)&1:
                filas[k] ^= filas[rango]
        rango += 1
    return n - rango, n

if __name__ == "__main__":
    w = 40
    rots = [(19,28),(21,39),(1,6),(10,17),(7,31)]
    print("Estado: 5 x %d = %d bits\n" % (w, 5*w))
    print(" chi         | R | dim. del nucleo | trayectorias de probabilidad 1")
    print("-------------+---+-----------------+-------------------------------")
    for frac, nom in ((1.0,"completo"), (0.5,"parcial 1/2"), (0.25,"parcial 1/4")):
        for R in (1,2,3):
            d, n = dim_nucleo(R, w, rots, frac)
            print(" %-11s | %d | %3d de %3d      | %s"
                  % (nom, R, d, n, "NO" if d == 0 else "SI: 2^%d diferencias" % d))

def matriz_lambda_paridad(w, rots, cs):
    """CORRECCION: P = XOR de los 5 carriles;  L_i <- L_i ^ rot(P,c_i) ^ rot(L_i,a_i).
    Introduce mezcla ENTRE carriles en la capa lineal, no solo dentro de cada uno."""
    n = 5*w
    idx = lambda i,z: i*w + z
    M = [0]*n
    for i in range(L):
        a,_ = rots[i]; c = cs[i]
        for z in range(w):
            fila = (1<<idx(i,z)) ^ (1<<idx(i,(z-a)%w))
            for k in range(L):
                fila ^= 1 << idx(k,(z-c)%w)
            M[idx(i,z)] = fila
    return M, n

def dim_nucleo_gen(R, w, fraccion, M, n):
    idx = lambda i,z: i*w + z
    filas = []
    cols = [1<<j for j in range(n)]
    for r in range(R):
        cols = [aplicar(M, c, n) for c in cols]
        act = columnas_activas(r, w, fraccion)
        for i in range(L):
            for z in act:
                fila = 0; p = idx(i,z)
                for j in range(n):
                    if (cols[j]>>p)&1: fila |= 1<<j
                filas.append(fila)
    rango = 0
    for c in range(n):
        piv = next((k for k in range(rango,len(filas)) if (filas[k]>>c)&1), None)
        if piv is None: continue
        filas[rango],filas[piv] = filas[piv],filas[rango]
        for k in range(len(filas)):
            if k!=rango and (filas[k]>>c)&1: filas[k] ^= filas[rango]
        rango += 1
    return n - rango
