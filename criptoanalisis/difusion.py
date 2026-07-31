"""
Medida de difusion por matriz de dependencia sobre GF(2).

D_r[i][j] = 1 si el bit de salida i depende del bit de entrada j tras r rondas.
Difusion = densidad de D_r = (numero de unos) / n^2.
Se compone ronda a ronda: D_{r+1} = M o D_r (producto booleano OR-AND).

Cada paso aporta sus dependencias exactas:
  theta : 11 entradas por bit    rho, pi : permutacion (1)
  chi   : 3 entradas por bit     capa lineal por carril (XOR+rot): 3
"""
ROT = [[0,36,3,41,18],[1,44,10,45,2],[62,6,43,15,61],[28,55,25,21,56],[27,20,39,8,14]]
L = 5

def mult(A, B, n):
    """Producto booleano de matrices dadas como listas de bitmasks por fila."""
    return [ (lambda m: __import__('functools').reduce(lambda a,b: a|b,
             (B[j] for j in range(n) if (m>>j)&1), 0))(A[i]) for i in range(n) ]

def densidad(D, n):
    return sum(bin(r).count("1") for r in D) / (n*n)

# ---------------- Keccak-f[25w] ----------------
def keccak_ronda(w):
    n = 25*w
    idx = lambda x,y,z: (x*L + y)*w + z
    M = [0]*n
    for x in range(L):
        for y in range(L):
            for z in range(w):
                # chi: depende de 3 bits tras rho,pi de theta
                dep = 0
                for dx in (0,1,2):
                    xx = (x+dx) % L
                    # deshacer pi: (xx,y) viene de (a,b) con b=xx, a=(xx+3y)%5... 
                    a = (xx + 3*y) % L; b = xx
                    sh = ROT[a][b] % w
                    zz = (z - sh) % w
                    # theta sobre (a,b,zz): 11 bits
                    dep |= 1 << idx(a,b,zz)
                    for yy in range(L):
                        dep |= 1 << idx((a-1)%L, yy, zz)
                        dep |= 1 << idx((a+1)%L, yy, (zz-1)%w)
                M[idx(x,y,z)] = dep
    return M, n

# ---------------- Variante: 5 carriles de w bits ----------------
def variante_ronda(w, rots):
    """rots[i] = (a_i, b_i): L_i <- L_i ^ rot(L_i,a_i) ^ rot(L_i,b_i); luego chi por columna."""
    n = 5*w
    idx = lambda i,z: i*w + z
    M = [0]*n
    for i in range(5):
        for z in range(w):
            dep = 0
            for di in (0,1,2):          # chi: 3 carriles
                j = (i+di) % 5
                a,b = rots[j]
                for sh in (0, a, b):    # capa lineal: 3 posiciones
                    dep |= 1 << idx(j, (z - sh) % w)
            M[idx(i,z)] = dep
    return M, n

def rondas_hasta(M, n, umbral=0.80, maxr=20):
    D = [1 << i for i in range(n)]      # identidad
    hist = []
    for r in range(1, maxr+1):
        D = mult(M, D, n)
        d = densidad(D, n)
        hist.append((r, d))
        if d >= umbral:
            return r, hist
    return None, hist

if __name__ == "__main__":
    print("=== LINEA BASE: Keccak-f[200] (5x5x8) ===")
    M, n = keccak_ronda(8)
    r, h = rondas_hasta(M, n)
    for rr, d in h[:8]:
        print("  ronda %2d: difusion %6.2f%%" % (rr, 100*d))
    print("  -> 80%% alcanzado en la ronda %s\n" % r)

# ---------------- Variante con paridad de carriles (theta simplificado) ----------------
def variante_paridad(w, rots, cs):
    """P = XOR de los 5 carriles;  L_i <- L_i ^ rot(P,c_i) ^ rot(L_i,a_i);  luego chi."""
    n = 5*w
    idx = lambda i,z: i*w + z
    M = [0]*n
    for i in range(5):
        for z in range(w):
            dep = 0
            for di in (0,1,2):                 # chi
                j = (i+di) % 5
                a,_ = rots[j]; c = cs[j]
                dep |= 1 << idx(j, z)                      # L_j
                dep |= 1 << idx(j, (z-a) % w)              # rot(L_j,a_j)
                for k in range(5):                         # rot(P,c_j): los 5 carriles
                    dep |= 1 << idx(k, (z-c) % w)
            M[idx(i,z)] = dep
    return M, n
