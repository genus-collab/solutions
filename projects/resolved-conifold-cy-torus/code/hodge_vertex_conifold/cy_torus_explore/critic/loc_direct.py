# Critic: direct virtual localisation of GW_{g,d}(Z) on the CY 4-torus, NO R-matrix, NO graph sum of the CohFT.
# Compares with P_g d^(2g-3)/(e1 e2). Hodge integrals from row_audit/hodge_engine via Omega. Exact in F_p.
# Usage: py -I loc_direct.py gmax dmax a0 a1 a2 e1
import sys, os, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'row_audit'))
from fractions import Fraction as Fr
from math import factorial
from cylib import P, inv, fp, BERN
import hodge_engine as he
from gsum import Omega

def inv_(x): return inv(x % P)

class Loc:
    def __init__(s, a0, a1, a2, e1, e2):
        s.a0, s.a1, s.a2, s.e1, s.e2 = [fp(x) for x in (a0, a1, a2, e1, e2)]
        A = (s.a0 + s.a2) % P; B = (s.a0 + s.a1) % P
        s.wt = {+1: [s.a0, s.a1, s.a2, -s.e1 % P, -s.e2 % P], -1: [-s.a0 % P, B, A, -s.e1 % P, -s.e2 % P]}
        s.om = {c: Omega(s.wt[c]) for c in s.wt}
        s.pw = {}
        for c in s.wt:
            r = 1
            for x in s.wt[c]: r = r * x % P
            s.pw[c] = r
        s.vc = {}

    def vertex(s, c, g, ws):
        ws = tuple(sorted(ws))
        k = (c, g, ws)
        if k in s.vc: return s.vc[k]
        n = len(ws)
        if g == 0 and n == 1:
            r = ws[0]
        elif g == 0 and n == 2:
            r = s.pw[c] * inv_(ws[0] + ws[1]) % P
        else:
            dim = 3 * g - 3 + n
            tot = 0
            for ds in itertools.product(range(dim + 1), repeat=n):
                if sum(ds) > dim: continue
                if g == 0 and sum(ds) != dim: continue
                if g == 0:
                    val = fp(he.wk(0, tuple(sorted(ds))))
                else:
                    val = s.om[c].integral(g, ds)
                if not val: continue
                t = val
                for d, w in zip(ds, ws):
                    t = t * inv_(pow(w, d + 1, P)) % P
                tot += t
            r = tot % P * pow(s.pw[c], n - 1 + g, P) % P
        s.vc[k] = r
        return r

    def edge(s, d):
        id_ = inv_(d)
        num = 1
        for j in range(1, d):
            num = num * (s.a1 + j * s.a0 * id_) % P * (s.a2 + j * s.a0 * id_) % P
        den = 1
        for j in range(1, d + 1):
            x = j * s.a0 * id_ % P
            den = den * x % P * (-x) % P
        den = den * s.e1 % P * s.e2 % P
        return id_ * num % P * inv_(den) % P

    def gw(s, g, d):
        tot = 0
        for m in range(1, d + 1):
            for comp in compositions(d, m):
                for vp in range(1, m + 1):
                    for vm in range(1, m + 1):
                        b1 = m - vp - vm + 1
                        if b1 < 0 or b1 > g: continue
                        G = g - b1
                        pref = inv_(factorial(vp) * factorial(vm) * factorial(m))
                        for assign in itertools.product(itertools.product(range(vp), range(vm)), repeat=m):
                            if not connected_all(assign, vp, vm): continue
                            for gs in compositions0(G, vp + vm):
                                t = pref
                                for e, dd in enumerate(comp):
                                    t = t * s.edge(dd) % P
                                ok = True
                                for i in range(vp):
                                    ws = [s.a0 * inv_(comp[e]) % P for e in range(m) if assign[e][0] == i]
                                    t = t * s.vertex(+1, gs[i], ws) % P
                                for j in range(vm):
                                    ws = [(-s.a0) * inv_(comp[e]) % P for e in range(m) if assign[e][1] == j]
                                    t = t * s.vertex(-1, gs[vp + j], ws) % P
                                tot += t
        return tot % P

def compositions(d, m):
    if m == 1:
        yield (d,); return
    for x in range(1, d - m + 2):
        for r in compositions(d - x, m - 1):
            yield (x,) + r

def compositions0(g, k):
    if k == 1:
        yield (g,); return
    for x in range(g + 1):
        for r in compositions0(g - x, k - 1):
            yield (x,) + r

def connected_all(assign, vp, vm):
    used_p = {a for a, b in assign}; used_m = {b for a, b in assign}
    if len(used_p) != vp or len(used_m) != vm: return False
    par = list(range(vp + vm))
    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for a, b in assign:
        par[f(a)] = f(vp + b)
    return len({f(x) for x in range(vp + vm)}) == 1

def Pg(e1, e2, g):
    # [t^(2g)] e1 e2 t^2/(4 sinh(e1 t/2) sinh(e2 t/2))
    def c(i):
        return Fr(2 - 2 ** (2 * i), 2 ** (2 * i)) * BERN[2 * i] / factorial(2 * i)
    r = 0
    for i in range(g + 1):
        r += fp(c(i)) * pow(e1, 2 * i, P) % P * fp(c(g - i)) % P * pow(e2, 2 * (g - i), P)
    return r % P

if __name__ == '__main__':
    gmax, dmax = int(sys.argv[1]), int(sys.argv[2])
    a0, a1, a2, e1 = [Fr(x) for x in sys.argv[3:7]]
    e2 = a0 + a1 + a2 - e1 if len(sys.argv) < 8 else Fr(sys.argv[7])
    L = Loc(a0, a1, a2, e1, e2)
    print('weights', a0, a1, a2, e1, e2, 'CY' if a0 + a1 + a2 == e1 + e2 else 'nonCY', flush=True)
    for g in range(0, gmax + 1):
        for d in range(1, dmax + 1):
            v = L.gw(g, d)
            if g == 0:
                ex = inv_(pow(d, 3, P) * L.e1 % P * L.e2 % P)
            elif g == 1:
                ex = Pg(L.e1, L.e2, 1) * inv_(d) % P * inv_(L.e1 * L.e2 % P) % P
            else:
                ex = Pg(L.e1, L.e2, g) * pow(d, 2 * g - 3, P) % P * inv_(L.e1 * L.e2 % P) % P
            print('g', g, 'd', d, 'localisation == closed form:', v == ex, flush=True)
