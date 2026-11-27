# Extra refutation tests (round 6d, backward read). Evidence only.
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import time
from itertools import product
from fractions import Fraction as Fr
from gs import *
from run_full import enumerate_graphs

e1, e2 = Fr(3), Fr(7)
E = e1 + e2
print('--- (a) values at t_c = 0 against draft s.6 table, (3,7)')
R = Rows(e1, e2, 12, 20)
print(' XX(0) =', ratrec(R.b.XX.coef(0)), ' tc =', [ratrec(R.b.tc.coef(k)) for k in range(1, 4)])
print(' E00(0,0) at 0 =', ratrec(R.edge(0, 0, 0, 0).coef(0)), ' expected', -e1 * e2 * E ** 2 / 27)
print(' val E00(a,b):', sorted(set(R.edge(0, 0, a, b).val for a in range(5) for b in range(5))))
print(' val E01(a,b):', sorted(set(R.edge(0, 1, a, b).val for a in range(5) for b in range(5))), ' E12:', sorted(set(R.edge(1, 2, a, b).val for a in range(5) for b in range(5))))
print(' Dil0(d) at 0 vs -E^(1-d):', all(R.Dil(0, d).coef(0) == fp(-E ** (1 - d)) for d in range(2, 9)))
om = OMEGA
for j in (1, 2):
    w = pow(om, j, P)
    val = fp(e1 * e2 * E ** 3) * pow(w, 2, P) % P * inv(3 * (1 - w) ** 2 % P) % P
    print(' Delta_%d(0) ok:' % j, R.colour_data[j]['Delta'].coef(0) == val)
# symplecticity at free F, sample
print(' symplectic N(-w,w) check (colour pairs, n<=6):')
ok = True
for i in range(3):
    for j in range(3):
        for n in range(1, 7):
            s = None
            for M in range(0, n + 1):
                b = n - M
                x = R.pair(i, j, M, b)
                x = x if b % 2 == 0 else -x
                s = x if s is None else s + x
            if s.c:
                ok = False
print('  all z^n coefficients of <P(z),P(-z)> vanish for n>=1:', ok)

print('--- (b) sensitivity: G_g[F_cf + c xi^k] - G_g[F_cf]: valuation in xi. Shows which Taylor coefficients of F_cf the polar part sees')
for g in (2, 3, 4):
    gr = enumerate_graphs(g)
    N = 2 * g + 8
    zmax = 3 * g + 3
    base = Base(N)

    def Gsum(FF):
        Rr = Rows(e1, e2, N, zmax, base=base, FF=FF)
        Ev = Evaluator(Rr, Vert(e1, e2))
        tot = None
        for gen, ed, aut in gr:
            if len(ed) == 0:
                continue
            for c in product(range(3), repeat=len(gen)):
                T = Ev.graph(gen, ed, c) * inv(aut)
                tot = T if tot is None else tot + T
        return tot
    if g == 4:
        break  # g=4 handled below with fewer perturbations
    G0 = Gsum(base.XX)
    for k in range(0, 6):
        pert = LS(k, [fp(Fr(1, 9))], N + 3)
        G1 = Gsum(base.XX + pert)
        d = G1 - G0
        print(' g', g, 'k', k, ' val(diff) =', d.val)

print('--- (c) mode A control, w -> -w (leaf scale not a weight): does an edge graph acquire a pole?')
for g in (2, 3):
    gr = enumerate_graphs(g)
    Rr = Rows(e1, e2, 2 * g + 8, 3 * g + 3)
    Ev = Evaluator(Rr, Vert(e1, e2, w=(-E, e1, e2)))
    nbad = 0
    nall = 0
    for gen, ed, aut in gr:
        if not ed:
            continue
        for c in product(range(3), repeat=len(gen)):
            T = Ev.graph(gen, ed, c)
            nall += 1
            if T.c and T.val < 0:
                nbad += 1
    print(' g', g, 'edge-graph colourings', nall, 'with a pole under mode A', nbad)

print('--- (d) individual decorated terms versus the sum, graph {2}-{1} at genus 3, and {2}-{2} at genus 4')


def decorated(genera, edges, col, Rr, Ev):
    nv = len(genera)
    hes = [[] for _ in range(nv)]
    for k, (u, v) in enumerate(edges):
        hes[u].append((k, 0))
        hes[v].append((k, 1))
    maxa = [3 * genera[i] - 3 + len(hes[i]) for i in range(nv)]
    out = []
    ne = len(edges)

    def rec(k, cur, assign, used):
        if k == ne:
            res = cur
            for i in range(nv):
                vv = Ev.vertex(genera[i], col[i], tuple(assign[h] for h in hes[i]))
                if not vv.c and vv.aprec >= 10 ** 6:
                    return
                res = res * vv
            out.append((dict(assign), res))
            return
        u, v = edges[k]
        for a in range(0, maxa[u] - used[u] + 1):
            for b in range(0, maxa[v] - used[v] + 1):
                if u == v and a + b > maxa[u] - used[u]:
                    continue
                assign[(k, 0)] = a
                assign[(k, 1)] = b
                used[u] += a
                used[v] += b
                rec(k + 1, cur * Rr.edge(col[u], col[v], a, b), assign, used)
                used[u] -= a
                used[v] -= b
    rec(0, LS(0, [1], 10 ** 6), {}, [0] * nv)
    return out


for genera, edges, cols in [((2, 1), [(0, 1)], [(0, 0), (0, 1), (0, 2)]), ((2, 2), [(0, 1)], [(0, 0), (0, 1)]), ((3, 1), [(0, 1)], [(0, 0), (0, 1)]), ((3,), [(0, 0)], [(0,)])]:
    g = sum(genera) + len(edges) - len(genera) + 1
    Rr = Rows(e1, e2, 2 * g + 10, 3 * g + 6)
    Ev = Evaluator(Rr, Vert(e1, e2))
    for col in cols:
        terms = decorated(genera, edges, col, Rr, Ev)
        vals = sorted(t[1].val for t in terms if t[1].c)
        s = None
        for _, x in terms:
            s = x if s is None else s + x
        print(' graph', genera, edges, 'colours', col, 'decorations', len(terms), 'min val of a decorated term', vals[0] if vals else None,
              ' val of the sum', s.val if s.c else None)

print('--- (e) leaf membership, colour 0 half-edge against colour 1 and colour 0: numerator degree test of Lambda_{k+1}')
# rho = -E v_0 / t_c ; [t_c^k][w^b] rho E_0j(psi,w) should be Sum_{l<=k+1} c_l lambda_l(psi), no constant.
ZM = 34
Rr = Rows(e1, e2, 10, ZM)
b_ = Rr.b
rho = Rr.colour_data[0]['v'] * b_.tc.inverse() * (P - Rr.E_)
tc = b_.tc
# convert a series in xi to coefficients in t_c: compute by successive subtraction


def tc_coefs(f, Kmax):
    out = {}
    cur = f
    v = cur.val
    tcinvpow = None
    k = min(0, cur.val)
    res = []
    while k <= Kmax:
        c = cur.coef(k)
        res.append((k, c))
        if c:
            cur = cur - tc.pow(k) * c
        k += 1
    return dict(res)


def numerator_degree(series, K):
    # series: list of F_p coefficients in psi. multiply by prod_{j=1}^{K}(E - j psi); return degree of the result among known coefficients
    s = list(series)
    for j in range(1, K + 1):
        ns = [0] * len(s)
        for i in range(len(s)):
            ns[i] = (Rr.E_ * s[i] - (j * s[i - 1] if i else 0)) % P
        s = ns
    d = -1
    for i, x in enumerate(s):
        if x:
            d = i
    return d, s


for (i, j) in [(0, 1), (0, 2), (0, 0)]:
    for b in (0, 1, 3):
        for k in range(0, 5):
            ser = []
            for a in range(0, 24):
                x = Rr.edge(i, j, a, b) * rho
                ser.append(x.coef(0) if k == 0 else None)
            # general k through tc_coefs
            ser = []
            for a in range(0, 24):
                x = Rr.edge(i, j, a, b) * rho
                ser.append(tc_coefs(x, k)[k])
            d, s = numerator_degree(ser, k + 1)
            # membership in Lambda_{k+1}, no constant: numerator degree <= k
            print(' edge', (i, j), 'b', b, 'k', k, 'numerator degree', d, 'expected <=', k, 'OK' if d <= k else 'FAIL')
print('done')
