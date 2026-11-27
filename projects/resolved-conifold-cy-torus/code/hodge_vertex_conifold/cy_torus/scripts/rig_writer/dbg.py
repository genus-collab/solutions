from fractions import Fraction as Fr
import random, math
exec(open('checks.py').read().split('# ---------- C. critical')[0].split('# ---------- B.')[0])  # helpers only
random.seed(1)
def rnd():
    while True:
        v = Fr(random.randint(-30, 30), random.randint(1, 12))
        if v != 0: return v
def padd(a, b):
    n = max(len(a), len(b)); r = [Fr(0)] * n
    for i, v in enumerate(a): r[i] += v
    for i, v in enumerate(b): r[i] += v
    return r
def pmul(a, b):
    r = [Fr(0)] * (len(a) + len(b) - 1)
    for i, v in enumerate(a):
        for j, w in enumerate(b): r[i + j] += v * w
    return r
def pscale(a, c): return [c * v for v in a]
def pmod(a, m):
    a = a[:]
    while len(a) >= len(m):
        if a[-1] == 0: a.pop(); continue
        c = a[-1] / m[-1]; d = len(a) - len(m)
        for i, v in enumerate(m): a[d + i] -= c * v
        a.pop()
    return a
def iszero(a): return all(v == 0 for v in a)
def Disc(Q, a0, a1, a2): return a0 ** 2 + (2 * a0 * (a1 + a2) + 4 * a1 * a2) * Q + (a1 - a2) ** 2 * Q ** 2
res = {}
for _ in range(30):
    a0, a1, a2, Q, e1 = rnd(), rnd(), rnd(), rnd(), rnd()
    if Q == 1 or a1 == a2: continue
    E, A, B = a0 + a1 + a2, a0 + a2, a0 + a1
    e2 = E - e1; tau = 1 - 1 / Q
    P = [-a2, -(a0 + Q * (a1 - a2)), B * Q]
    one_U = [Fr(1), Fr(-1)]; one_QU = [Fr(1), -Q]
    expr = padd(padd(pscale(pmul(one_U, one_QU), -a2), pscale(pmul([Fr(0), Fr(1)], one_QU), -A)), pscale(pmul([Fr(0), Fr(1)], one_U), -a1 * Q))
    res['c1'] = res.get('c1', True) and iszero(padd(expr, pscale(P, -1)))
    res['c2'] = res.get('c2', True) and ((a0 + Q * (a1 - a2)) ** 2 + 4 * B * Q * a2 == Disc(Q, a0, a1, a2))
    L_num = [Fr(0), -a1 * Q]
    num = padd(padd(pscale(pmul(L_num, L_num), tau), pscale(pmul(L_num, one_QU), -(E - a0 * tau))), pscale(pmul(one_QU, one_QU), a1 * a2))
    res['c3'] = res.get('c3', True) and iszero(pmod(num, P))
    twoL = padd(pscale(L_num, 2), pscale(one_QU, -(E / tau - a0)))
    a1a2EL = padd(pscale(one_QU, a1 * a2), pscale(L_num, -E))
    lhs = pscale(pmul(pmul(twoL, a1a2EL), pmul(one_U, one_U)), e1 * e2)
    W2num = padd(pscale(pmul([Fr(0), Fr(1)], pmul(one_QU, one_QU)), -A), pscale(pmul([Fr(0), Fr(1)], pmul(one_U, one_U)), -a1 * Q))
    rhs = pscale(W2num, a1 * A * e1 * e2)
    res['c4'] = res.get('c4', True) and iszero(pmod(padd(lhs, pscale(rhs, -1)), P))
    alt = padd(padd(pscale(pmul(one_QU, one_QU), a0 * a1 * a2), pscale(pmul(L_num, one_QU), 2 * a1 * a2)), pscale(pmul(L_num, L_num), -E))
    alt = pscale(pmul(alt, pmul(one_U, one_U)), e1 * e2)
    res['c5'] = res.get('c5', True) and iszero(pmod(padd(lhs, pscale(alt, -1)), P))
    res['c6'] = res.get('c6', True) and (sum(c * Fr(1) ** i for i, c in enumerate(P)) == A * (Q - 1))
    res['c7'] = res.get('c7', True) and (sum(c * (1 / Q) ** i for i, c in enumerate(P)) == B * (1 - Q) / Q)
    P0 = [-a2, -(Q * (a1 - a2)), a1 * Q]; xx = a1 / a2
    res['c8'] = res.get('c8', True) and (P0 == pscale([Fr(-1), Q * (1 - xx), xx * Q], a2))
print(res)
# D debug
Nn = 7
def ser_add(a, b): return [x + y for x, y in zip(a, b)]
def ser_scal(a, c): return [c * x for x in a]
a1, a2 = Fr(3), Fr(9,5)
a0s = [Fr(0), Fr(1)] + [Fr(0)] * (Nn - 1); one = [Fr(1)] + [Fr(0)] * Nn
def DiscS(Qs):
    c1 = ser_add(ser_scal(a0s, 2 * (a1 + a2)), ser_scal(one, 4 * a1 * a2))
    return ser_add(ser_add(smul(a0s, a0s, Nn), smul(c1, Qs, Nn)), ser_scal(smul(Qs, Qs, Nn), (a1 - a2) ** 2))
def DiscpS(Qs):
    c1 = ser_add(ser_scal(a0s, 2 * (a1 + a2)), ser_scal(one, 4 * a1 * a2))
    return ser_add(c1, ser_scal(Qs, 2 * (a1 - a2) ** 2))
Qs = [Fr(0)] * (Nn + 1)
for it in range(8):
    Qs = ser_add(Qs, ser_scal(smul(DiscS(Qs), sinv(DiscpS(Qs), Nn), Nn), -1))
print("Disc(Q1)=", [str(v) for v in DiscS(Qs)])
print("Q1=", [str(v) for v in Qs])
print("pred c2,c3,c4:", -1/(4*a1*a2), (a1+a2)/(8*a1**2*a2**2), -(a1**2+3*a1*a2+a2**2)/(16*a1**3*a2**3))
# exact quadratic formula check: Q1 = (-(2a0(a1+a2)+4a1a2) + sqrt(16 a1 a2 (a0+a1)(a0+a2)))/(2 (a1-a2)^2) -> series of sqrt
# E debug
def Li_poly(n):
    Apoly = [Fr(0), Fr(1)]
    for m in range(n):
        d = [i * Apoly[i] for i in range(1, len(Apoly))]
        t1 = pmul([Fr(1), Fr(-1)], d) if d else [Fr(0)]
        t2 = pscale(Apoly, m + 1)
        Apoly = pmul([Fr(0), Fr(1)], padd(t1, t2))
    return Apoly
def Li_rat(n, Q):
    Apoly = Li_poly(n)
    return sum(c * Q ** i for i, c in enumerate(Apoly)) / (1 - Q) ** (n + 1)
for n in range(1, 5):
    Ap = Li_poly(n); print(n, [str(v) for v in Ap], Li_rat(n, Fr(1,2)), float(sum(Fr(d)**n * Fr(1,2)**d for d in range(1,120))))
    Q = Fr(3,7); print("  inversion:", Li_rat(n, 1/Q) == (-1)**(n+1)*Li_rat(n, Q))
