# Referee U05_B: independent exact test of Lemma HPF (lem:gau-HPF), parts (i), (ii).
# Ring: Laurent polynomials in s (s^5 = u), polynomials in z.  D s = -(1-s^5)/(5 s^4), D z = 0.
from fractions import Fraction as Fr
import random

random.seed(7)


def add(a, b, cb=1):
    d = dict(a)
    for k, x in b.items():
        d[k] = d.get(k, 0) + cb * x
        if d[k] == 0:
            del d[k]
    return d


def mul(a, b):
    d = {}
    for (i1, j1), x1 in a.items():
        for (i2, j2), x2 in b.items():
            k = (i1 + i2, j1 + j2)
            d[k] = d.get(k, 0) + x1 * x2
    return {k: x for k, x in d.items() if x != 0}


def D(a):  # q d/dq on s-exponent c: D s^c = -(c/5)(s^{c-5} - s^c)
    d = {}
    for (i, j), x in a.items():
        if i != 0:
            for k, y in (((i - 5, j), -Fr(i, 5) * x), ((i, j), Fr(i, 5) * x)):
                d[k] = d.get(k, 0) + y
    return {k: x for k, x in d.items() if x != 0}


z = {(0, 1): Fr(1)}
L = {(-1, 0): Fr(1)}
one = {(0, 0): Fr(1)}
u = {(5, 0): Fr(1)}
q = {(0, 0): Fr(1, 5**5), (5, 0): -Fr(1, 5**5)}  # q = (1-u)/5^5


def M0(a):  # z D + L
    return add(mul(z, D(a)), mul(L, a))


def PF5_0(a):  # (zD+L)^5 - 1 - q prod_{k=1}^5 (5(zD+L) + k z), from eq. (set-PF5)
    t = a
    for _ in range(5):
        t = M0(t)
    r = add(t, a, -1)
    w = a
    for k in range(5, 0, -1):  # apply (5M + k z) for k = 5,...,1 (they commute; z is a constant)
        w = add({kk: 5 * x for kk, x in M0(w).items()}, mul({(0, 1): Fr(k)}, w))
    return add(r, mul(q, w), -1)


# Hatted side: series in (s, zh): dict keyed (s-exp, zh-exp); theta_u s^i = (i/5) s^i
def O_e(e, F):
    # 1 + zh (1-u) (4/5 theta_zh - e/5 - theta_u)
    G = {}
    for (i, j), x in F.items():
        c = Fr(4 * j, 5) - Fr(e, 5) - Fr(i, 5)
        if c != 0:
            for k, y in (((i, j + 1), c * x), ((i + 5, j + 1), -c * x)):
                G[k] = G.get(k, 0) + y
    return add(F, {k: x for k, x in G.items() if x != 0})


def hat_to_flat(F):  # zh^j -> z^j s^{-4j}
    d = {}
    for (i, j), x in F.items():
        d[(i - 4 * j, j)] = d.get((i - 4 * j, j), 0) + x
    return d


sig = [None, Fr(3), Fr(17, 5), Fr(9, 5), Fr(274, 625), Fr(24, 625)]


def HPF_lhs(Phi):
    Ph = [Phi]
    for j in range(1, 6):
        Ph.append(O_e(-j, Ph[-1]))
    res = add(Ph[5], Phi, -1)
    one_minus_u = {(0, 0): Fr(1), (5, 0): Fr(-1)}
    br = {}
    for k in range(1, 6):
        term = mul({(5 * (k - 1), k): sig[k]}, Ph[5 - k])
        br = add(br, term)
    return add(res, mul(one_minus_u, br), -1), Ph


def rand_ph(nmax, frac_powers):
    Phi = {}
    for n in range(nmax + 1):
        exps = range(-3, 9) if frac_powers else range(0, 30, 5)
        for i in exps:
            if random.random() < 0.6:
                Phi[(i, n)] = Fr(random.randint(-9, 9), random.randint(1, 7))
    return {k: x for k, x in Phi.items() if x != 0}


ok = True
for frac in (False, True):
    for trial in range(3):
        Phi = rand_ph(4, frac)
        lhs, Ph = HPF_lhs(Phi)
        # (i): PF5^{(0)}(L Phi) = s^{-1} * lhs
        LPhi = mul(L, hat_to_flat(Phi))
        a = PF5_0(LPhi)
        b = mul(L, hat_to_flat(lhs))
        if add(a, b, -1):
            ok = False
            print("MISMATCH (i)", frac, trial)
        # (ii): coefficient of zh^{n+1} of lhs = (1-u)[(4n - 5 theta_u) ph_n + r_n]
        for n in range(0, 5):
            Phi_flat0 = {k: x for k, x in Phi.items() if k[1] != n}  # ph_n replaced by 0
            _, Phb = HPF_lhs(Phi_flat0)

            def cz(F, m):
                return {(i, 0): x for (i, j), x in F.items() if j == m}

            def lin(c, F):  # (c - theta_u) F for F in s only
                return {(i, 0): (c - Fr(i, 5)) * x for (i, j), x in F.items() if (c - Fr(i, 5)) * x != 0}

            rn = {}
            for j in range(1, 6):
                rn = add(rn, lin(Fr(4 * n + j, 5), cz(Phb[j - 1], n)))
            rn = add(rn, cz(Phb[4], n), -3)
            for k in range(2, 6):
                rn = add(rn, mul({(5 * (k - 1), 0): sig[k]}, cz(Ph[5 - k], n - (k - 1))), -1)
            phn = cz(Phi, n)
            inner = add(lin(Fr(4 * n, 5), phn), {})
            inner = {k: 5 * x for k, x in inner.items()}  # (4n - 5 theta_u) ph_n
            inner = add(inner, rn)
            pred = mul({(0, 0): Fr(1), (5, 0): Fr(-1)}, inner)
            if add(pred, cz(lhs, n + 1), -1):
                ok = False
                print("MISMATCH (ii)", frac, trial, n)
print("HPF (i),(ii) exact checks passed:", ok)

# r_1 = 3/5 when ph_0 = 1 (ph_1 arbitrary), and ph_1 = -3/20
Phi = {(0, 0): Fr(1), (0, 1): Fr(5, 3), (5, 1): Fr(2, 7)}
lhs, Ph = HPF_lhs(Phi)
c2 = {i: x for (i, j), x in lhs.items() if j == 2}
# predicted: (1-u)[(4 - 5 theta_u) ph_1 + 3/5]
ph1 = {0: Fr(5, 3), 5: Fr(2, 7)}
inner = {i: (4 - i) * x for i, x in ph1.items()}
inner[0] = inner.get(0, 0) + Fr(3, 5)
pred = {}
for i, x in inner.items():
    pred[i] = pred.get(i, 0) + x
    pred[i + 5] = pred.get(i + 5, 0) - x
pred = {i: x for i, x in pred.items() if x != 0}
print("r_1 = 3/5 check:", pred == c2)
