"""Stage-2 critic, core of cr2_legdict.py: exact legs on x = p^4/4 + v/p, y = p, the closed form of PiC, and residue
pairings at the five branch points modulo a prime P with 5 | P - 1 (2^61-1 and 2^89-1 both qualify).

Legs (dict worker's convention): Phi_0^b = p, log p, v/p, p^3, p^2 (b = 0..4); Phi_(k+1)^b = -Int Phi_k^b dx, termwise in p,
integration constant 0.  A leg is a dict (j, l, m) -> Fraction meaning coeff * v^m p^j (log p)^l, l in {0, 1}.
Functions xif_a of arXiv:1509.06954 (2-34) (TR_MODEL s.2.2 at v = 1; the v-powers alpha are forced by homogeneity):
  xif_a = eps_a v^alpha_a p^m_a / (p^5 - v), (m, eps, alpha) = (4,-1,0), (0,1,1), (1,1,0), (2,-1,0), (3,-1,0);
  xif_{a,d} = D^d xif_a, D = d/dx = p^2/(p^5 - v) d/dp.
Pairing <Phi, w> = Sum_j Res_{c_j} Phi w; at each c only the principal part of D^d xif_a is needed:
  Res_c(phi d F) = -Sum_k k phi_k [t^-k] F  (t = p - c; exact identity for Laurent series).
log p near c is expanded as log c + log(1 + t/c); the constant log c is dropped, which is legitimate because
Res_c f(x) d xif_{a,d} = 0 at each c for every polynomial f (checked as L2 in cr2_legdict.py).
Imported read-only by cr2_legdict.py; stdlib only.
"""
from fractions import Fraction as Fr
from math import factorial

CB = (-1, 1, -1, -3, -2)
XI = ((4, -1, 0), (0, 1, 1), (1, 1, 0), (2, -1, 0), (3, -1, 0))
FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg)
        print("FAIL", msg)
    return cond


def leg0(b, ctrl=None):
    if ctrl == "C2" and b == 2:
        return {(4, 0, 0): Fr(1, 4), (-1, 0, 1): Fr(1)}          # the function x itself, instead of v/p
    return [{(1, 0, 0): Fr(1)}, {(0, 1, 0): Fr(1)}, {(-1, 0, 1): Fr(1)}, {(3, 0, 0): Fr(1)}, {(2, 0, 0): Fr(1)}][b]


def step(F, sign=-1):
    """sign * Int F dx, dx = (p^3 - v p^-2) dp, termwise; a (log p)^2 term must cancel."""
    integrand = {}
    for (j, l, m), c in F.items():
        for jj, mm, cc in ((j + 3, m, c), (j - 2, m + 1, -c)):
            integrand[(jj, l, mm)] = integrand.get((jj, l, mm), 0) + cc
    G = {}

    def add(key, c):
        G[key] = G.get(key, 0) + c
        if G[key] == 0:
            del G[key]
    for (n, l, m), c in integrand.items():
        if c == 0:
            continue
        if l == 0:
            add((0, 1, m) if n == -1 else (n + 1, 0, m), c if n == -1 else c / (n + 1))
        elif check(n != -1, "a (log p)^2 term would appear"):
            add((n + 1, 1, m), c / (n + 1))
            add((n + 1, 0, m), -c / (n + 1) ** 2)
    return {k: sign * c for k, c in G.items()}


def legs(kmax, ctrl=None):
    T = {}
    for b in range(5):
        cur = leg0(b, ctrl)
        for k in range(kmax + 1):
            T[(b, k)] = cur
            cur = step(cur, +1 if ctrl == "C1" else -1)
    return T


def split(F):
    check(all(l <= 1 for (_, l, _) in F), "log power > 1")
    return ({(j, m): c for (j, l, m), c in F.items() if l == 0}, {(j, m): c for (j, l, m), c in F.items() if l == 1})


def pic(b, n):
    """PiC(tau_n H^b) = (-1)^(n+1) [z^-n] Y^C_b, Y^C = (I_0, -I_-1, -I_3, I_2, I_1) (Lemma V): (coeff, v-power) or None."""
    mm, sg = [(0, 1), (-1, -1), (3, -1), (2, 1), (1, 1)][b]
    if n < 0 or (n - mm - 1) % 5:
        return None
    i = (n - mm - 1) // 5
    e = mm + 4 * i + 1
    if i < 0 or e < 0:
        return None
    return Fr((-1) ** i, 4 ** i * factorial(i)) * Fr((-1) ** e, factorial(e)) * sg * (-1) ** (n + 1), e


def as_poly_in_x(G):
    """G(p) (dict (j,m) -> c) = f(x(p)); returns f as dict (n,m) -> c, or None if G is not a polynomial in x."""
    G, f = dict(G), {}
    while G:
        j, m = max(G)
        c = G[(j, m)]
        if j < 0 or j % 4:
            return None
        n = j // 4
        f[(n, m)] = f.get((n, m), 0) + c * 4 ** n
        for t in range(n + 1):        # subtract c 4^n v^m x^n, x^n = Sum_t C(n,t) 4^(t-n) v^t p^(4n-5t)
            key = (4 * n - 5 * t, m + t)
            G[key] = G.get(key, 0) - c * 4 ** n * Fr(factorial(n), factorial(t) * factorial(n - t)) * Fr(4) ** (t - n)
            if G[key] == 0:
                del G[key]
    return f


def fr(c, P):
    return c.numerator % P * pow(c.denominator % P, P - 2, P) % P


def zeta5(P):
    for g in range(2, 10 ** 4):
        z = pow(g, (P - 1) // 5, P)
        if z != 1:
            return z


def mulser(a, b, n, P):
    out = [0] * n
    for i, x in enumerate(a[:n]):
        if x:
            for j, y in enumerate(b[:n - i]):
                out[i + j] = (out[i + j] + x * y) % P
    return out


def invser(a, n, P):
    inv0 = pow(a[0], P - 2, P)
    out = [inv0] + [0] * (n - 1)
    for k in range(1, n):
        s = sum(a[i] * out[k - i] for i in range(1, min(k, len(a) - 1) + 1)) % P
        out[k] = -s * inv0 % P
    return out


def binser(c, j, n, P):
    """(c + t)^j as a power series in t, n terms (j any integer)."""
    out, bk, ci, cj = [], 1, pow(c, P - 2, P), pow(c, j % (P - 1), P)
    for k in range(n):
        out.append(cj * bk % P)
        bk = bk * ((j - k) % P) % P * pow(k + 1, P - 2, P) % P * ci % P
    return out


def principal_parts(c, v0, P, dmax, ctrl=None):
    """PP[(a, d)][k] = [t^-k] D^d xif_a at the branch point c, k = 1..2d+1 (index 0 unused)."""
    L0 = 2 * dmax + 3
    h = [5 * pow(c, 4, P) % P, 10 * pow(c, 3, P) % P, 10 * c * c % P, 5 * c % P, 1]     # ((c+t)^5 - c^5)/t
    hinv = invser(h, L0, P)
    Q = mulser(binser(c, 2, 3, P), hinv, L0, P)                                        # D = t^-1 Q(t) d/dt
    PP = {}
    for a, (m, eps, al) in enumerate(XI):
        pref = eps * pow(v0, 0 if (ctrl == "C3" and a == 1) else al, P) % P
        low, ser = -1, [pref * x % P for x in mulser(binser(c, m, L0, P), hinv, L0, P)]
        for d in range(dmax + 1):
            PP[(a, d)] = [0] + [ser[-k - low] if 0 <= -k - low < len(ser) else 0 for k in range(1, 2 * d + 2)]
            der = [ser[i] * ((low + i) % P) % P for i in range(len(ser))]
            ser, low = mulser(der, Q, len(der), P), low - 2
    return PP


def local_leg(F, c, v0, n, P, ctrl=None):
    """phi(t) = Phi(c + t) to t^(n-1) with log c dropped, and g(t) = the log coefficient of Phi at c."""
    A, G = split(F)
    phi, g = [0] * n, [0] * n
    for (j, m), co in A.items():
        s = fr(co, P) * pow(v0, m, P) % P
        phi = [(x + s * y) % P for x, y in zip(phi, binser(c, j, n, P))]
    for (j, m), co in G.items():
        s = fr(co, P) * pow(v0, m, P) % P
        g = [(x + s * y) % P for x, y in zip(g, binser(c, j, n, P))]
    if ctrl != "C4":
        ci = pow(c, P - 2, P)
        lser = [0] + [(-1) ** (k + 1) * pow(ci, k, P) * pow(k, P - 2, P) % P for k in range(1, n)]
        phi = [(x + y) % P for x, y in zip(phi, mulser(g, lser, n, P))]
    return phi, g


def ancestor(F, data, v0, P, dmax, ctrl=None):
    """q[(a, d)] = Sum_j Res_{c_j} F d xif_{a,d}; also the number of non-zero L2 residues Res_c f(x) d xif_{a,d}."""
    q, l2 = {}, 0
    for c, PP in data:
        phi, g = local_leg(F, c, v0, 2 * dmax + 3, P, ctrl)
        for (a, d), pp in PP.items():
            r = -sum(k * phi[k] * pp[k] for k in range(1, len(pp))) % P
            q[(a, d)] = (q.get((a, d), 0) + r) % P
            if any(g):
                l2 += (sum(k * g[k] * pp[k] for k in range(1, len(pp))) % P) != 0
    return q, l2


def branch_data(s, P, dmax, ctrl=None):
    """Branch points c_j = s zeta^j of p^5 = v0 = s^5, with their principal parts."""
    z, v0 = zeta5(P), pow(s, 5, P)
    return v0, [(s * pow(z, j, P) % P, principal_parts(s * pow(z, j, P) % P, v0, P, dmax, ctrl)) for j in range(5)]
