# check1.py -- evidence only. Exact checks of Prop. MV1 ("(1.1)"), of the literal LLZ and OP statements, with controls.
from fractions import Fraction as Fr
from math import comb, factorial
from itertools import combinations
import sys, time
from engine import *
from schur import *


def sigma(t):
    return t * (1 + t)


def rhs(g, mu, tau):
    """sigma^(g+l-1) Pi_mu M_g(mu)/|Aut mu|"""
    l = len(mu)
    return sigma(tau) ** (g + l - 1) * PiMu(mu, tau) * Mg(g, mu, tau) / aut(mu)


def cases():
    out = []
    def plist(n, lmax):
        return [m for N in range(1, n + 1) for m in partitions(N) if len(m) <= lmax]
    out += [(0, m) for m in plist(8, 4)]
    out += [(1, m) for m in plist(8, 3)]
    out += [(2, m) for m in plist(7, 3)]
    out += [(3, m) for m in plist(7, 2)]
    out += [(3, m) for m in plist(5, 3) if len(m) == 3]
    out += [(4, m) for m in plist(5, 1)]
    return out


def main():
    t0 = time.time()
    cs = cases()
    Dmax = max(sum(m) + len(m) + 2 * g - 2 for g, m in cs)
    Pmax = max(2 * g - 2 + len(m) + sum(m) for g, m in cs)
    taus = [Fr(k + 2, 7) for k in range(Dmax + 1)]
    print(f"T2: {len(cs)} (g,mu) cases, tau sample {len(taus)} points, max degree in tau {Dmax}, series precision {Pmax}")
    conns = {}
    nok = nbad = nid = 0
    for g, m in cs:
        D = sum(m) + len(m) + 2 * g - 2
        k = 2 * g - 2 + len(m)
        good = True
        for tau in taus[:D + 1]:
            if tau not in conns:
                conns[tau] = ConnSeries(tau, Pmax)
            lhs = conns[tau].coeff(m, k)
            r = rhs(g, m, tau)
            if lhs == r:
                nok += 1
            else:
                nbad += 1
                good = False
                print("MISMATCH", g, m, tau, lhs, r)
        nid += good
    print(f"T2 (1.1): {nok} matches, {nbad} mismatches; {nid}/{len(cs)} cases hold as identities in tau (deg+1 points). t={time.time()-t0:.1f}s")
    # valuation: no power of x below 2*0-2+l, and only powers of parity l, in [p_mu]C
    nv = 0
    for m in [m for N in range(1, 8) for m in partitions(N)]:
        l = len(m)
        tau = taus[0]
        for k in range(-sum(m), 2 * 3 - 2 + l + 1):
            if k + sum(m) > Pmax:
                continue
            c = conns[tau].coeff(m, k)
            if k < l - 2 or (k - l) % 2:
                if c != 0:
                    nv += 1
                    print("valuation/parity violation", m, k, c)
    print(f"T2b valuation >= l-2 and parity x^(l mod 2): violations {nv}")
    # symmetries of the Schur side, which force the sign pattern of the Pi's and the x-parity
    nsym = 0
    tot = 0
    for tau in taus[:3]:
        C1 = ConnSeries(tau, 18)
        C2 = ConnSeries(-1 - tau, 18)
        for m in [m for N in range(1, 7) for m in partitions(N)]:
            sg = (-1) ** (sum(m) - len(m))
            for i in range(19):
                tot += 1
                if C2.C(m)[i] != sg * C1.C(m)[i]:
                    nsym += 1
    print(f"T2c [p_mu]C(x;-1-tau) = (-1)^(|mu|-l)[p_mu]C(x;tau): {tot} coefficients, violations {nsym}")

    # controls
    print("Controls (each must produce mismatches):")
    sub = [(g, m) for g, m in cs if g <= 2 and sum(m) <= 6]
    ctaus = [Fr(2, 7), Fr(3, 2), Fr(-7, 3)]
    conn = {t: ConnSeries(t, Pmax) for t in ctaus}
    conn0 = {t: ConnSeries(t, Pmax, frame_shift=Fr(0)) for t in ctaus}
    def run(name, f, want=None):
        bad = tot = 0
        for g, m in sub:
            if want and not want(g, m):
                continue
            for t in ctaus:
                tot += 1
                k = 2 * g - 2 + len(m)
                if f(g, m, t, k):
                    bad += 1
        print(f"  {name}: {bad} mismatches of {tot}")
    run("true formula (expected 0)", lambda g, m, t, k: conn[t].coeff(m, k) != rhs(g, m, t))
    run("K1 drop sigma^g", lambda g, m, t, k: conn[t].coeff(m, k) != sigma(t) ** (len(m) - 1) * PiMu(m, t) * Mg(g, m, t) / aut(m), lambda g, m: g >= 1)
    def flipped(g, m, t):
        if g == 0:
            return Mg(g, m, t)
        return leg_integral(g, m, (Fr(-1), 1 / t, -1 / (1 + t)))
    run("K2 weights (1,-tau,1+tau)", lambda g, m, t, k: conn[t].coeff(m, k) != sigma(t) ** (g + len(m) - 1) * PiMu(m, t) * flipped(g, m, t) / aut(m), lambda g, m: g >= 1)
    def flip1(g, m, t):
        if g == 0:
            return Mg(g, m, t)
        return leg_integral(g, m, (Fr(1), -1 / t, 1 / (1 + t)))
    run("K2b weight (-1,tau,-1-tau) i.e. Lambda(+1)", lambda g, m, t, k: conn[t].coeff(m, k) != sigma(t) ** (g + len(m) - 1) * PiMu(m, t) * flip1(g, m, t) / aut(m), lambda g, m: g >= 1)
    run("K3 framing e^(tau kappa x/2), no 1/2", lambda g, m, t, k: conn0[t].coeff(m, k) != rhs(g, m, t))
    run("K4 Pi_mu(-1-tau)", lambda g, m, t, k: conn[t].coeff(m, k) != sigma(t) ** (g + len(m) - 1) * PiMu(m, -1 - t) * Mg(g, m, t) / aut(m))
    def leafneg(g, m, t):
        if g == 0 and len(m) <= 2:
            return Mg(g, m, t)
        return leg_integral(g, tuple(-v for v in m), ts_of(t))
    run("K5 leaf 1/(1+mu psi)", lambda g, m, t, k: conn[t].coeff(m, k) != sigma(t) ** (g + len(m) - 1) * PiMu(m, t) * leafneg(g, m, t) / aut(m))
    run("K6 drop 1/|Aut mu|", lambda g, m, t, k: conn[t].coeff(m, k) != sigma(t) ** (g + len(m) - 1) * PiMu(m, t) * Mg(g, m, t), lambda g, m: aut(m) > 1)
    run("K7 M_0(m) = 1/m instead of 1/m^2", lambda g, m, t, k: conn[t].coeff(m, k) != sigma(t) ** (len(m) - 1) * PiMu(m, t) * Fr(1, m[0]) / aut(m), lambda g, m: g == 0 and len(m) == 1)
    run("K8 extra (-1)^g", lambda g, m, t, k: conn[t].coeff(m, k) != (-1) ** g * rhs(g, m, t), lambda g, m: g % 2 == 1)

    t3()
    t4()
    print(f"total time {time.time()-t0:.1f}s")


# ---------------- literal LLZ, in lambda, with Gaussian rationals --------------------------------------
class G:
    __slots__ = ("a", "b")
    def __init__(self, a=0, b=0):
        self.a = Fr(a); self.b = Fr(b)
    def __add__(self, o):
        return G(self.a + o.a, self.b + o.b)
    def __sub__(self, o):
        return G(self.a - o.a, self.b - o.b)
    def __mul__(self, o):
        if isinstance(o, G):
            return G(self.a * o.a - self.b * o.b, self.a * o.b + self.b * o.a)
        return G(self.a * o, self.b * o)
    def __eq__(self, o):
        return self.a == o.a and self.b == o.b
    def iszero(self):
        return self.a == 0 and self.b == 0
    def __repr__(self):
        return f"({self.a}+{self.b}i)"


def ipow(n):
    return [G(1, 0), G(0, 1), G(-1, 0), G(0, -1)][n % 4]


def gmul(a, b, P):
    c = [G() for _ in range(P + 1)]
    for i, x in enumerate(a):
        if x.iszero():
            continue
        for j, y in enumerate(b):
            if i + j > P:
                break
            if not y.iszero():
                c[i + j] = c[i + j] + x * y
    return c


def sin_series(n, P):
    """sin(n l/2)/(n l /2) truncated at l^P"""
    s = [Fr(0)] * (P + 1)
    for k in range(P // 2 + 1):
        s[2 * k] = (-1) ** k * Fr(n, 2) ** (2 * k) / factorial(2 * k + 1)
    return s


def sinv(a, P):
    b = [Fr(0)] * (P + 1)
    b[0] = 1 / a[0]
    for n in range(1, P + 1):
        b[n] = -sum(a[i] * b[n - i] for i in range(1, n + 1)) / a[0]
    return b


def V_literal(nu, P):
    """LLZ (3) literally: lambda^{|nu|} V_nu(lambda), truncated at P; product of sines, no hooks"""
    l = len(nu)
    s = [Fr(0)] * (P + 1)
    s[0] = Fr(1)
    const = Fr(1)
    for a in range(1, l + 1):
        for b in range(a + 1, l + 1):
            n1 = nu[a - 1] - nu[b - 1] + b - a
            n2 = b - a
            s = smul(s, sin_series(n1, P), P)
            s = smul(s, sinv(sin_series(n2, P), P), P)
            const *= Fr(n1, n2)
    for i in range(1, l + 1):
        for v in range(1, nu[i - 1] + 1):
            n = v - i + l
            s = smul(s, sinv(sin_series(n, P), P), P)
            const *= Fr(1, n)          # 1/(2 sin(n l/2)) = 1/(n l) * [..]^{-1}
    return [c * const for c in s]


def bracket_V(nu, P):
    """the report's line-177 V: prod_box 1/[h], [m] = q^{m/2}-q^{-m/2} = 2 i sin(m lambda/2); returned as lambda^{|nu|} V"""
    k = len(hooks(nu))
    # 1/(2 i sin) = (-i)/(2 sin), so prod = (-i)^{|nu|} V_literal
    fac = ipow(-k)
    return [fac * c for c in V_literal(nu, P)]


def R_conn_gauss(tau, P, mus, Vfun):
    """[p_mu] log of sum_nu e^{i (tau+1/2) kappa l/2} V_nu s_nu(p), as lambda^{-|mu|} series, Gaussian"""
    Zc = {}
    Cc = {}
    def Zmu(mu):
        if mu in Zc:
            return Zc[mu]
        n = sum(mu)
        tot = [G() for _ in range(P + 1)]
        zz = zee(mu)
        for nu in partitions(n):
            ch = chi(nu, mu)
            if ch == 0:
                continue
            e = (tau + Fr(1, 2)) * kappa(nu) / 2
            ex = []
            cur = G(1, 0)
            for k in range(P + 1):
                ex.append(cur)
                cur = cur * G(0, e) * Fr(1, k + 1)
            V = Vfun(nu, P)
            V = [v if isinstance(v, G) else G(v) for v in V]
            w = gmul(V, ex, P)
            for i in range(P + 1):
                tot[i] = tot[i] + w[i] * Fr(ch, zz)
        Zc[mu] = tot
        return tot
    def Cmu(mu):
        if mu in Cc:
            return Cc[mu]
        n = sum(mu)
        tot = list(Zmu(mu))
        for rho in sub_multisets(mu):
            if not rho or rho == mu:
                continue
            t = msub(mu, rho)
            f = Fr(sum(rho), n)
            pr = gmul(Cmu(rho), Zmu(t), P)
            for i in range(P + 1):
                tot[i] = tot[i] - pr[i] * f
        Cc[mu] = tot
        return tot
    return Cmu


def I_literal(g, mu, tau):
    """LLZ's integral: Int Lambda^vee(1) Lambda^vee(-tau-1) Lambda^vee(tau) / prod (1 - mu_i psi_i),
       Lambda^vee(u) = u^g - lambda_1 u^(g-1) + ... ; expanded in lambda_i lambda_j lambda_k literally."""
    l = len(mu)
    if g == 0:
        return Fr(sum(mu)) ** (l - 3)
    dimtot = 3 * g - 3 + l
    tot = Fr(0)
    us = (Fr(1), -tau - 1, tau)
    for i in range(g + 1):
        for j in range(g + 1):
            for k in range(g + 1):
                if i + j + k > dimtot:
                    continue
                c = (-1) ** (i + j + k) * us[0] ** (g - i) * us[1] ** (g - j) * us[2] ** (g - k)
                mi = omega_monos(g, i, (Fr(1),))
                mj = omega_monos(g, j, (Fr(1),))
                mk = omega_monos(g, k, (Fr(1),))
                for a in compositions(dimtot - i - j - k, l):
                    w = Fr(1)
                    for m, e in zip(mu, a):
                        w *= Fr(m) ** e
                    for c1, x1 in mi:
                        for c2, x2 in mj:
                            for c3, x3 in mk:
                                tot += c * w * x1 * x2 * x3 * hodge(g, a, c1 + c2 + c3)
    return tot


def t3():
    print("T3: LLZ eqs (1)-(3) literally, in lambda over Q(i), against -(i)^(|mu|+l)/|Aut| sigma^(l-1) Pi_mu I_g:")
    cs = [(g, m) for g in (0, 1, 2) for N in range(1, 6) for m in partitions(N) if len(m) <= 3]
    P = 2 * 2 - 2 + 3 + 5
    for name, Vfun, expo in (("printed LLZ: V with 2 sin, prefactor i^(|mu|+l)", V_literal, lambda m: sum(m) + len(m)),
                             ("A_mv line 177 variant: V = prod 1/[h], prefactor i^l", bracket_V, lambda m: len(m)),
                             ("CONTROL: V with 2 sin, prefactor i^l (mixed)", V_literal, lambda m: len(m))):
        ok = bad = 0
        for tau in (Fr(2, 7), Fr(-5, 3)):
            C = R_conn_gauss(tau, P, None, Vfun)
            for g, m in cs:
                l = len(m)
                k = 2 * g - 2 + l + sum(m)
                lhs = C(m)[k]
                r = ipow(expo(m)) * (Fr(-1) * sigma(tau) ** (l - 1) * PiMu(m, tau) * I_literal(g, m, tau) / aut(m))
                if lhs == r:
                    ok += 1
                else:
                    bad += 1
        print(f"  {name}: {ok} matches, {bad} mismatches")


def t4():
    print("T4: OP Theorem 1 (0.8) literally: disconnected series, parameters (-1,-1/a,1/(a+1)), sqrt(a(a+1))u:")
    from itertools import product
    def setparts(items):
        if not items:
            yield []
            return
        first, rest = items[0], items[1:]
        for sp in setparts(rest):
            yield [[first]] + sp
            for i in range(len(sp)):
                yield sp[:i] + [[first] + sp[i]] + sp[i + 1:]
    def Hconn(g, z, a):
        """H^o_g(z; -1,-1/a,1/(a+1)) = prod z_i Int ... with OP's unstable conventions"""
        pr = Fr(1)
        for v in z:
            pr *= v
        if g == 0 and len(z) == 1:
            return Fr(1, z[0])
        if g == 0 and len(z) == 2:
            return pr / (z[0] + z[1])
        return pr * leg_integral(g, tuple(z), (Fr(-1), -1 / a, 1 / (a + 1)))
    ok = bad = 0
    ok28 = bad28 = autfac = 0
    ctrl_ok = ctrl_bad = 0
    K = 2
    for a in (Fr(2, 7), Fr(3, 2)):
        for N in range(1, 7):
            for mu in partitions(N):
                l = len(mu)
                if l > 3:
                    continue
                pref = Fr(1)
                for m in mu:
                    pref *= Fr(comb(0, 0))
                    # binomial((a+1)m, m) = prod_{b=1}^m ((a+1) m - m + b)/m!
                    b_ = Fr(1)
                    for b in range(1, m + 1):
                        b_ *= ((a + 1) * m - m + b)
                    pref *= b_ / factorial(m)
                # left side, coefficient of u^k, k <= K
                L = {}
                for sp in setparts(list(range(l))):
                    blocks = [tuple(mu[i] for i in B) for B in sp]
                    def rec(i, pw, val):
                        if i == len(blocks):
                            L[pw] = L.get(pw, Fr(0)) + val
                            return
                        z = blocks[i]
                        remaining = len(blocks) - i - 1
                        g = 0
                        while True:
                            e = 2 * g - 2
                            if pw + e - 2 * remaining > K:
                                break
                            if not (g == 0 and len(z) >= 3 and False):
                                v = Hconn(g, z, a) * (a * (a + 1)) ** (g - 1)
                                rec(i + 1, pw + e, val * v)
                            g += 1
                    rec(0, 0, pref)
                # right side: (a u)^(-l)/|Aut| sum_lambda hook(u) e^{(a+1/2) u f2} chi; hook(u) = u^{-|mu|} series
                P = K + l + N
                S = [Fr(0)] * (P + 1)
                for nu in partitions(N):
                    ch = chi(nu, mu)
                    if ch == 0:
                        continue
                    w = Wseries(nu, a, P)          # uses f2 = kappa/2
                    for i in range(P + 1):
                        S[i] += ch * w[i]
                c0 = a ** (-l) / aut(mu)
                for k in range(-l - N, K + 1):
                    r = c0 * S[k + l + N]
                    lv = L.get(k, Fr(0))
                    if r == lv:
                        ok += 1
                    else:
                        bad += 1
                        ratio = lv / r
                        if ratio != aut(mu):
                            print("OP mismatch not explained by |Aut mu|", a, mu, k, r, lv)
                        else:
                            autfac += 1
                    # the operator form OP (1.28) has no 1/|Aut mu|
                    if a ** (-l) * S[k + l + N] == lv:
                        ok28 += 1
                    else:
                        bad28 += 1
                    # control: without the sqrt(a(a+1)) rescaling, compare only where it matters
                # control run: drop the rescaling factor (a(a+1))^(g-1) in one test coefficient
                Lc = {}
                for sp in setparts(list(range(l))):
                    blocks = [tuple(mu[i] for i in B) for B in sp]
                    def rec2(i, pw, val):
                        if i == len(blocks):
                            Lc[pw] = Lc.get(pw, Fr(0)) + val
                            return
                        z = blocks[i]
                        remaining = len(blocks) - i - 1
                        g = 0
                        while True:
                            e = 2 * g - 2
                            if pw + e - 2 * remaining > K:
                                break
                            rec2(i + 1, pw + e, val * Hconn(g, z, a))
                            g += 1
                    rec2(0, 0, pref)
                for k in range(-l - N, K + 1):
                    r = c0 * S[k + l + N]
                    if r == Lc.get(k, Fr(0)):
                        ctrl_ok += 1
                    else:
                        ctrl_bad += 1
    print(f"  OP (0.8) as printed, labelled disconnected series: {ok} coefficients match, {bad} mismatches, of which {autfac} are exactly the factor |Aut mu| (repeated parts only)")
    print(f"  OP (0.8) with no 1/|Aut mu|, as in OP (1.24), (1.28): {ok28} match, {bad28} mismatches")
    print(f"  CONTROL, no rescaling sqrt(a(a+1)): {ctrl_bad} mismatches of {ctrl_ok+ctrl_bad} coefficients")


if __name__ == "__main__":
    main()
