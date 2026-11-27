# -*- coding: utf-8 -*-
"""rf_g2flat.py -- the genus-2 graph sum of P.Omega^mod_{v(u)} computed in the FLAT frame (the form in which Theorem R of
REDUCTION_DRAFT_2.md s.3.2 is stated), graph by graph, exactly over Q, with the canonical Theorem-U gauge P of the formal quintic.

Nobody computed the flat-frame sum before: the project and the univ worker compute the HATTED outer sum R'.(M.omega^).  Agreement
of the seven genus-2 graph values with the hatted values on record (UNIFORM_REDUCTION.md s.6.2: vertex u^-2(1/1200 - 7u/6000 + ..),
edge graphs O(u^2) with [u^2] = 479/756000, 1/2250, -127/2016000, 0, 0, 0) tests the identification of s.1.5 of REDUCTION_DRAFT.md
(flat = hatted, graph by graph) and Theorem R's conclusion at genus 2 in the flat frame.
Ingredients (own code): P (rf_flat_core: ur_gauge read-only + canonical normalisation), vector form, flat edge E_P, leaf point
q_T = Pv^-1(-z phi_0) = -(z/d) phi_0 + T, vertices d^(2h-2+n) Sum_j <x.. (dT)^j>_h / j! with model correlators C v(u)^e, C from
MVertex6 fed with the referee's own LG M_0 (rf_mla_core.Corr) at a 128-bit prime, rationally reconstructed and re-checked at 2^89-1.
Controls: (X1) translation T dropped (leaves removed), (X2) vertex evaluated at v(u)(1 + u) instead of v(u), (X3), (X4) a left
dressing P -> (1 + u^m z E_21) P (lowered entry (2,1), m = 0, 1): flat-regular, graded, symplectic, but (P Y^C)_+ = -u^m z e_2 + ..: (A) fails.
E3-type test: P -> exp(u^m z N) P, N = E_04 + E_43 (Pi-type, infinitesimally symplectic), m = 0, 1: flat-regular, graded, symplectic and
annihilating but not a gauge; Theorem R predicts regular edge graphs and vertex (1/1200) T_P^-2 + regular, T_P = -[z^-1](P Y^C)_0.
Run: py -I -B rf_g2flat.py [a b]   Log: rf_g2flat.log (appended).
"""
import os, sys, time, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rf_flat_core as C
from rf_flat_core import yc_series, sinv, spow, smul, sadd, sscale, sc, s0, zinv, eta_conj, mT, mmul, mscale, madd, mzero, tau_v, load_gauge, canonical
from rf_mla_core import Corr, w
from rf_mla import is_prime
UMIN, UMAX = -6, 5          # Laurent window in u


def lser(s, shift=0):       # u-series list -> Laurent dict
    return {i + shift: x for i, x in enumerate(s) if x and UMIN <= i + shift <= UMAX}
def lmul(a, b):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            if i + j <= UMAX: out[i + j] = out.get(i + j, 0) + x * y
    return {k: x for k, x in out.items() if x}
def ladd(a, b, c=1):
    out = dict(a)
    for k, x in b.items(): out[k] = out.get(k, 0) + c * x
    return {k: x for k, x in out.items() if x}
def lsc(a, c): return {k: c * x for k, x in a.items() if c * x}


def ratrec(a, p):
    a %= p; r0, r1, s0_, s1 = p, a, 0, 1; B = int((p // 2) ** 0.5)
    while r1 > B:
        q = r0 // r1; r0, r1 = r1, r0 - q * r1; s0_, s1 = s1, s0_ - q * s1
    assert s1 != 0 and abs(s1) <= B, "ratrec failed"
    return Fr(r1, s1)


class Setup(object):
    def __init__(self, a, b, log, vshift=False, noT=False, dress=None, dressN=None):
        C.NU[0] = 10; K = 6
        tau, v, itau = tau_v()
        Ps, tc, var0, res = load_gauge(K, a, b)
        P, _ = canonical(Ps, K)
        if dress is not None:            # P -> (1 + u^m z E_21) P
            m = dress; P2 = {}
            for k in range(K + 1):
                M = [[list(x) for x in r] for r in P[k]]
                if k >= 1:
                    for j in range(5):
                        M[2][j] = sadd(M[2][j], [Fr(0)] * m + P[k - 1][1][j][:C.NU[0] - m])
                P2[k] = M
            P = P2
            from rf_flat import sympl_orders
            log("  dressed gauge: symplectic defects %s" % sympl_orders(P, K)[:3])
        if dressN is not None:           # P -> exp(u^m z N) P, N = E_04 + E_43 (lowered; N^2 = E_03, N^3 = 0): annihilating, NOT a gauge
            m = dressN; P2 = {}
            for k in range(K + 1):
                M = [[list(x) for x in r] for r in P[k]]
                for j in range(5):
                    if k >= 1:
                        M[0][j] = sadd(M[0][j], [Fr(0)] * m + P[k - 1][4][j][:C.NU[0] - m])
                        M[4][j] = sadd(M[4][j], [Fr(0)] * m + P[k - 1][3][j][:C.NU[0] - m])
                    if k >= 2:
                        M[0][j] = sadd(M[0][j], [Fr(0)] * (2 * m) + [x / 2 for x in P[k - 2][3][j][:C.NU[0] - 2 * m]])
                P2[k] = M
            P = P2
            from rf_flat import sympl_orders
            log("  dressed gauge exp(u^%d z N): symplectic defects %s" % (m, sympl_orders(P, K)[:3]))
        Y = yc_series(v, 12); te = s0()
        for k in range(K + 1):
            for j in range(5): te = sadd(te, smul(P[k][0][j], Y[k + 1][j]))
        self.teff = sscale(te, Fr(-1))
        Pi = zinv(P, K); Pvi = {k: eta_conj(Pi[k]) for k in Pi}
        self.d = P[0][3][3]
        self.dT = {}
        for k in range(1, K + 1):                       # q_T at z^(k+1) = -(Pv^-1)_k e_0;  T = q_T + (z/d) e_0 at z^(>=2)
            for c in range(5):
                s = sscale(Pvi[k][c][0], Fr(-1))
                if any(s): self.dT[(k + 1, c)] = lser(smul(self.d, s))
        EI = mzero()
        for i in range(5): EI[i][(3 - i) % 5] = sc(Fr(1, 5))
        Num = {}
        for m in range(K + 1):
            for n in range(K + 1 - m): Num[(m, n)] = mscale(mmul(mmul(Pvi[m], EI), mT(Pvi[n])), Fr(-1))
        Num[(0, 0)] = madd(Num[(0, 0)], EI)
        self.E = {}
        for m in range(4):
            for n in range(4):
                if m + n > K - 1: continue
                acc = mzero()
                for t in range(n + 1): acc = madd(acc, Num[(m + 1 + t, n - t)], (-1) ** t)
                for x in range(5):
                    for y in range(5):
                        if any(acc[x][y]): self.E[((m, x), (n, y))] = lser(acc[x][y])
        vv = v if not vshift else smul(v, [Fr(1), Fr(1)] + [Fr(0)] * (C.NU[0] - 2))
        self.w = vv[1:] + [Fr(0)]                       # v = u w
        self.dl = {e: lser(spow(self.d, Fr(e))) for e in range(-3, 7)}
        if noT: self.dT = {}
        self.vp = {}
        p = (1 << 127) + 1
        while not (p % 10 == 1 and is_prime(p)): p += 1
        self.Cp = Corr(p); self.Cq = Corr((1 << 89) - 1); self.cache = {}
        self.tc = tc
        log("  setup: d = %s ..., dT components %d, edge entries %d, prime %d" % (self.d[:3], len(self.dT), len(self.E), p))

    def vpow(self, e):
        if e not in self.vp: self.vp[e] = lser(spow(self.w, Fr(e)), e)
        return self.vp[e]

    def corr(self, h, ins):
        key = (h, tuple(sorted(ins)))
        if key not in self.cache:
            n = len(ins)
            if 2 * h - 2 + n <= 0 or 3 * h - 3 + n - sum(k for k, a in ins) < 0: val = Fr(0)
            else:
                x = self.Cp.flat(h, list(ins)); val = ratrec(x, self.Cp.p) if x else Fr(0)
                y = self.Cq.flat(h, list(ins)); assert (self.Cq.R.of(val) - y) % self.Cq.p == 0, "ratrec cross-check failed"
            self.cache[key] = val
        return self.cache[key]

    def vertex(self, h, labels):
        """d^(2h-2+n) Sum_j (1/j!) Sum_{ordered j-tuples of dT} <labels, tuple>_h (v)"""
        n = len(labels); out = {}
        Tk = list(self.dT.items()); jmax = 3 * h - 3 + n - sum(k for k, a in labels)
        for j in range(0, max(jmax, 0) + 1):
            fj = Fr(1, 1)
            for i in range(2, j + 1): fj /= i
            for tup in itertools.product(Tk, repeat=j):
                ins = list(labels) + [kc for kc, _ in tup]
                if 3 * h - 3 + len(ins) - sum(k for k, a in ins) < 0: continue
                cval = self.corr(h, ins)
                if not cval: continue
                e = 2 - 2 * h + sum(w(k, a) for k, a in ins)
                assert e.denominator == 1
                term = lsc(self.vpow(int(e)), cval * fj)
                for _, ser in tup: term = lmul(term, ser)
                out = ladd(out, term)
        return lmul(out, self.dl[2 * h - 2 + n])


def graphs(S, log):
    E = S.E; lab = sorted(set(x for (x, y) in E))
    def ev(X, Y): return E.get((X, Y), {})
    out = {}
    out['0 vertex g=2'] = S.vertex(2, [])
    acc = {}
    for (X, Y), e in E.items(): acc = ladd(acc, lmul(e, S.vertex(1, [X, Y])), Fr(1, 2))
    out['1 loop on g=1'] = acc
    V1 = {X: S.vertex(1, [X]) for X in lab}; acc = {}
    for (X, Y), e in E.items(): acc = ladd(acc, lmul(lmul(e, V1[X]), V1[Y]), Fr(1, 2))
    out['2 bridge g1-g1'] = acc
    acc = {}
    for (X, Y), e in E.items():
        for (Z, W), f in E.items(): acc = ladd(acc, lmul(lmul(e, f), S.vertex(0, [X, Y, Z, W])), Fr(1, 8))
    out['3 two loops on g=0'] = acc
    acc = {}
    for (X, Y), e in E.items():
        for (Z, W), f in E.items():
            if not V1[W]: continue
            acc = ladd(acc, lmul(lmul(lmul(e, f), S.vertex(0, [X, Y, Z])), V1[W]), Fr(1, 2))
    out['4 g0(loop)-g1'] = acc
    V0 = {}
    def v0(t):
        k = tuple(sorted(t))
        if k not in V0: V0[k] = S.vertex(0, list(k))
        return V0[k]
    acc = {}
    items = list(E.items())
    for (X1, Y1), e1 in items:
        for (X2, Y2), e2 in items:
            e12 = lmul(e1, e2)
            if not e12: continue
            for (X3, Y3), e3 in items:
                a = v0((X1, X2, X3))
                if not a: continue
                b = v0((Y1, Y2, Y3))
                if not b: continue
                acc = ladd(acc, lmul(lmul(e12, e3), lmul(a, b)), Fr(1, 12))
    out['5 g0=g0 triple edge'] = acc
    acc = {}
    for (X1, Y1), e1 in items:
        for (Z, W), f in items:
            a = v0((X1, Y1, Z))
            if not a: continue
            for (X2, Y2), e2 in items:
                b = v0((W, X2, Y2))
                if b: acc = ladd(acc, lmul(lmul(lmul(e1, f), e2), lmul(a, b)), Fr(1, 8))
    out['6 g0(loop)-g0(loop)'] = acc
    return out


def main():
    a = Fr(sys.argv[1]) if len(sys.argv) > 1 else Fr(0); b = Fr(sys.argv[2]) if len(sys.argv) > 2 else Fr(0)
    LOG = open(os.path.join(HERE, "rf_g2flat.log"), "a", encoding="ascii")
    def log(s): print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    log("=" * 110); log("rf_g2flat: frame (%s, %s) [%s]" % (a, b, time.strftime("%Y-%m-%dT%H:%M:%S")))
    for tag, kw in (("canonical gauge", {}), ("control X1: no translation", dict(noT=True)), ("control X2: vertex at v(1+u)", dict(vshift=True)),
                    ("control X3: (1 + z E_21) P", dict(dress=0)), ("control X4: (1 + u z E_21) P", dict(dress=1)),
                    ("E3-type: exp(z N) P (annihilating non-gauge)", dict(dressN=0)), ("E3-type: exp(u z N) P", dict(dressN=1))):
        t0 = time.time(); S = Setup(a, b, log, **kw); G = graphs(S, log); tot = {}
        log("  --- %s: graph values u^-2 (G_0 + G_1 u + G_2 u^2 + G_3 u^3)" % tag)
        for name, ser in G.items():
            tot = ladd(tot, ser)
            log("    %-22s G = %s" % (name, [str(ser.get(k, 0)) for k in range(-2, 2)]))
        log("    %-22s G = %s" % ("TOTAL", [str(tot.get(k, 0)) for k in range(-2, 2)]))
        tcm2 = lser(spow(S.tc[1:] + [Fr(0)], Fr(-2)), -2); tem2 = lser(spow(S.teff[1:] + [Fr(0)], Fr(-2)), -2)
        log("    (1/1200) t_c^-2       G = %s" % [str(Fr(1, 1200) * tcm2.get(k, 0)) for k in range(-2, 2)])
        log("    (1/1200) T_P^-2       G = %s   T_P = -[z^-1](P Y^C)_0 = %s ...  (%.0f s)"
            % ([str(Fr(1, 1200) * tem2.get(k, 0)) for k in range(-2, 2)], [str(x) for x in S.teff[:4]], time.time() - t0))
    LOG.close()


if __name__ == "__main__":
    main()
