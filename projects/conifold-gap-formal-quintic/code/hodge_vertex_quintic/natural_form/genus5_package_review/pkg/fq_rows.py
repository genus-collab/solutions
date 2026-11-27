# -*- coding: utf-8 -*-
"""fq_rows.py -- WKB rows of the formal quintic at the zeta point (PLAN step 2), in the differential ring N.

Part 1 (this file head): operator algebra over F_p[L^+-1] and the WKB recursion with residue checks.
  Laurent polynomials in L: dict {Lpow: coeff};  D L^k = (k/5)(L^(k+5) - L^k)  (DL = (L^6-L)/5).
  Operators: dict {(zpow, Dord): LP}, composition by Leibniz (D o f = f D + Df).
  PF at colour i, conjugated by e^{u_i/z} (D u_i = l := zeta^i L):
     Pcal = (M+l)^5 - 1 - q [3125 (M+l)^5 + 9375 z (M+l)^4 + 10625 z^2 (M+l)^3 + 5625 z^3 (M+l)^2 + 1370 z^4 (M+l) + 120 z^5],
     M = zD, q = (1 - L^-5)/3125;  Pcal = Sum_{m>=1} z^m Pcal_m  (Pcal_0 = 0 is the symbol equation l^5 (1-3125q) = 1).
  WKB:  I_i = e^{u_i/z} z psi_i,  psi_i = Sum_n z^n psi_{i,n},  psi_{i,0} = L (from Pcal_1 psi_0 = 0, a_1 = 5 zeta^-i / L),
     and for n >= 1:  D(psi_n/psi_0) = r_n := -(zeta^i/5) Sum_{m=2}^{n+1} Pcal_m psi_{n+1-m}   (Laurent polynomial in L).
     Integration: d(phi_n)/dL = 5 r_n/(L (L^5-1)); r_n = (L^5-1) s_n + rem_n: residues at the fifth roots of unity are the
     values rem_n(zeta^j) (all zero iff rem_n = 0), residue at L = 0 is 5 [L^0] s_n.  Both must vanish (else a log).
     phi_n = 5 Sum_{k != 0} s_{n,k} L^k / k + const, const fixed by phi_n(L=1) = [z^n] h_{wt_i}(z)  (q = 0 <-> L = 1,
     R|_{q=0} = diag(h_{wt_i}): the WKB solution from the I-function has psi(q=0) = 1, and R = R~ h_{wt_i}).
Part 2 (appended): rows P_a^{(i)} in N, Tier-1 q-series rows, pairing/edge helpers, Dil.
"""
import os, sys, math
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fq_series import Ser
from fq_ring import Ring, NEl, IDX, NAMES
from fq_hodge import hv, wt, W0


# ------------------------------------------------------------------ Laurent polynomials in L over F_p
def lp_add(a, b, p):
    out = dict(a)
    for k, v in b.items():
        x = (out.get(k, 0) + v) % p
        if x: out[k] = x
        elif k in out: del out[k]
    return out
def lp_scale(a, c, p):
    c %= p
    return {k: v * c % p for k, v in a.items() if v * c % p}
def lp_mul(a, b, p):
    out = {}
    for k1, v1 in a.items():
        for k2, v2 in b.items():
            k = k1 + k2
            x = (out.get(k, 0) + v1 * v2) % p
            if x: out[k] = x
            elif k in out: del out[k]
    return out
def lp_D(a, p, inv5):
    out = {}
    for k, v in a.items():
        if k == 0: continue
        c = v * k % p * inv5 % p
        for kk, cc in ((k + 5, c), (k, (-c) % p)):
            x = (out.get(kk, 0) + cc) % p
            if x: out[kk] = x
            elif kk in out: del out[kk]
    return out
def lp_eval(a, x, R):
    return sum((v * R.pow(x, k) for k, v in a.items()), 0) % R.p
def lp_fmt(a, R):
    return " + ".join("(%s) L^%d" % (R.fmt(v), k) for k, v in sorted(a.items())) if a else "0"


# ------------------------------------------------------------------ operators sum_{m,r} z^m a_{m,r}(L) D^r
def op_add(A, B, p):
    out = dict(A)
    for k, v in B.items():
        s = lp_add(out.get(k, {}), v, p)
        if s: out[k] = s
        elif k in out: del out[k]
    return out
def op_scale(A, c, p): return {k: lp_scale(v, c, p) for k, v in A.items() if lp_scale(v, c, p)}
def op_mul(A, B, p, inv5):
    """A o B (operator composition)"""
    out = {}
    for (m, r), a in A.items():
        for (n, s), b in B.items():
            Dtb = b
            for t in range(r + 1):
                coef = lp_scale(lp_mul(a, Dtb, p), math.comb(r, t), p)
                key = (m + n, r - t + s)
                if coef:
                    x = lp_add(out.get(key, {}), coef, p)
                    if x: out[key] = x
                    elif key in out: del out[key]
                Dtb = lp_D(Dtb, p, inv5)
    return out
def op_pow(A, n, p, inv5):
    out = {(0, 0): {0: 1}}
    for _ in range(n): out = op_mul(out, A, p, inv5)
    return out
def op_apply(A, f, p, inv5):
    """apply to a Laurent polynomial f(L) (z-free); returns {zpow: LP}"""
    out = {}
    Dk = {}
    maxr = max((r for (m, r) in A), default=0)
    g = f
    for r in range(maxr + 1):
        Dk[r] = g; g = lp_D(g, p, inv5)
    for (m, r), a in A.items():
        t = lp_mul(a, Dk[r], p)
        if t:
            x = lp_add(out.get(m, {}), t, p)
            if x: out[m] = x
            elif m in out: del out[m]
    return out


class WKB(object):
    """WKB solution psi_i = Sum_n z^n psi_{i,n}(L) at colour i (zeta point), with residue checks.
    .psi[n]: Laurent polynomials in L (dict), n = 0..Nz;  .phi[n] = psi_n/L;  .residues[n] = (rem_n dict, res0);
    .h: the QRR constants h_{wt_i}(z) coefficients used as the q = 0 values (phi_n(1) = h_n)."""
    def __init__(self, R, i, Nz, lams=None):
        self.R = R; self.i = i; self.Nz = Nz; p = R.p; inv5 = R.inv(5)
        if lams is None: lams = [R.z(k) for k in range(5)]
        self.lams = lams
        l = lams[i]                                 # l = zeta^i (times L)
        # operator Ml = zD + l L
        Ml = {(1, 1): {0: 1}, (0, 0): {1: l}}
        q = {0: R.inv(3125), -5: (-R.inv(3125)) % p}
        Pw = [op_pow(Ml, k, p, inv5) for k in range(6)]
        P = op_add(Pw[5], {(0, 0): {0: p - 1}}, p)
        for k, c, zp in ((5, 3125, 0), (4, 9375, 1), (3, 10625, 2), (2, 5625, 3), (1, 1370, 4), (0, 120, 5)):
            t = {(m + zp, r): lp_mul(lp_scale(a, c, p), q, p) for (m, r), a in Pw[k].items()}
            P = op_add(P, op_scale(t, p - 1, p), p)
        self.P = P
        # split by z-power
        self.Pm = {}
        for (m, r), a in P.items():
            self.Pm.setdefault(m, {})[(0, r)] = a
        assert 0 not in self.Pm, "symbol equation fails: Pcal_0 != 0: %s" % self.Pm.get(0)
        P1 = self.Pm[1]
        assert set(P1.keys()) <= {(0, 0), (0, 1)}
        self.a1 = P1[(0, 1)]
        assert self.a1 == {-1: 5 * R.inv(l) % p}, "a_1 != 5 zeta^-i / L: %s" % self.a1
        # psi_0 = L: check Pcal_1 L = 0
        self.psi = [{1: 1}]
        assert not op_apply(P1, {1: 1}, p, inv5), "Pcal_1 psi_0 != 0"
        self.h = hv(R, wt(R, lams, i), Nz)
        self.phi = [{0: 1}]
        self.residues = {}
        self.r = {}
        self.sfac = R.mul(l, inv5)                  # zeta^i / 5
        for n in range(1, Nz + 1):
            rhs = {}
            for m in range(2, n + 2):
                Pm = self.Pm.get(m)
                if not Pm: continue
                t = op_apply(Pm, self.psi[n + 1 - m], p, inv5).get(0, {})
                rhs = lp_add(rhs, t, p)
            r = lp_scale(rhs, (-self.sfac) % p, p)      # D phi_n = r
            self.r[n] = r
            # r = (L^5 - 1) s + rem, rem supported on 5 consecutive exponents [kmin, kmin+4]
            s, rem = self._divide(r, p)
            res0 = 5 * s.get(0, 0) % p
            self.residues[n] = (rem, res0)
            if rem or res0:
                self.failed = n
                self.phi.append(None); self.psi.append(None)
                continue
            phi = {}
            for k, v in s.items():
                if k == 0: continue
                phi[k] = 5 * v % p * R.inv(k % p) % p
            val1 = sum(phi.values()) % p
            phi[0] = (self.h[n] - val1) % p
            phi = {k: v for k, v in phi.items() if v}
            self.phi.append(phi)
            self.psi.append(lp_mul(phi, {1: 1}, p))
        self.failed = getattr(self, 'failed', None)

    @staticmethod
    def _divide(r, p):
        """r = (L^5 - 1) s + rem with rem supported in [kmin, kmin+4]"""
        if not r: return {}, {}
        kmin = min(r); kmax = max(r)
        r = dict(r); s = {}
        for k in range(kmax, kmin + 4, -1):
            v = r.get(k, 0)
            if not v: continue
            s[k - 5] = v
            r[k] = 0
            r[k - 5] = (r.get(k - 5, 0) + v) % p
        rem = {k: v for k, v in r.items() if v}
        return s, rem


# ------------------------------------------------------------------ Part 2: rows in N
LKEY = lambda k: (k, 0, 0, 0, 0, 0, 0, 0, 0, 0)          # key of L^k
def kmul(num, key, coeff, p):
    """numerator dict times the monomial coeff * key"""
    return {tuple(x + y for x, y in zip(k, key)): v * coeff % p for k, v in num.items()}


class Rows(object):
    """Rows of the formal quintic at the zeta point in the ring N (PLAN step 2), flat basis H^a, a = 0..4.
    Conventions (FRAMEB s.0.1-0.2, A2-A3; loc_check.py signs):
      R_{ai}(z) = e^{-u_i/z} [S(H^a)]_i  normalised with R|_{q=0} = diag(h_{wt_i}) in the f-basis
      (f_i = Delta_i^{1/2} e_i, eta(f_i,f_j) = delta_ij), i.e. R_{ai}(q=0) = h_{wt_i}(z) zeta^{ia} (W_i^0)^{-1/2}.
      P_a^{(i)}(z) := Delta_i^{1/2} R_{ai}(z),  Delta_i^{1/2} = C0 zeta^{3i/2}/L with zeta^{1/2} := zeta^3, so that
          P_0^{(i)} = psi_i/psi_{i,0} = psi_i/L,     P_{a+1}^{(i)} = (1/(C1 kappa_a)) (zD + zeta^i L) P_a^{(i)} * (L/C0) ...
      precisely  P_a^{(i)} = (C0/L) X_a,  X_0 = psi_i/C0,  X_{a+1} = (1/(C1 kappa_a)) (zD + zeta^i L) X_a,
      kappa = (1, L^5/(C0^2 C1^3), 1, C0/C1, C0/C1)  (H*H^a = kappa_a H^{a+1}; genus0_check.py).
      P_a^{(i)}(q=0) = zeta^{ia} h_{wt_i}(z).   Classical pairing eta_{ab} = 5 delta_{a+b = 3 mod 5},
      eta^{ab} = (1/5) delta_{a+b = 3 mod 5} (pairs (0,3),(1,2),(2,1),(3,0),(4,4)).
      Symplecticity:  Sum_{ab} P_a^{(i)}(z) eta^{ab} P_b^{(j)}(-z) = delta_ij Delta_i,  Delta_i = zeta^{3i} C0^2/L^2.
      Dil_i(d) = (-1)^d [z^{d-1}] P_0^{(i)}  (d >= 2), Laurent polynomial in L.
      Frame change: rows_times(hser) multiplies every P_a^{(i)} by a z-series (colour-independent list of coefficients,
      or a dict i -> list), e.g. hv(R, [lam], Nz, -1) for R' = R h_{(lam)}^{-1}.
    API:  Rows(R, Nz=6, reduce=True).P[a][i] = list over z-order n of NEl (reduced modulo the relation);
          .Delta(i) -> NEl;  .eta_inv -> dict (a,b) -> field element;  .Dil(i, d) -> dict Lpow -> coeff;
          .pair(i, j) -> list over z-order of NEl (the symplecticity form);  .wkb[i] -> WKB objects;
          .q0_value(el) -> field element at q = 0;  .rows_times(h) -> new Rows-like dict P."""
    def __init__(self, R, Nz=6, reduce=True, ring=None, colours=range(5)):
        self.R = R; self.Nz = Nz; p = R.p
        self.ring = ring if ring is not None else Ring(R)
        if self.ring.i3num is None: self.ring.derive_relation()
        ring = self.ring
        self.lams = [R.z(k) for k in range(5)]
        self.wkb = {i: WKB(R, i, Nz) for i in colours}
        for i in colours: assert self.wkb[i].failed is None, "WKB residue failure at colour %d order %d" % (i, self.wkb[i].failed)
        # 1/(C1 kappa_a) as NEl
        c0 = ring.gen('c0'); W = ring.W()
        self.kinv = [ring.C1inv(),                                             # a=0: 1/C1 = c0^2/W
                     NEl(R, kmul(W.mul(W).num, LKEY(-5), 1, p), 2, 0),          # a=1: W^2 L^-5 / c0^2
                     ring.C1inv(),                                             # a=2: 1/C1
                     ring.C0inv()]                                             # a=3: 1/c0
        self.P = [[None] * 5 for _ in range(5)]
        for i in colours:
            X = [NEl(R, {LKEY(k): v for k, v in self.wkb[i].psi[n].items()}, 1, 0) for n in range(Nz + 1)]   # psi_{i,n}/c0
            l = self.lams[i]
            for a in range(5):
                # P_a = (c0/L) X_a
                row = [NEl(R, kmul(x.num, (-1, 1, 0, 0, 0, 0, 0, 0, 0, 0), 1, p), x.a, x.b) for x in X]
                self.P[a][i] = [ring.reduce(x) for x in row] if reduce else row
                if a == 4: break
                # X_{a+1} = kinv[a] (zD + l L) X_a
                Xn = []
                for n in range(Nz + 1):
                    t = X[n].mul(NEl(R, {LKEY(1): l}))
                    if n >= 1: t = t.add(X[n - 1].D(ring), ring)
                    t = t.mul(self.kinv[a])
                    if reduce: t = ring.reduce(t)
                    Xn.append(t)
                X = Xn
        self.eta_inv = {(0, 3): R.inv(5), (1, 2): R.inv(5), (2, 1): R.inv(5), (3, 0): R.inv(5), (4, 4): R.inv(5)}

    def Delta(self, i):
        return self.ring.mono({'L': -2, 'c0': 2}, self.R.z(3 * i))

    def Dil(self, i, d):
        """(-1)^d [z^{d-1}] P_0^{(i)} as a Laurent polynomial in L (dict)"""
        el = self.P[0][i][d - 1]                  # stored as c0^a (...) / c0^a: cancel
        assert el.b == 0 and all(k[1] == el.a and k[2:] == (0,) * 8 for k in el.num), "P_0 not a pure L-Laurent polynomial"
        s = (-1) ** d
        return {k[0]: v * s % self.R.p for k, v in el.num.items()}

    def pair(self, i, j, Nmax=None):
        """Sum_{ab} P_a^{(i)}(z) eta^{ab} P_b^{(j)}(-z) as a list over z-order of NEl (reduced)"""
        ring = self.ring; R = self.R
        Nz = self.Nz if Nmax is None else Nmax
        out = [ring.const(0) for _ in range(Nz + 1)]
        for (a, b), e in self.eta_inv.items():
            for n in range(Nz + 1):
                for m in range(n + 1):
                    t = self.P[a][i][m].mul(self.P[b][j][n - m]).scale(e * (1 if (n - m) % 2 == 0 else -1))
                    out[n] = out[n].add(t, ring)
        return [ring.reduce(x) for x in out]

    def q0_value(self, el):
        return el.evalpt(self.ring.q0_values(), self.ring)

    def rows_times(self, h):
        """multiply the rows by a z-series h (list of coefficients) or colour-dependent dict i -> list; returns P' """
        ring = self.ring
        out = [[None] * 5 for _ in range(5)]
        for a in range(5):
            for i in range(5):
                if self.P[a][i] is None: continue
                hs = h[i] if isinstance(h, dict) else h
                row = []
                for n in range(self.Nz + 1):
                    acc = ring.const(0)
                    for m in range(n + 1):
                        if m < len(hs) and hs[m]: acc = acc.add(self.P[a][i][n - m].scale(hs[m]), ring)
                    row.append(acc)
                out[a][i] = row
        return out


# ------------------------------------------------------------------ Tier-1: exact q-series rows from the I-function
class QZ(object):
    """q-series (to O(q^(D+1))) with Laurent-in-z Ser coefficients: list over d"""
    def __init__(self, R, lst): self.R = R; self.c = lst
    def mul(self, o):
        R = self.R; D = len(self.c) - 1
        out = []
        for d in range(D + 1):
            acc = None
            for e in range(d + 1):
                t = self.c[e].mul(o.c[d - e])
                acc = t if acc is None else acc.add(t)
            out.append(acc)
        return QZ(R, out)
    def mulq(self, qs):
        """multiply by a z-free q-series (Ser in q)"""
        R = self.R; D = len(self.c) - 1
        out = []
        for d in range(D + 1):
            acc = None
            for e in range(d + 1):
                t = self.c[e].scale(qs.coef(d - e))
                acc = t if acc is None else acc.add(t)
            out.append(acc)
        return QZ(R, out)
    def zD_plus(self, c):
        """(zD + c) with D = q d/dq: on q^d: (d z + c)"""
        R = self.R
        from fq_ifun import lin
        return QZ(R, [f.mul(lin(R, c, d, f.prec + 2)) for d, f in enumerate(self.c)])
    def mulz(self, zser):
        return QZ(self.R, [f.mul(zser) for f in self.c])
    def coef_series(self, n):
        """[z^n] as a q-Ser"""
        R = self.R
        return Ser.from_list(R, [f.coef(n) for f in self.c], 0, len(self.c))


def tier1_rows(R, G, Nz=6, zprec=24):
    """P_a^{(i)} as QZ from the I-function at the zeta point (Genus0 object G to q^D):
       f_0 = F_i/C0, f_{a+1} = (1/(C1 kappa_a)) (zD + zeta^i) f_a,  P_a = (C0/L) h_{wt_i}(z) e^{-zeta^i utilde/z} f_a,
       utilde = Int (L-1) dq/q.  Returns dict (a, i) -> QZ."""
    D = G.D; N = D + 1
    lams = [R.z(k) for k in range(5)]
    ut = Ser(R, 0, [R.zero] + [R.mul(G.L.coef(d), R.inv(d)) for d in range(1, N)])     # utilde
    C0inv = G.C0.inv(); C1inv = G.C1.inv()
    kinv = [G.kappa[a].inv().mul(C1inv) for a in range(4)]
    out = {}
    for i in range(5):
        F = QZ(R, [G.IF.c_zero(d, i, zprec) for d in range(N)])
        f = F.mulq(C0inv)
        # e^{-zeta^i utilde / z}
        E = QZ(R, [Ser.const(R, 1, zprec)] + [Ser.zero(R, zprec) for _ in range(D)])
        term = QZ(R, [Ser.const(R, 1, zprec)] + [Ser.zero(R, zprec) for _ in range(D)])
        x = QZ(R, [Ser.const(R, R.mul((-lams[i]) % R.p, ut.coef(d)), zprec).shift(-1) for d in range(N)])   # -zeta^i utilde/z
        for n in range(1, D + 1):
            term = term.mul(x).mulq(Ser.const(R, Fr(1, n), N))
            E = QZ(R, [E.c[d].add(term.c[d]) for d in range(N)])
        pref = E.mulq(G.C0.mul(G.L.inv())).mulz(Ser.from_list(R, hv(R, wt(R, lams, i), Nz), 0, zprec))
        for a in range(5):
            out[(a, i)] = pref.mul(f)
            if a < 4: f = f.zD_plus(lams[i]).mulq(kinv[a])
    return out


def edge_coeffs(rows, i, j, nmax):
    """Edge coefficients E^{(ij)}_{mn} (m, n <= nmax) of  E_ij(z,w) = [delta_ij Delta_i - Sum_ab P_a^{(i)}(-z) eta^{ab} P_b^{(j)}(-w)]/(z+w)
    (FRAMEB C Step 1, draft (1.1)):  with G_{mn} = delta_ij Delta_i delta_{m0} delta_{n0} - (-1)^{m+n} Sum_ab eta^{ab} P_{a,m}^{(i)} P_{b,n}^{(j)},
    E_{mn} = Sum_{k=0}^{n} (-1)^k G_{m+1+k, n-k}.  Needs rows to z-order m+n+1.  Returns dict (m,n) -> reduced NEl."""
    ring = rows.ring
    def G(m, n):
        if m > rows.Nz or n > rows.Nz: raise ValueError("z-order too low for the edge")
        acc = ring.const(0)
        for (a, b), e in rows.eta_inv.items():
            acc = acc.add(rows.P[a][i][m].mul(rows.P[b][j][n]).scale(e), ring)
        acc = acc.scale(-((-1) ** (m + n)))
        if i == j and m == 0 and n == 0: acc = acc.add(rows.Delta(i), ring)
        return acc
    out = {}
    for m in range(nmax + 1):
        for n in range(nmax + 1):
            if m + n + 1 > rows.Nz: continue
            acc = ring.const(0)
            for k in range(n + 1):
                t = G(m + 1 + k, n - k)
                acc = acc.add(t if k % 2 == 0 else t.neg(), ring)
            out[(m, n)] = ring.reduce(acc)
    return out
