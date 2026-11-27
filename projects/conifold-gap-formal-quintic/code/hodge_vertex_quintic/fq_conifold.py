# -*- coding: utf-8 -*-
"""fq_conifold.py -- Frobenius basis of the quintic Picard-Fuchs equation at the conifold point u := 1 - 3125 q = 0
and the conifold frame map Cf (PLAN step 1, 'Conifold frame Cf').

In the variable u:  D = q d/dq = -(1-u) d/du,  the PF operator  P = p_4 D^4 + p_3 D^3 + p_2 D^2 + p_1 D + p_0
with p_4 = u, p_3 = -2(1-u), p_2 = -7/5 (1-u), p_1 = -2/5 (1-u), p_0 = -24/625 (1-u)   (p_k = 3125 q * (...) of the
large-radius form (1-3125q) D^4 - 6250q D^3 - 4375q D^2 - 1250q D - 120q).
Exponents at u = 0: {0, 1, 1, 2} (indicial polynomial P(n) = [u^(n-3)] P(u^n), checked = n (n-1)^2 (n-2)).
Analytic solutions: f_k = u^k + O(u^3), k = 0, 1, 2 (three-dimensional); log solution  g log u + h  with
g = f_1 + gamma f_2 (gamma fixed by the vanishing of the u^-1 obstruction; h analytic with zero u^0, u^1, u^2
coefficients).  g is the vanishing period varpi^cf_1 = u + gamma u^2 + ....
Pairing: the Lagrange concomitant B (fq_ring.Ring.concomitant, same differential polynomial, p_k in u) is
D-constant on pairs of solutions; B(f_a, f_b) are recorded.  varpi^cf_0 := f_0 + alpha g + beta f_2.  Since g is the
vanishing period, B(f_a, g) = 0 for every analytic f_a (checked: .pairing_vacuous), so the requirement
B(varpi^cf_0, varpi^cf_1) = 0 = B(varpi_0, varpi_1) holds automatically and fixes nothing: (alpha, beta) and the
scale of g remain free (defaults alpha = beta = 0, g = u + O(u^2)).   t_c := g/varpi^cf_0.
Cf: c_j -> D^j varpi^cf_0,  i_j -> D^j varpi^cf_1 - j D^(j-1) varpi^cf_0 - Lam D^j varpi^cf_0  (from
I_1 = C0 Lam + Itilde_1, D^j I_1 = c_j Lam + j c_(j-1) + i_j),  L -> 1/s with u = s^5, Lam -> a chosen value
(elements that are Lam-free in N give Lam-independent images).  Images are Laurent series in s (Ser over R).

Conifold(R, Nu=40, alpha=0, beta=0): .f (list of the three analytic Ser in u, QRing exact), .h (log companion), .Bmat (pairing
constants), .beta, .varpi0 (u-series), .tc (u-series), .cf_images(Lam) -> dict for NEl.subst, .Ds(ser) the
operator D on s-series, .to_s(user) u-series -> s-series.
"""
import os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fq_fp import QRing
from fq_series import Ser


class LogPair(object):
    """A log u + B with Laurent Ser A, B; D acts with D log u = -(1-u)/u"""
    def __init__(self, Q, A, B): self.Q = Q; self.A = A; self.B = B
    def D(self, one_m):
        A = self.A; B = self.B
        dA = one_m.mul(A.deriv()).neg()
        dB = one_m.mul(B.deriv()).neg().add(A.mul(one_m).shift(-1).neg())
        return LogPair(self.Q, dA, dB)
    def mulser(self, f): return LogPair(self.Q, self.A.mul(f), self.B.mul(f))
    def add(self, o): return LogPair(self.Q, self.A.add(o.A), self.B.add(o.B))


class Conifold(object):
    def __init__(self, R, Nu=40, alpha=0, beta=0):
        self.R = R; self.Nu = Nu
        Q = QRing(); self.Q = Q
        self.N = Nu + 12
        N = self.N
        self.one_m = Ser.from_list(Q, [1, -1], 0, N + 8)        # 1 - u
        self.pk = {4: Ser.x(Q, N + 8), 3: self.one_m.scale(-2), 2: self.one_m.scale(Fr(-7, 5)),
                   1: self.one_m.scale(Fr(-2, 5)), 0: self.one_m.scale(Fr(-24, 625))}
        # indicial polynomial and the operator on monomials
        self.P = lambda n: Fr(n) * Fr(n - 1) ** 2 * Fr(n - 2)
        for n in range(0, 8):
            Lun = self.Lop(Ser.mono(Q, n, 1, n + 12))
            assert Lun.coef(n - 3) == self.P(n) if n >= 3 else True
            if n < 3: assert all(Lun.coef(k) == 0 for k in range(n - 3, 0))
            assert all(Lun.coef(k) == 0 for k in range(max(n - 3, 0)) )
        # analytic solutions f_k = u^k + O(u^3)
        self.f = [self.solve_analytic(k) for k in range(3)]
        # log solution g log u + h, g = f_1 + gamma f_2 fixed by the vanishing of the u^-1 obstruction
        h1, ob1 = self.solve_log(self.f[1]); h2, ob2 = self.solve_log(self.f[2])
        assert ob1[0] == ob1[1] == ob2[0] == ob2[1] == 0, (ob1, ob2)
        self.gamma = -ob1[2] / ob2[2]
        self.g = self.f[1].add(self.f[2].scale(self.gamma))        # the vanishing period varpi^cf_1 = u + gamma u^2 + ...
        self.h, self.log_obstruction = self.solve_log(self.g)
        assert all(x == 0 for x in self.log_obstruction)
        # pairings in the basis (f_0, g, f_2)
        self.basis = [self.f[0], self.g, self.f[2]]
        self.Bmat = [[None] * 3 for _ in range(3)]
        self.Bconst_ok = True
        for a in range(3):
            for b in range(3):
                Bab = self.concomitant(self.basis[a], self.basis[b])
                c0 = Bab.coef(0)
                if not all(Bab.coef(k) == 0 for k in range(1, Bab.prec)): self.Bconst_ok = False
                self.Bmat[a][b] = c0
        self.Blog = self.concomitant(self.g, LogPair(Q, self.g, self.h))     # B(g, log solution): B-part, constant
        # varpi0 = f0 + alpha f1 + beta f2.  B(f_a, f_1) = 0 for ALL analytic f_a (f_1 is the vanishing period,
        # orthogonal to the monodromy-invariant 3-space), so the pairing normalisation B(varpi0, f1) = 0 is automatic
        # and leaves (alpha, beta) free; defaults 0 (zero u^1, u^2 coefficients).  Recorded by the driver.
        self.pairing_vacuous = all(self.Bmat[a][1] == 0 for a in range(3))
        self.alpha = Fr(alpha); self.beta = Fr(beta)
        self.varpi0 = self.f[0].add(self.g.scale(self.alpha)).add(self.f[2].scale(self.beta))
        self.varpi1 = self.g
        self.tc = self.varpi1.mul(self.varpi0.inv())
        # D-jets in u
        self.j0 = [self.varpi0]; self.j1 = [self.varpi1]
        for j in range(4):
            self.j0.append(self.Dop(self.j0[-1])); self.j1.append(self.Dop(self.j1[-1]))

    # ---- operators on u-series (QRing)
    def Dop(self, f): return self.one_m.mul(f.deriv()).neg()
    def Lop(self, f):
        """P f for a Ser f in u (or a LogPair: returns the B-part of P(A log u + B) - log u P(A))"""
        if isinstance(f, LogPair):
            jets = [f]
            for k in range(4): jets.append(jets[-1].D(self.one_m))
            acc = None
            for k in range(5):
                t = jets[k].B.mul(self.pk[k])
                acc = t if acc is None else acc.add(t)
            return acc
        jets = [f]
        for k in range(4): jets.append(self.Dop(jets[-1]))
        acc = None
        for k in range(5):
            t = jets[k].mul(self.pk[k])
            acc = t if acc is None else acc.add(t)
        return acc

    def solve_analytic(self, k):
        """power series solution with b_k = 1, b_j = 0 for j in {0,1,2}, j != k"""
        Q = self.Q; N = self.N
        b = [Fr(0)] * N; b[k] = Fr(1)
        for n in range(3, N):
            f = Ser.from_list(Q, b[:n], 0, n + 12)
            r = self.Lop(f).coef(n - 3)
            b[n] = -r / self.P(n)
        f = Ser.from_list(Q, b, 0, N)
        chk = self.Lop(f)
        assert all(chk.coef(m) == 0 for m in range(chk.v, N - 5)), "analytic solution check"
        return f

    def solve_log(self, g):
        """h with P(g log u + h) = 0; h_0 = h_1 = h_2 = 0; returns (h, obstruction = u^-3..u^-1 coefficients)"""
        Q = self.Q; N = self.N
        K = self.Lop(LogPair(Q, g, Ser.zero(Q, N)))          # P(g log u) minus log u P(g)
        obst = [K.coef(m) for m in range(-3, 0)]
        b = [Fr(0)] * N
        for n in range(3, N):
            f = Ser.from_list(Q, b[:n], 0, n + 12)
            r = self.Lop(f).coef(n - 3) + K.coef(n - 3)
            b[n] = -r / self.P(n)
        h = Ser.from_list(Q, b, 0, N)
        return h, obst

    def concomitant(self, f, g):
        """B(f,g) = Sum_{k=1}^4 Sum_{j=0}^{k-1} (-1)^j D^j(p_k g) D^{k-1-j} f  (g may be a LogPair: B-part)"""
        fj = [f]
        for k in range(3): fj.append(self.Dop(fj[-1]))
        acc = None
        for k in range(1, 5):
            t = g.mulser(self.pk[k]) if isinstance(g, LogPair) else self.pk[k].mul(g)
            for j in range(k):
                tb = t.B if isinstance(t, LogPair) else t
                term = tb.mul(fj[k - 1 - j])
                if j % 2: term = term.neg()
                acc = term if acc is None else acc.add(term)
                t = t.D(self.one_m) if isinstance(t, LogPair) else self.Dop(t)
        return acc

    # ---- s-series (over R) and Cf
    def to_s(self, user):
        """u-series (QRing) -> Ser in s over R via u = s^5"""
        R = self.R
        c = [R.of(x) for x in user.c]
        f = Ser(R, user.v, c)
        return f.compose(Ser.mono(R, 5, 1, 5 * f.prec))

    def Ds(self, f):
        """D = -((1 - s^5)/5) s^-4 d/ds on an s-series over R"""
        R = self.R
        one_m = Ser.from_list(R, [1, 0, 0, 0, 0, -1], 0, f.prec + 10)
        return one_m.mul(f.deriv().shift(-4)).scale(Fr(-1, 5))

    def cf_images(self, Lam=0):
        """images of the generators of N under Cf, as Laurent Ser in s over R"""
        R = self.R
        img = {}
        c = [self.to_s(x) for x in self.j0]
        i1 = [self.to_s(x) for x in self.j1]
        N = min(x.prec for x in c + i1)
        for j in range(4):
            img['c%d' % j] = c[j]
            t = i1[j].sub(c[j].scale(Lam))
            if j >= 1: t = t.sub(c[j - 1].scale(j))
            img['i%d' % j] = t
        img['L'] = Ser.mono(R, -1, 1, N)
        img['Lam'] = Ser.const(R, Lam, N)
        return img


if __name__ == "__main__":
    from fq_fp import FpRing
    import time
    t0 = time.time()
    C = Conifold(FpRing(), Nu=20)
    print("time", time.time() - t0)
    print("f0:", C.f[0].fmt(4)); print("f1:", C.f[1].fmt(4)); print("f2:", C.f[2].fmt(4)); print("h:", C.h.fmt(5))
    print("log obstruction (u^-3..u^-1):", C.log_obstruction)
    print("B constant:", C.Bconst_ok, "Bmat:", C.Bmat, "pairing vacuous:", C.pairing_vacuous)
    print("gamma:", C.gamma, "g:", C.g.fmt(4)); print("B(g, log sol) B-part:", C.Blog.fmt(3))
    print("varpi0:", C.varpi0.fmt(4)); print("tc:", C.tc.fmt(4))
    img = C.cf_images(0)
    print("Cf(c1) == Ds Cf(c0):", C.Ds(img['c0']).eq(img['c1'], img['c1'].prec - 6))
    print("Cf(i1) == Ds Cf(i0) - Cf(c0) [D Lam = 1]:", C.Ds(img['i0']).sub(img['c0']).eq(img['i1'], img['i1'].prec - 6))
