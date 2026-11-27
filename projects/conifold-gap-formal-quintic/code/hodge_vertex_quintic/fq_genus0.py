# -*- coding: utf-8 -*-
"""fq_genus0.py -- genus-0 data of the formal quintic from the z -> infinity expansion of the I-function.

Mirror theorem (PLAN Key facts): S(1) = e^{-tau0/z} I/(z C0), S the fundamental solution with
z d_T S(x) = S(H * x) for CONSTANT x.  Rows are built in the FLAT basis H^a (a = 0..4; H^5 = 1 at the zeta
point, since prod_j (H - zeta^j) = H^5 - 1):  reading off the z^0 term of z d_T S(H^a) = S(H*H^a) gives the
structure constants, and at the zeta point Z/5-covariance forces H*H^a = kappa_a H^{a+1 mod 5} (checked), so
   S(H^{a+1}) = z d_T S(H^a) / kappa_a,   d_T = (1/C1) D,  D = q d/dq.
(The naive recursion ((z/C1) D)^k S(1) is NOT S(H*^k) for k >= 3: H*^k depends on T.)
Classical pairing eta(H^a, H^b) = Sum_i lambda_i^{a+b}/W_i^0 = 5 if a+b = 3 mod 5, else 0 (zeta point).

Genus0(R, lams, D, Ns=6):  all q-series to O(q^(D+1)) (Ser in q over R); lams must be the zeta point for the
Z/5 structure (the I-function / tau0 / mirror-map part works for generic weights: see genus0_check.py).
  .C0, .It1   = Sum_d (5d)!/(d!)^5 q^d,  5 Sum_d q^d (5d)!/(d!)^5 Sum_{r=d+1}^{5d} 1/r
  .comp1      H^k-components (Vandermonde) of the log-free s^1 term of I/(z C0); .tau0 = comp1[0];
  .b          = C0 * comp1[1]  (mirror map T = log q + b/C0; b = It1 checked by the driver)
  .Tq, .C1    T = log q + Tq, C1 = D T
  .S[a]       S(H^a) at the fixed points (list over a of lists over i of SL)
  .kappa[a]   H*H^a = kappa_a H^{a+1 mod 5};  .offdiag_ok: other components vanish;  .lfree: z^0 terms log-free
  .closure_ok z d_T S(H^4) = kappa_4 S(1) (checked on the whole SL to the available s-precision)
  .v[k]       cyclic vectors H*^k = (prod_{a<k} kappa_a) H^k as fixed-point vectors, k = 0..5
  .Y(m)       Gram eta(v_k, v_l), k + l = m;  .eta(a,b)
  .mu[i], .idem[i], .Delta[i] = 1/eta(e_i, e_i);  .L = (1 - 3125 q)^(-1/5)
  Qseries(f): re-expand a q-series in Q = q exp(It1/C0).
"""
import math, os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fq_series import Ser
from fq_ifun import IFun, SL, vandermonde, components
from fq_hodge import W0


def quintic_C0(R, D):
    return Ser.from_list(R, [Fr(math.factorial(5 * d), math.factorial(d) ** 5) for d in range(D + 1)], 0, D + 1)


def quintic_It1(R, D):
    lst = [Fr(0)]
    for d in range(1, D + 1):
        a = Fr(math.factorial(5 * d), math.factorial(d) ** 5)
        lst.append(5 * a * sum((Fr(1, r) for r in range(d + 1, 5 * d + 1)), Fr(0)))
    return Ser.from_list(R, lst, 0, D + 1)


class Genus0(object):
    def __init__(self, R, lams, D, Ns=6):
        assert Ns >= 5, 'Ns >= 5 needed for the flat recursion'
        self.R = R; self.lams = [R.of(l) for l in lams]; self.D = D; self.Ns = Ns
        N = D + 1; self.N = N
        self.IF = IFun(R, lams, D)
        self.V, self.Vinv = vandermonde(R, self.lams)
        self.W0 = [W0(R, self.lams, i) for i in range(5)]
        self.C0 = quintic_C0(R, D)
        self.It1 = quintic_It1(R, D)
        C0inv = self.C0.inv()
        # J_i = q^{lam_i s} Sum_d q^d c_d(s) / C0
        J = []
        for i in range(5):
            acc = None
            for d in range(D + 1):
                t = SL.from_sser(R, self.IF.c_inf(d, i, Ns + 1), N).mulq(Ser.mono(R, d, 1, N))
                acc = t if acc is None else acc.add(t)
            J.append(SL.exp_ls(R, self.lams[i], Ns, N).mul(acc, Ns).mulq(C0inv))
        self.J = J
        self.comp1 = components(R, self.Vinv, [J[i].part(1, 0, N) for i in range(5)])
        self.tau0 = self.comp1[0]; self.b = self.comp1[1].mul(self.C0)
        self.comp1l = components(R, self.Vinv, [J[i].part(1, 1, N) for i in range(5)])
        # S(1) = e^{-tau0 s} J
        ex = SL(R, {(0, 0): Ser.const(R, 1, N)})
        pw = Ser.const(R, 1, N)
        for n in range(1, Ns + 1):
            pw = pw.mul(self.tau0.neg())
            ex = ex.add(SL(R, {(n, 0): pw.scale(Fr(1, math.factorial(n)))}))
        self.S1 = [ex.mul(J[i], Ns) for i in range(5)]
        # mirror map
        self.Tq = self.comp1[1]                            # T = log q + Tq
        self.C1 = Ser.const(R, 1, N).add(self.Tq.theta())
        C1inv = self.C1.inv()
        # flat rows S(H^a) and structure constants
        self.S = [self.S1]; self.kappa = []; self.lfree = True; self.offdiag_ok = True
        X = self.S1
        for a in range(5):
            Xd = [x.theta().s_shift(-1).mulq(C1inv).trunc_s(Ns - a - 1) for x in X]   # S(H * H^a)
            vec = []
            for i in range(5):
                sp = Xd[i].s_part(0)
                if any(m >= 1 and not ser.iszero() for m, ser in sp.items()): self.lfree = False
                vec.append(sp.get(0, Ser.zero(R, N)))
            comp = components(R, self.Vinv, vec)
            for c in range(5):
                if c != (a + 1) % 5 and not comp[c].iszero(): self.offdiag_ok = False
            self.kappa.append(comp[(a + 1) % 5])
            if a < 4:
                X = [x.mulq(self.kappa[a].inv()) for x in Xd]
                self.S.append(X)
            else:
                dif = [Xd[i].sub(self.S1[i].mulq(self.kappa[4]).trunc_s(Ns - 5)) for i in range(5)]
                self.closure_ok = all(not d.d for d in dif)
        # cyclic vectors v_k = (prod_{a<k} kappa_a) H^k
        self.v = []
        pr = Ser.const(R, 1, N)
        for k in range(6):
            self.v.append([pr.scale(R.pow(self.lams[i], k)) for i in range(5)])
            pr = pr.mul(self.kappa[k % 5])
        self.gram = {}; self.gram_ok = True
        for k in range(5):
            for l in range(5):
                y = self.eta(self.v[k], self.v[l]); m = k + l
                if m in self.gram:
                    if not self.gram[m].eq(y): self.gram_ok = False
                else: self.gram[m] = y
        # eigenvalues: x^5 = prod kappa_a ;  mu_i = lam_i (prod kappa)^(1/5)
        kp = Ser.const(R, 1, N)
        for a in range(5): kp = kp.mul(self.kappa[a])
        self.kprod = kp
        rt = kp.root(5)
        self.mu = [rt.scale(self.lams[i]) for i in range(5)]
        # idempotents e_i = (1/5) Sum_k mu_i^{-k} H*^k
        self.idem = []; self.Delta = []
        for i in range(5):
            muinv = self.mu[i].inv()
            e = [Ser.zero(R, N) for _ in range(5)]
            pw = Ser.const(R, Fr(1, 5), N)
            for k in range(5):
                for p in range(5): e[p] = e[p].add(pw.mul(self.v[k][p]))
                pw = pw.mul(muinv)
            self.idem.append(e)
            self.Delta.append(self.eta(e, e).inv())
        self.L = Ser.from_list(R, [1, -3125], 0, N).root(5).inv()

    def eta(self, a, b):
        R = self.R
        acc = Ser.zero(R, self.N)
        for i in range(5):
            acc = acc.add(a[i].mul(b[i]).scale(R.inv(self.W0[i])))
        return acc

    def eta_flat(self, a, b):
        """eta(H^a, H^b) (field element)"""
        R = self.R
        return sum((R.mul(R.pow(self.lams[i], a + b), R.inv(self.W0[i])) for i in range(5)), R.zero) % R.p

    def Y(self, m): return self.gram.get(m)

    def Qseries(self, f):
        """re-expand the q-series f in Q = q exp(It1/C0)"""
        R = self.R
        Qq = Ser.x(R, self.N).mul(self.Tq.exp())
        return f.compose(Qq.reverse())
