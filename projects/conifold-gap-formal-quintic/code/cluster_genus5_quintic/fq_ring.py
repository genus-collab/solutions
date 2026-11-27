# -*- coding: utf-8 -*-
"""fq_ring.py -- the differential ring N of the formal quintic (PLAN step 1).

Generators (all functions of q; D = q d/dq; L = (1 - 3125 q)^(-1/5), DL = (L^6 - L)/5, q = (1 - L^-5)/3125):
  c_j := D^j C0 (j = 0..3),  i_j := D^j Itilde_1 (j = 0..3),  Lam := log q (D Lam = 1; physical elements are Lam-free),
  with D^4 reduced by the Picard-Fuchs equation  [(1-3125q) D^4 - 6250q D^3 - 4375q D^2 - 1250q D - 120q] f = 0:
     c_4 = (L^5-1) (2 c_3 + 7/5 c_2 + 2/5 c_1 + 24/625 c_0),
     i_4 = (L^5-1) (2 i_3 + 7/5 i_2 + 2/5 i_1 + 24/625 i_0) - 4 c_3 + (L^5-1)(6 c_2 + 14/5 c_1 + 2/5 c_0)
  (inhomogeneous: I_1 = C0 Lam + Itilde_1 is a solution).  Tracked denominators: c_0 and  W := c_0^2 C1 =
  c_0^2 + c_0 i_1 - c_1 i_0  (so C1 = D T = W/c_0^2, 1/C1 = c_0^2/W).
Element:  NEl(R, num, a, b) = num / (c_0^a W^b),  num = {key: coeff mod p},
  key = (Lpow, e_c0, e_c1, e_c2, e_c3, e_i0, e_i1, e_i2, e_i3, e_Lam)  (Lpow may be negative).
Operations: add, sub, mul, neg, scale, D, pow, iszero (numerator zero), eq; Ldeg() -> (min, max) L-power;
  subst(img) -> Ser (images of the generators: dict with keys 'c0'..'c3','i0'..'i3','L','Lam'); evalpt(vals) -> scalar.
The quadratic Wronskian relation (derived in Ring.derive_relation, antisymmetric Lagrange concomitant of the
self-adjoint PF operator, B(varpi_0, varpi_1) = 0) is linear in i_3:  i_3 = L^5 * Rest / c_0.  Ring.reduce(el)
eliminates i_3 (normal form in F_p[L^+-1, c_0..c_3, i_0..i_2, Lam][c_0^-1, W^-1]).
Ring(R) holds the constants; use Ring.gen(name), Ring.const(x), Ring.Lpow(k), Ring.C1inv(), Ring.W().
"""
import os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fq_series import Ser

NG = 10
IDX = {'L': 0, 'c0': 1, 'c1': 2, 'c2': 3, 'c3': 4, 'i0': 5, 'i1': 6, 'i2': 7, 'i3': 8, 'Lam': 9}
NAMES = ['L', 'c0', 'c1', 'c2', 'c3', 'i0', 'i1', 'i2', 'i3', 'Lam']
ZERO_KEY = (0,) * NG


def _unit(j, k=1):
    e = [0] * NG; e[j] = k; return tuple(e)


class NEl(object):
    __slots__ = ('R', 'num', 'a', 'b')

    def __init__(self, R, num, a=0, b=0):
        self.R = R; self.num = num; self.a = a; self.b = b

    # ---- helpers
    @staticmethod
    def _addto(d, k, v, p):
        x = (d.get(k, 0) + v) % p
        if x: d[k] = x
        elif k in d: del d[k]

    @staticmethod
    def _mulnum(A, B, p):
        out = {}
        for k1, v1 in A.items():
            for k2, v2 in B.items():
                k = tuple(x + y for x, y in zip(k1, k2))
                x = (out.get(k, 0) + v1 * v2) % p
                if x: out[k] = x
                elif k in out: del out[k]
        return out

    def _raise(self, a, b, ring):
        """numerator rewritten over the denominator c_0^a W^b (a >= self.a, b >= self.b)"""
        num = self.num; p = self.R.p
        if a > self.a:
            num = NEl._mulnum(num, {_unit(IDX['c0'], a - self.a): 1}, p)
        for _ in range(b - self.b):
            num = NEl._mulnum(num, ring.Wnum, p)
        return num

    # ---- arithmetic
    def add(self, o, ring):
        a = max(self.a, o.a); b = max(self.b, o.b); p = self.R.p
        n1 = self._raise(a, b, ring); n2 = o._raise(a, b, ring)
        out = dict(n1)
        for k, v in n2.items(): NEl._addto(out, k, v, p)
        return NEl(self.R, out, a, b)

    def neg(self): return NEl(self.R, {k: (-v) % self.R.p for k, v in self.num.items()}, self.a, self.b)
    def sub(self, o, ring): return self.add(o.neg(), ring)
    def scale(self, c):
        c = self.R.of(c); p = self.R.p
        return NEl(self.R, {k: v * c % p for k, v in self.num.items() if v * c % p}, self.a, self.b)
    def mul(self, o):
        return NEl(self.R, NEl._mulnum(self.num, o.num, self.R.p), self.a + o.a, self.b + o.b)
    def pow(self, n, ring):
        r = ring.const(1)
        for _ in range(n): r = r.mul(self)
        return r
    def iszero(self): return not self.num
    def eq(self, o, ring): return self.sub(o, ring).iszero()
    def nterms(self): return len(self.num)
    def Ldeg(self):
        if not self.num: return None
        ks = [k[0] for k in self.num]
        return (min(ks), max(ks))
    def has(self, name):
        j = IDX[name]
        return any(k[j] for k in self.num)
    def maxexp(self, name):
        j = IDX[name]
        return max((k[j] for k in self.num), default=0)

    def D(self, ring):
        """derivation D = q d/dq (see module docstring); result over c_0^(a+1) W^(b+1) when a or b > 0"""
        R = self.R; p = R.p
        dn = ring._Dnum(self.num)
        if self.a == 0 and self.b == 0: return NEl(R, dn, 0, 0)
        a2 = self.a + (1 if self.a else 0); b2 = self.b + (1 if self.b else 0)
        out = {}
        t = dn
        if self.a: t = NEl._mulnum(t, {_unit(IDX['c0']): 1}, p)
        if self.b: t = NEl._mulnum(t, ring.Wnum, p)
        for k, v in t.items(): NEl._addto(out, k, v, p)
        if self.a:
            t = NEl._mulnum(self.num, {_unit(IDX['c1']): (-self.a) % p}, p)
            if self.b: t = NEl._mulnum(t, ring.Wnum, p)
            for k, v in t.items(): NEl._addto(out, k, v, p)
        if self.b:
            t = NEl._mulnum(self.num, ring.DWnum, p)
            t = NEl._mulnum(t, {ZERO_KEY: (-self.b) % p}, p)
            if self.a: t = NEl._mulnum(t, {_unit(IDX['c0']): 1}, p)
            for k, v in t.items(): NEl._addto(out, k, v, p)
        return NEl(R, out, a2, b2)

    def subst(self, img, ring):
        """substitute Ser images for the generators (img: name -> Ser; 'c0'.. 'i3', 'L', 'Lam'); returns Ser"""
        R = self.R
        gens = [img[n] for n in NAMES]
        prec = min(g.prec for g in gens if g is not None)
        cache = {}
        def genpow(j, e):                      # e != 0; returns Ser
            key = (j, e)
            r = cache.get(key)
            if r is None:
                if e > 0: r = genpow(j, e - 1).mul(gens[j]) if e > 1 else gens[j]
                else:
                    inv = cache.get((j, -1))
                    if inv is None:
                        inv = gens[j].inv(); cache[(j, -1)] = inv
                    r = inv if e == -1 else genpow(j, e + 1).mul(inv)
                cache[key] = r
            return r
        pcache = {}
        def prod(key):                         # None stands for 1
            if len(key) == 0: return None
            if key in pcache: return pcache[key]
            r = prod(key[:-1]); e = key[-1]
            if e:
                g = genpow(len(key) - 1, e)
                r = g if r is None else r.mul(g)
            pcache[key] = r
            return r
        acc = None
        for k, v in self.num.items():
            pr = prod(k)
            t = Ser.const(R, v, prec) if pr is None else pr.scale(v)
            acc = t if acc is None else acc.add(t)
        if acc is None: acc = Ser.zero(R, prec)
        if self.a: acc = acc.mul(img['c0'].inv().pow(self.a))
        if self.b:
            W = ring.W().subst(img, ring)
            acc = acc.mul(W.inv().pow(self.b))
        return acc

    def evalpt(self, vals, ring):
        """evaluate at a point (vals: name -> field element); returns a field element"""
        R = self.R; p = R.p
        v = [vals[n] for n in NAMES]
        tot = 0
        for k, c in self.num.items():
            t = c
            for j in range(NG):
                if k[j]: t = t * R.pow(v[j], k[j]) % p
            tot = (tot + t) % p
        if self.a: tot = tot * R.pow(R.inv(v[IDX['c0']]), self.a) % p
        if self.b:
            W = ring.W().evalpt(vals, ring)
            tot = tot * R.pow(R.inv(W), self.b) % p
        return tot

    def __repr__(self):
        items = sorted(self.num.items())[:6]
        return "NEl(%d terms, /c0^%d W^%d; first: %s)" % (len(self.num), self.a, self.b, items)


class Ring(object):
    """constants of N: D-images of the generators, W, the relation; factory methods"""
    def __init__(self, R):
        self.R = R; p = R.p
        self.Dgen = [None] * NG
        # D c_j = c_{j+1}, D i_j = i_{j+1}, D Lam = 1
        for j in range(3):
            self.Dgen[IDX['c%d' % j]] = {_unit(IDX['c%d' % (j + 1)]): 1}
            self.Dgen[IDX['i%d' % j]] = {_unit(IDX['i%d' % (j + 1)]): 1}
        self.Dgen[IDX['Lam']] = {ZERO_KEY: 1}
        # c_4 = (L^5-1)(2 c3 + 7/5 c2 + 2/5 c1 + 24/625 c0)
        cf = {'c3': Fr(2), 'c2': Fr(7, 5), 'c1': Fr(2, 5), 'c0': Fr(24, 625)}
        d4 = {}
        for n, x in cf.items():
            e = list(_unit(IDX[n])); x = R.of(x)
            e5 = list(e); e5[0] = 5
            NEl._addto(d4, tuple(e5), x, p); NEl._addto(d4, tuple(e), (-x) % p, p)
        self.Dgen[IDX['c3']] = d4
        # i_4 = (L^5-1)(2 i3 + 7/5 i2 + 2/5 i1 + 24/625 i0) - 4 c3 + (L^5-1)(6 c2 + 14/5 c1 + 2/5 c0)
        cf = {'i3': Fr(2), 'i2': Fr(7, 5), 'i1': Fr(2, 5), 'i0': Fr(24, 625), 'c2': Fr(6), 'c1': Fr(14, 5), 'c0': Fr(2, 5)}
        d4 = {}
        for n, x in cf.items():
            e = list(_unit(IDX[n])); x = R.of(x)
            e5 = list(e); e5[0] = 5
            NEl._addto(d4, tuple(e5), x, p); NEl._addto(d4, tuple(e), (-x) % p, p)
        NEl._addto(d4, _unit(IDX['c3']), (-4) % p, p)
        self.Dgen[IDX['i3']] = d4
        # W = c0^2 + c0 i1 - c1 i0 ; D W = 2 c0 c1 + c0 i2 - c2 i0
        self.Wnum = {_unit(IDX['c0'], 2): 1, self._key({'c0': 1, 'i1': 1}): 1, self._key({'c1': 1, 'i0': 1}): p - 1}
        self.DWnum = self._Dnum(self.Wnum)
        assert self.DWnum == {self._key({'c0': 1, 'c1': 1}): 2, self._key({'c0': 1, 'i2': 1}): 1,
                              self._key({'c2': 1, 'i0': 1}): p - 1}
        self.rel = None; self.i3num = None; self.i3_c0pow = None
        self.inv5 = R.inv(5)

    @staticmethod
    def _key(d):
        e = [0] * NG
        for n, k in d.items(): e[IDX[n]] = k
        return tuple(e)

    # ---- factories
    def const(self, x):
        x = self.R.of(x)
        return NEl(self.R, {ZERO_KEY: x} if x else {})
    def gen(self, name, k=1): return NEl(self.R, {_unit(IDX[name], k): 1})
    def Lpow(self, k): return NEl(self.R, {_unit(0, k): 1})
    def mono(self, d, coeff=1):
        """monomial from {name: exponent} ('L' allowed)"""
        return NEl(self.R, {self._key(d): self.R.of(coeff)})
    def W(self): return NEl(self.R, dict(self.Wnum))
    def C1(self): return NEl(self.R, dict(self.Wnum), 2, 0)              # W/c0^2
    def C1inv(self): return NEl(self.R, {_unit(IDX['c0'], 2): 1}, 0, 1)   # c0^2/W
    def C0inv(self): return NEl(self.R, {ZERO_KEY: 1}, 1, 0)
    def q(self):
        """q = (1 - L^-5)/3125"""
        p = self.R.p; i = self.R.inv(3125)
        return NEl(self.R, {ZERO_KEY: i, _unit(0, -5): (-i) % p})

    # ---- derivation on numerators
    def _Dnum(self, num):
        p = self.R.p; out = {}
        for k, v in num.items():
            m = k[0]
            if m:                                    # D L^m = m/5 (L^(m+5) - L^m)
                c = v * m % p * self.inv5 % p
                NEl._addto(out, (m + 5,) + k[1:], c, p); NEl._addto(out, k, (-c) % p, p)
            for j in range(1, NG):
                e = k[j]
                if not e: continue
                base = list(k); base[j] = e - 1; base = tuple(base)
                c = v * e % p
                for k2, v2 in self.Dgen[j].items():
                    kk = tuple(x + y for x, y in zip(base, k2))
                    NEl._addto(out, kk, c * v2, p)
        return out

    # ---- the relation
    def concomitant(self, f_jets, g):
        """B(f,g) = Sum_{k=1}^4 Sum_{j=0}^{k-1} (-1)^j D^j(p_k g) D^{k-1-j} f, the Lagrange concomitant of the
        (self-adjoint) PF operator p_4 D^4 + ... + p_0 with p_4 = L^-5, p_3 = -2(1-L^-5), p_2 = -7/5 (1-L^-5),
        p_1 = -2/5 (1-L^-5);  f_jets[m] = D^m f (m = 0..3), g an element (D applied inside)."""
        one_m = self.const(1).sub(self.Lpow(-5), self)          # 1 - L^-5 = 3125 q
        pk = {4: self.Lpow(-5), 3: one_m.scale(-2), 2: one_m.scale(Fr(-7, 5)), 1: one_m.scale(Fr(-2, 5))}
        B = self.const(0)
        for k in range(1, 5):
            t = pk[k].mul(g)
            for j in range(k):
                term = t.mul(f_jets[k - 1 - j])
                if j % 2: term = term.neg()
                B = B.add(term, self)
                t = t.D(self)
        return B

    def derive_relation(self):
        """Rel := B(varpi_0, varpi_1), varpi_0 = c0, varpi_1 = c0 Lam + i0 (Lam-part stored separately, expected 0).
        Stores the i3-elimination  i3 = i3num / c0^i3_c0pow."""
        R = self.R; p = R.p
        f_jets = [self.gen('c%d' % m) for m in range(4)]
        g = self.gen('c0').mul(self.gen('Lam')).add(self.gen('i0'), self)
        B = self.concomitant(f_jets, g)
        assert B.a == 0 and B.b == 0
        self.rel = B
        self.rel_lam = {k: v for k, v in B.num.items() if k[IDX['Lam']]}
        rel0 = {k: v for k, v in B.num.items() if not k[IDX['Lam']]}
        self.rel0 = NEl(R, rel0)
        ji = IDX['i3']
        A = {k: v for k, v in rel0.items() if k[ji]}
        rest = {k: v for k, v in rel0.items() if not k[ji]}
        assert len(A) == 1, "relation not linear in i3 with a monomial coefficient: %s" % A
        (kA, vA), = A.items()
        assert kA[ji] == 1 and all(kA[j] == 0 for j in range(1, NG) if j not in (ji, IDX['c0']))
        m = kA[0]; n = kA[IDX['c0']]
        c = (-R.inv(vA)) % p                       # i3 = -rest L^-m /(vA c0^n)
        self.i3num = {(k[0] - m,) + k[1:]: v * c % p for k, v in rest.items()}
        self.i3_c0pow = n
        return B

    def reduce(self, el):
        """eliminate i3 with the relation; result over c_0^(a + n*kmax) W^b"""
        if self.i3num is None: self.derive_relation()
        ji = IDX['i3']; p = self.R.p; n = self.i3_c0pow
        kmax = el.maxexp('i3')
        if kmax == 0: return el
        groups = {}
        for k, v in el.num.items():
            e = k[ji]; kk = k[:ji] + (0,) + k[ji + 1:]
            groups.setdefault(e, {})[kk] = v
        powc = {0: {ZERO_KEY: 1}}
        for e in range(1, kmax + 1): powc[e] = NEl._mulnum(powc[e - 1], self.i3num, p)
        out = {}
        for e, g in groups.items():
            t = NEl._mulnum(g, powc[e], p)
            if n * (kmax - e): t = NEl._mulnum(t, {_unit(IDX['c0'], n * (kmax - e)): 1}, p)
            for k, v in t.items(): NEl._addto(out, k, v, p)
        return NEl(self.R, out, el.a + n * kmax, el.b)

    def q0_values(self):
        """generator values at q = 0: L = 1, c0 = 1, others 0 (Lam must be absent)"""
        v = {n: 0 for n in NAMES}; v['L'] = 1; v['c0'] = 1
        return v

    def images_from_qseries(self, G):
        """Ser images (in q) of the generators from a Genus0 object: c_j = D^j C0, i_j = D^j It1, L, Lam = 0"""
        img = {'L': G.L, 'Lam': Ser.zero(self.R, G.N)}
        c = G.C0; i = G.It1
        for j in range(4):
            img['c%d' % j] = c; img['i%d' % j] = i
            c = c.theta(); i = i.theta()
        return img
