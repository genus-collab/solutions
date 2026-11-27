# -*- coding: utf-8 -*-
"""gsum_q.py -- Hodge-vertex graph sums of the formal quintic in the ring N (PLAN steps 3-4).

Conventions (FRAMEB_PROOF.md s.0.1, C Steps 1-2; cylib.GraphSum bookkeeping):
  F_g = Sum_Gamma |Aut|^-1 Sum_{colourings c} Sum_{labelled decorations}
        prod_v [ Delta_{c_v}^{g_v-1} (1/m_v!) Int_{Mbar_{g_v, n_v+m_v}} Class_{c_v} prod psi^{a_h} prod psi^{d_k} prod_k Dil_{c_v}(d_k) ]
        prod_e E_{c c'}(a_e, b_e),
  dilaton leaves d_k >= 2 (ordered tuples, weight 1/m!), E from fq_rows.edge_coeffs, Dil_i(d) = (-1)^d [z^(d-1)] P_0^{(i)}.
  A frame is (rows P' = P h, vertex class Class_i): Frame(rows, hfac, classes), hfac = None | list | {i: list} of z-coefficients
  multiplying the rows (fq_rows.Rows.rows_times), classes = {i: weight multiset or None (point class)}; Class_i = Hodge_w.
  Genus-1: the one-leg derivative  d_T F_1 = <H>_1 = Sum_d d N_{1,d} Q^d  (divisor equation), with the leg factor
  (-1)^k [z^k] P_1^{(i)} for the insertion H = H^1 with psi^k at a vertex of colour i (R^-1(psi) H = Sum_i P_1^{(i)}(-psi) e_i).
  Per coloured graph: (1/|Aut|) Sum over ALL labelled decorations (not orbit representatives: for a non-symmetric colouring
  the members of an Aut-orbit have different values; the colour sum is unaffected).
API: Frame; term_value(G, dec, col, frame) -> NEl or None; coloured_sum(G, col, frame) -> NEl;
     genus2_sums(frame, deco) -> {(idx, col): NEl}; genus1_dT_sums(frame) -> {(name, col): NEl}; total(dict) -> NEl.
"""
import itertools, math, os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_ring import Ring, NEl, IDX, NAMES, ZERO_KEY
from fq_rows import Rows, edge_coeffs
from fq_hodge import VertexTable, hv, wt


class RowsView(object):
    """the attributes edge_coeffs needs, with replaced rows"""
    def __init__(self, rows, P):
        self.R = rows.R; self.ring = rows.ring; self.Nz = rows.Nz; self.P = P; self.eta_inv = rows.eta_inv
        self._rows = rows
    def Delta(self, i): return self._rows.Delta(i)


class Frame(object):
    def __init__(self, rows, name, hfac=None, classes=None, nmax=2):
        self.rows = rows; self.R = rows.R; self.ring = rows.ring; self.name = name; self.nmax = nmax
        self.P = rows.P if hfac is None else rows.rows_times(hfac)
        self.view = RowsView(rows, self.P)
        if classes is None: classes = {i: None for i in range(5)}
        self.classes = classes
        self.vt = {i: VertexTable(self.R, classes[i]) for i in range(5)}
        self._edges = {}; self._dil = {}; self._dp = {}; self._leg = {}

    def edge(self, i, j, m, n):
        tab = self._edges.get((i, j))
        if tab is None:
            tab = edge_coeffs(self.view, i, j, self.nmax); self._edges[(i, j)] = tab
        return tab[(m, n)]

    def dil(self, i, d):
        r = self._dil.get((i, d))
        if r is None:
            r = self.P[0][i][d - 1].scale((-1) ** d); self._dil[(i, d)] = r
        return r

    def deltapow(self, i, k):
        r = self._dp.get((i, k))
        if r is None:
            ring = self.ring; R = self.R
            if k >= 0: r = self.rows.Delta(i).pow(k, ring)
            else:
                # Delta_i^-1 = zeta^-3i L^2 / c0^2
                r = NEl(R, {(2, 0, 0, 0, 0, 0, 0, 0, 0, 0): R.z(-3 * i)}, 2, 0).pow(-k, ring)
            self._dp[(i, k)] = r
        return r

    def leg(self, i, a, k):
        """(-1)^k [z^k] P_a^{(i)}: insertion H^a with psi^k at a vertex of colour i"""
        return self.P[a][i][k].scale((-1) ** k)


def term_value(G, dec, col, frame, legs=()):
    """value of one labelled decoration at the colouring col (NEl), or None if the vertex integrals vanish.
    legs: tuple of (vertex, a, k) for primary insertions H^a with psi^k (appended after half-edges and dilaton leaves)."""
    R = frame.R; ring = frame.ring; p = R.p
    V = G['V']; ed = G['ed']; hsAt = G['hsAt']; gv = G['gv']
    scal = 1; facs = []
    hp = [0] * (2 * len(ed))
    legs_at = {}
    for (v, a, k) in legs: legs_at.setdefault(v, []).append((a, k))
    for v in range(V):
        c = col[v]; hps, dils = dec[v]
        lk = tuple(k for (a, k) in legs_at.get(v, []))
        val = frame.vt[c].val(gv[v], tuple(hps) + tuple(dils) + lk)
        if val == 0: return None
        scal = scal * val % p * R.inv(math.factorial(len(dils))) % p
        for r, h in enumerate(hsAt[v]): hp[h] = hps[r]
        if gv[v] != 1: facs.append(frame.deltapow(c, gv[v] - 1))
        for d in dils: facs.append(frame.dil(c, d))
        for (a, k) in legs_at.get(v, []): facs.append(frame.leg(c, a, k))
    for e in range(len(ed)):
        i, j = col[ed[e][0]], col[ed[e][1]]
        facs.append(frame.edge(i, j, hp[2 * e], hp[2 * e + 1]))
    res = ring.const(scal)
    for f in facs: res = res.mul(f)
    return res


def coloured_sum(G, col, frame, legs=()):
    """(1/|Aut|) Sum over all labelled decorations at the colouring col"""
    ring = frame.ring; R = frame.R
    acc = ring.const(0)
    for key, decs in G['orbs'].items():
        for dec in decs:
            t = term_value(G, dec, col, frame, legs)
            if t is not None: acc = acc.add(t, ring)
    return ring.reduce(acc.scale(R.inv(G['aut'])))


def genus2_sums(frame, deco, verbose=None):
    """{(graph idx, colouring): NEl} over all 5^V colourings of the leg-free genus-2 graphs (deco = Decorated(2))"""
    out = {}
    for G in deco.graphs:
        for col in itertools.product(range(5), repeat=G['V']):
            out[(G['idx'], col)] = coloured_sum(G, col, frame)
        if verbose: verbose("   graph %d done" % G['idx'])
    return out


# ---- genus 1: one-leg graphs for d_T F_1 = <H>_1 (built by hand; StableGraph-1.bin does not exist)
def _dils_options(g, hs, nleg_k):
    """ordered dilaton tuples (parts >= 2) and leg exponents with  k + sum d <= dim = 3g-3+hs+1+m"""
    opts = []
    for m in range(0, 4):
        dim = 3 * g - 3 + hs + 1 + m
        if dim < 2 * m: break
        tuples = [()] if m == 0 else [c for s in range(2 * m, dim + 1) for c in _compositions(s, m, 2)]
        for dils in tuples:
            for k in range(0, dim - sum(dils) + 1):
                opts.append((k, dils))
    return opts

def _compositions(s, nd, minp):
    if nd == 0:
        if s == 0: yield ()
        return
    for f in range(minp, s - minp * (nd - 1) + 1):
        for r in _compositions(s - f, nd - 1, minp):
            yield (f,) + r

def genus1_graphs():
    """graph A: one vertex g=1, leg; graph B: one vertex g=0 with a loop, leg.  Same dict format as Decorated.graphs,
    orbs holding ALL labelled decorations under a single key (per-colouring sums iterate over all of them)."""
    A = dict(idx='A', gv=(1,), ed=[], H=[], hsAt=[[]], V=1, aut=1, bound=0)
    A['orbs'] = {('all',): [(((), dils),) for (k, dils) in _dils_options(1, 0, None)]}
    A['legk'] = {(): [k for (k, dils) in _dils_options(1, 0, None)]}
    B = dict(idx='B', gv=(0,), ed=[(0, 0)], H=[(0, 0), (0, 1)], hsAt=[[0, 1]], V=1, aut=2, bound=0)
    B['orbs'] = {('all',): [(((0, 0), ()),)]}
    return A, B

def genus1_dT_sums(frame):
    """{('A' or 'B', (i,)): NEl}: the one-leg genus-1 graph sum for <H>_1 = d_T F_1, per colour"""
    ring = frame.ring; R = frame.R
    A, B = genus1_graphs()
    out = {}
    for i in range(5):
        acc = ring.const(0)
        for (k, dils) in _dils_options(1, 0, None):
            t = term_value(A, (((), dils),), (i,), frame, legs=((0, 1, k),))
            if t is not None: acc = acc.add(t, ring)
        out[('A', (i,))] = ring.reduce(acc)
        t = term_value(B, (((0, 0), ()),), (i,), frame, legs=((0, 1, 0),))
        out[('B', (i,))] = ring.reduce(t.scale(R.inv(2))) if t is not None else ring.const(0)
    return out


def total(d, ring):
    acc = ring.const(0)
    for v in d.values(): acc = acc.add(v, ring)
    return ring.reduce(acc)


# ---- standard frames
def make_frames(rows, which=None):
    """dict name -> Frame.  R := the Teleman rows (QRR constants h_{wt_i} included, fq_rows.Rows);  R~ := R h_{wt_i}^-1.
    F0  (R~, Hodge_wt)   reference (true vertex class)        F0b (R, Hodge_wt)     control (constants counted twice)
    F1  (R h_(1)^-1, Hodge_(1))  target                        F2  (R h_(1)^+1, Hodge_(1))  control
    F3  (R h_(5)^-1, Hodge_(5))                                F4  (R h_(zeta^i)^-1, Hodge_(zeta^i))  colour-dependent
    F5  (R, 1) point frame                                     F5b (R~, 1) control (no Hodge class with the untwisted rows)
    F4b (R h_(zeta^i)^+1, Hodge_(zeta^i))  wrong-sign control with colour-dependent (lambda-covariant) weights.
    NOTE (zeta point, genus 2): for a colour-independent, lambda-independent h and vertex class (F1, F2, F3) every decorated
    coloured term has lambda-degree n = (total z-order of the h factors) + (extra Hodge degree) <= 3g-3+... < 5, and the Z/5
    colour sum multiplies it by Sum_i zeta^{in} = 0 unless n = 0; so the TOTALS of F1, F2, F3 equal F0 for a structural reason
    whatever the sign of the exponent.  The genuine tests of the frame-change identity are F0 = F5 and F4 (= F0), with F4b
    (wrong sign) as the control that can fail.
    Hodge_w = prod_a Lambda^vee(w_a)/w_a^g (FRAMEB A3.2), so Hodge_(1) = Lambda^vee(1), Hodge_(5) = Lambda^vee(5)/5^g."""
    R = rows.R; Nz = rows.Nz
    lams = [R.z(k) for k in range(5)]
    wts = {i: wt(R, lams, i) for i in range(5)}
    hwt_inv = {i: hv(R, wts[i], Nz, -1) for i in range(5)}
    specs = {
        'F0': (hwt_inv, wts),
        'F0b': (None, wts),
        'F1': (hv(R, [1], Nz, -1), {i: [1] for i in range(5)}),
        'F2': (hv(R, [1], Nz, +1), {i: [1] for i in range(5)}),
        'F3': (hv(R, [5], Nz, -1), {i: [5] for i in range(5)}),
        'F4': ({i: hv(R, [lams[i]], Nz, -1) for i in range(5)}, {i: [lams[i]] for i in range(5)}),
        'F4b': ({i: hv(R, [lams[i]], Nz, +1) for i in range(5)}, {i: [lams[i]] for i in range(5)}),   # wrong sign, covariant
        'F5': (None, None),
        'F5b': (hwt_inv, None),
    }
    names = which if which is not None else list(specs)
    return {n: Frame(rows, n, specs[n][0], specs[n][1]) for n in names}


# ---- substitution of series images, split by the L-power (L is handled separately by the callers)
from fq_series import Ser

def subst_split(el, ring, img):
    """el: NEl; img: name -> Ser (c0..c3, i0..i3, Lam; all in the same variable x, L NOT used).
    Returns dict Lpow -> Ser: the L-free numerator parts, each divided by the images of c0^a W^b."""
    R = el.R
    names = [n for n in NAMES if n != 'L']
    gens = [img[n] for n in names]
    prec = min(g.prec for g in gens)
    cache = {}
    def genpow(j, e):
        key = (j, e); r = cache.get(key)
        if r is None:
            r = gens[j] if e == 1 else genpow(j, e - 1).mul(gens[j])
            cache[key] = r
        return r
    pcache = {}
    def prod(key):
        if len(key) == 0: return None
        if key in pcache: return pcache[key]
        r = prod(key[:-1]); e = key[-1]
        if e:
            g = genpow(len(key) - 1, e)
            r = g if r is None else r.mul(g)
        pcache[key] = r
        return r
    groups = {}
    for k, v in el.num.items():
        pr = prod(k[1:])
        t = Ser.const(R, v, prec) if pr is None else pr.scale(v)
        acc = groups.get(k[0])
        groups[k[0]] = t if acc is None else acc.add(t)
    den = None
    if el.a:
        den = img['c0'].inv().pow(el.a)
    if el.b:
        W = img['c0'].mul(img['c0']).add(img['c0'].mul(img['i1'])).sub(img['c1'].mul(img['i0']))
        Wi = W.inv().pow(el.b)
        den = Wi if den is None else den.mul(Wi)
    if den is not None:
        groups = {k: v.mul(den) for k, v in groups.items()}
    return groups


def to_qseries(el, ring, img):
    """q-series of el from the images img (fq_ring.Ring.images_from_qseries: includes 'L')"""
    R = el.R
    groups = subst_split(el, ring, img)
    L = img['L']; Linv = None
    acc = None
    for k, U in groups.items():
        if k == 0: t = U
        elif k > 0: t = U.mul(L.pow(k))
        else:
            if Linv is None: Linv = L.inv()
            t = U.mul(Linv.pow(-k))
        acc = t if acc is None else acc.add(t)
    return acc if acc is not None else Ser.zero(R, min(g.prec for g in img.values()))


def u_to_s(U):
    """u-series -> s-series via u = s^5 (exact spreading of coefficients)"""
    R = U.R
    if not U.c: return Ser(R, 5 * U.prec, [])
    c = []
    for x in U.c: c.extend([x, R.zero, R.zero, R.zero, R.zero])
    return Ser(R, 5 * U.v, c)                       # precision O(s^(5 prec))


def to_sseries(el, ring, img_u):
    """Cf image as a Laurent series in s: Sum_k s^-k U_k(s^5), img_u = images of c, i (u-series over R), L -> 1/s"""
    R = el.R
    groups = subst_split(el, ring, img_u)
    acc = None
    for k, U in groups.items():
        t = u_to_s(U).shift(-k)
        acc = t if acc is None else acc.add(t)
    if acc is None:
        prec = min(g.prec for g in img_u.values())
        return Ser(R, 5 * prec, [])
    return acc


def cf_images_u(C, R, Lam=0):
    """u-series images (over R) of c_j, i_j (j = 0..3), Lam, from a fq_conifold.Conifold object C"""
    def conv(ser): return Ser(R, ser.v, [R.of(x) for x in ser.c])
    j0 = [conv(x) for x in C.j0]; j1 = [conv(x) for x in C.j1]
    N = min(x.prec for x in j0 + j1)
    img = {}
    for j in range(4):
        img['c%d' % j] = j0[j].truncate(N)
        t = j1[j].sub(j0[j].scale(Lam))
        if j >= 1: t = t.sub(j0[j - 1].scale(j))
        img['i%d' % j] = t.truncate(N)
    img['Lam'] = Ser.const(R, Lam, N)
    return img


def is_u_series(S):
    """all s-exponents of the Laurent series S divisible by 5"""
    return all(R_iszero(S.R, x) or (S.v + i) % 5 == 0 for i, x in enumerate(S.c))

def R_iszero(R, x): return R.iszero(x)


def s_to_u(S):
    """s-Laurent series with exponents = 0 mod 5 -> u-series (u = s^5)"""
    R = S.R
    assert is_u_series(S)
    v = S.v
    lo = -((-v) // 5)                  # ceil(v/5)
    c = [S.coef(5 * m) for m in range(lo, S.prec // 5 + (1 if S.prec % 5 else 0)) if 5 * m < S.prec]
    return Ser(R, lo, c)


def polar_in_tc(Fu, tc_u, nmax=2):
    """Fu: Laurent series in u;  tc_u: t_c(u) = u + O(u^2) (Ser over the same R).  Returns the Ser in t = t_c of Fu re-expanded
    (exact to the common precision): Fu(u(t)), u(t) the compositional inverse of t_c."""
    R = Fu.R
    ut = tc_u.reverse()                                     # u(t) = t + ...
    v = ut.shift(-1)                                        # v(t) = u/t, v(0) = 1
    m = -Fu.v if Fu.v < 0 else 0
    G = Fu.shift(m)                                         # power series: Fu = u^-m G
    Gt = G.compose(ut)
    res = Gt
    if m: res = Gt.mul(v.inv().pow(m)).shift(-m)
    return res


def ser_fmt(S, R, kmax=None, kmin=None):
    terms = []
    for i, x in enumerate(S.c):
        k = S.v + i
        if R.iszero(x): continue
        if kmin is not None and k < kmin: continue
        if kmax is not None and k > kmax: break
        terms.append("(%s) s^%d" % (R.fmt(x), k))
    return (" + ".join(terms) if terms else "0") + " + O(s^%d)" % S.prec
