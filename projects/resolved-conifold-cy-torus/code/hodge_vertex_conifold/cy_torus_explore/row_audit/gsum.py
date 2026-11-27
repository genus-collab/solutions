# Graph sums at generic weights near Q = 1. Evidence only.
# Usage: py -I gsum.py <genus> <sample index or a0,a1,a2,e1,e2> [frame: B|PT|WRONG] [Lextra]
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cylib import *
import hodge_engine as he
from fractions import Fraction as Fr
from itertools import permutations, product
from functools import lru_cache

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')


def load_graphs(g):
    jf = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'stable_graphs_%d.json' % g)
    if os.path.exists(jf):
        return [(a, [0] * len(a), b) for a, b in json.load(open(jf))]
    raw = open(os.path.join(BASE, 'StableGraph-%d.bin' % g), 'rb').read()
    raw = [x for x in raw[1:-1] if x != 13]
    out = []
    i = 0
    while i < len(raw):
        K = raw[i]
        gv = raw[i + 1:i + 1 + K]
        mv = raw[i + 1 + K:i + 1 + 2 * K]
        lv = raw[i + 1 + 2 * K:i + 1 + 3 * K]
        i2 = i + 1 + 3 * K
        A = [[0] * K for _ in range(K)]
        for j in range(K):
            A[j][j] = lv[j]
        for j in range(K):
            for q in range(j + 1, K):
                A[j][q] = A[q][j] = raw[i2]
                i2 += 1
        out.append((list(gv), list(mv), A))
        i = i2
    return out


def aut_count(gv, A):
    V = len(gv)
    cnt = 0
    for p in permutations(range(V)):
        if all(gv[p[i]] == gv[i] for i in range(V)) and all(A[p[i]][p[j]] == A[i][j] for i in range(V) for j in range(V)):
            cnt += 1
    fac = 1
    from math import factorial
    for i in range(V):
        fac *= 2 ** A[i][i] * factorial(A[i][i])
        for j in range(i + 1, V):
            fac *= factorial(A[i][j])
    return cnt * fac


def edges_of(A):
    V = len(A)
    ed = []
    for i in range(V):
        for _ in range(A[i][i]):
            ed.append((i, i))
    for i in range(V):
        for j in range(i + 1, V):
            for _ in range(A[i][j]):
                ed.append((i, j))
    return ed


@lru_cache(maxsize=None)
def ch_data(g, ds):
    """list of (chs, rational integral) for int prod ch prod psi^ds"""
    n = len(ds)
    D = 3 * g - 3 + n - sum(ds)
    if D < 0:
        return ()
    res = []
    for chs in he.ch_partitions(D):
        v = he.hodge(g, tuple(chs), ds)
        if v != 0:
            res.append((tuple(chs), v))
    return tuple(res)


class Omega:
    """Omega_v = prod Lambda^vee(v_a)/v_a^g; cls None = 1"""

    def __init__(self, cls):
        self.cls = cls
        self.cache = {}
        self.pk = {}

    def p(self, c):
        # c odd degree of ch_c: coefficient -(c-1)! Sum_a v_a^-c
        if c not in self.pk:
            s = 0
            for v in self.cls:
                s += inv(pow(v, c, P))
            from math import factorial
            self.pk[c] = (-factorial(c - 1) % P) * s % P
        return self.pk[c]

    def integral(self, g, ds):
        ds = tuple(sorted(ds))
        key = (g, ds)
        if key in self.cache:
            return self.cache[key]
        n = len(ds)
        if 2 * g - 2 + n <= 0:
            return 0
        if self.cls is None:
            r = fp(he.wk(g, ds)) if sum(ds) == 3 * g - 3 + n else 0
        else:
            r = 0
            from math import factorial
            for chs, val in ch_data(g, ds):
                co = 1
                from collections import Counter
                for c, m in Counter(chs).items():
                    co = co * pow(self.p(c), m, P) % P * inv(factorial(m)) % P
                r += co * fp(val)
            r %= P
        self.cache[key] = r
        return r


class Sum:
    def __init__(self, w, g, frame, Lx=6):
        self.w = w
        self.g = g
        E = w.Ea
        Ee = w.Ee
        self.Mz = 3 * g - 3 + 1
        N = self.Mz + 2
        L = 2 * g - 2 + Lx + 4
        self.L = L
        self.R = build_rows(w, N, L)
        if frame == 'B':
            cls = [Ee, (-w.e1) % P, (-w.e2) % P]
        elif frame == 'BA':       # class (Ea,-e1,-e2): the off-CY singular-row constant
            cls = [E, (-w.e1) % P, (-w.e2) % P]
        elif frame == 'WRONG':
            cls = [E, (-w.a1) % P, (-w.A) % P]
        else:
            cls = None
        self.cls = cls
        self.Om = Omega(cls)
        self.rows, _ = framed_rows(self.R, cls)
        M = self.Mz
        self.E = {(i, j): edge_table(self.R, self.rows, i, j, M) for i in (1, 2) for j in (1, 2)}
        self.Dl = {1: self.R.D1, 2: self.R.D2}
        vD, sD = self.R.D1
        self.Dinv = {1: (-vD, sinv(sD, L)), 2: (0, sinv(self.R.D2[1], L))}
        self.vcache = {}

    def dil(self, c, d):
        return lscal((-1) ** d % P, self.rows[c][0][d - 1])

    def dil_sum(self, c, g, a):
        """Sum_m 1/m! Sum_{d_i >= 2} Int Omega psi^a prod psi^d Dil(d): Laurent"""
        L = self.L
        n = len(a)
        total = (0, [0] * L)
        dim0 = 3 * g - 3 + n - sum(a)

        def rec(prev, ds, cur, multfac):
            nonlocal total
            # current term
            val = self.Om.integral(g, tuple(a) + tuple(ds))
            if val:
                total = ladd(total, lscal(val * multfac % P, cur), L)
            # add further dilaton degree >= prev
            m = len(ds)
            # dimension: sum(d) - m <= dim0 - ... next: need sum(ds)+d <= dim0 + m + 1
            for d in range(prev, dim0 + m + 1 - sum(ds) + 1):
                if d - 1 > self.R.N:
                    break
                # multiplicity factor for multiset: 1/(mult!)
                mu = ds.count(d) + 1
                rec(d, ds + [d], lmul(cur, self.dil(c, d), L), multfac * inv(mu) % P)
        one = (0, const(1, L))
        rec(2, [], one, 1)
        return total

    def vertex(self, c, g, a):
        a = tuple(sorted(a))
        key = (c, g, a)
        if key in self.vcache:
            return self.vcache[key]
        L = self.L
        s = self.dil_sum(c, g, a)
        if g >= 2:
            D = self.Dl[c]
            for _ in range(g - 1):
                s = lmul(s, D, L)
        elif g == 0:
            s = lmul(s, self.Dinv[c], L)
        self.vcache[key] = s
        return s

    def graph_colour(self, G, col):
        """sum over half-edge powers at a fixed colouring: Laurent"""
        gv, mv, A = G
        V = len(gv)
        ed = edges_of(A)
        L = self.L
        H = []
        at = [[] for _ in range(V)]
        for e, (x, y) in enumerate(ed):
            at[x].append((e, 0))
            at[y].append((e, 1))
        dims = [3 * gv[v] - 3 + len(at[v]) for v in range(V)]
        pw = {}  # (e,end)->power
        total = (0, [0] * L)
        # order vertices so edges close early
        order = list(range(V))

        def assign(vidx, hidx, rem, cur):
            nonlocal total
            if vidx == V:
                total = ladd(total, cur, L)
                return
            v = order[vidx]
            if hidx == 0:
                rem = dims[v]
            if hidx == len(at[v]):
                val = self.vertex(col[v], gv[v], [pw[h] for h in at[v]])
                if not any(val[1]):
                    return
                nxt = lmul(cur, val, L)
                assign(vidx + 1, 0, 0, nxt)
                return
            h = at[v][hidx]
            for p in range(rem + 1):
                pw[h] = p
                c2 = cur
                e = h[0]
                if h[1] == 1:
                    x, y = ed[e]
                    ee = self.E[(col[x], col[y])]
                    c2 = lmul(cur, ee[(pw[(e, 0)], p)], L)
                    if not any(c2[1]):
                        continue
                assign(vidx, hidx + 1, rem - p, c2)
        assign(0, 0, 0, (0, const(1, L)))
        return total

    def run(self, verbose=True):
        g = self.g
        L = self.L
        graphs = load_graphs(g)
        tot = (0, [0] * L)
        rows = []
        for idx, G in enumerate(graphs):
            gv, mv, A = G
            V = len(gv)
            au = aut_count(gv, A)
            iau = inv(au)
            gtot = (0, [0] * L)
            percol = []
            for col in product((1, 2), repeat=V):
                t = self.graph_colour(G, col)
                t = lscal(iau, t)
                gtot = ladd(gtot, t, L)
                percol.append((col, t))
            tot = ladd(tot, gtot, L)
            rows.append((idx, gv, A, percol, gtot))
        return tot, rows


def predicted(w, g, L):
    """P_g Li_{3-2g}(e^t)/(e1 e2) as Laurent in t"""
    from cylib import tau_series
    tau = tau_series(L + 2 * g + 4)
    reg = sinv(tau[1:] + [0], L + 2 * g + 4)
    A = (-1, [(-x) % P for x in reg])  # -1/tau
    for _ in range(2 * g - 3):
        A = lder(A)
    # P_g
    Ln = 2 * g + 2
    S = lambda e: [(0 if k % 2 else _Sc(k // 2, e)) for k in range(Ln)]
    a = S(w.e1)
    b = S(w.e2)
    Pg = 0
    for i in range(2 * g + 1):
        Pg += a[i] * b[2 * g - i]
    Pg %= P
    return lscal(Pg * inv(w.e1 * w.e2 % P) % P, A), Pg


def _Sc(i, e):
    # (y/2)/sinh(y/2) = Sum_i c_i y^(2i), c_i = -(2^(2i)-2)/2^(2i) B_{2i}/(2i)! ... computed from BERN
    from math import factorial
    c = Fr(2 - 2 ** (2 * i), 2 ** (2 * i)) * BERN[2 * i] / factorial(2 * i)
    return fp(c) * pow(e, 2 * i, P) % P


def constant_term(w, g):
    """Int_{M_g} [e(wt_+)^(g-1) Hodge_{wt_+} + e(wt_-)^(g-1) Hodge_{wt_-}]"""
    wp = [w.a0, w.a1, w.a2, (-w.e1) % P, (-w.e2) % P]
    wm = [(-w.a0) % P, w.B, w.A, (-w.e1) % P, (-w.e2) % P]
    r = 0
    for wt in (wp, wm):
        Om = Omega(wt)
        e = 1
        for x in wt:
            e = e * x % P
        r += pow(e, g - 1, P) * Om.integral(g, ())
    return r % P


if __name__ == '__main__':
    pass
