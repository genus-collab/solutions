# Fresh frame-B graph evaluator (round 6d, backward read). Evidence only.
# Arithmetic in F_p, p = 2^127-1, with omega in F_p. Laurent series in xi = 1+27y with tracked precision.
# Rows: closed form S(v,z) of draft s.1.3 only, then T = z theta + Y(1+z/E). Edges from the normalised pairing.
from fractions import Fraction as Fr
from functools import lru_cache
from itertools import permutations, product
from math import factorial, comb
import hodge

P = 2 ** 127 - 1


def inv(a):
    return pow(a % P, P - 2, P)


def fp(x):
    x = Fr(x)
    return x.numerator % P * inv(x.denominator) % P


# element of order 3
def _omega():
    e = (P - 1) // 3
    for g in range(2, 100):
        w = pow(g, e, P)
        if w != 1:
            assert (w * w + w + 1) % P == 0
            return w


OMEGA = _omega()


def ratrec(a, bound=2 ** 62):
    """rational reconstruction of a mod P"""
    r0, r1 = P, a % P
    s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0:
        return None
    return Fr(r1, s1)


class LS:
    """Laurent series a_{val} x^val + ... known mod x^aprec."""
    __slots__ = ('val', 'c', 'aprec')

    def __init__(self, val, c, aprec):
        # strip
        c = [x % P for x in c]
        c = c[:max(0, aprec - val)]
        i = 0
        while i < len(c) and c[i] == 0:
            i += 1
        self.val = val + i
        self.c = c[i:]
        self.aprec = aprec
        if not self.c:
            self.val = aprec

    @staticmethod
    def const(a, N):
        return LS(0, [a], N)

    def coef(self, k):
        if k >= self.aprec:
            raise ValueError('precision exhausted: asked %d, aprec %d' % (k, self.aprec))
        i = k - self.val
        if i < 0 or i >= len(self.c):
            return 0
        return self.c[i]

    def __add__(self, o):
        if not isinstance(o, LS):
            o = LS(0, [o], self.aprec if self.aprec > 0 else 1)
        ap = min(self.aprec, o.aprec)
        v = min(self.val, o.val)
        if v >= ap:
            return LS(ap, [], ap)
        L = ap - v
        c = [0] * L
        for i, x in enumerate(self.c):
            j = self.val + i - v
            if j < L:
                c[j] = x
        for i, x in enumerate(o.c):
            j = o.val + i - v
            if j < L:
                c[j] = (c[j] + x) % P
        return LS(v, c, ap)

    __radd__ = __add__

    def __neg__(self):
        return LS(self.val, [-x for x in self.c], self.aprec)

    def __sub__(self, o):
        return self + (-o if isinstance(o, LS) else (-o) % P)

    def __mul__(self, o):
        if not isinstance(o, LS):
            o = o % P
            return LS(self.val, [x * o for x in self.c], self.aprec)
        ap = min(self.val + o.aprec, o.val + self.aprec)
        v = self.val + o.val
        if v >= ap:
            return LS(ap, [], ap)
        L = ap - v
        c = [0] * L
        for i, x in enumerate(self.c):
            if x == 0 or i >= L:
                continue
            for j, y in enumerate(o.c):
                if i + j >= L:
                    break
                c[i + j] += x * y
        return LS(v, [t % P for t in c], ap)

    __rmul__ = __mul__

    def inverse(self):
        if not self.c:
            raise ZeroDivisionError
        L = self.aprec - self.val
        a0i = inv(self.c[0])
        b = [0] * L
        b[0] = a0i
        for n in range(1, L):
            s = 0
            for k in range(1, min(n, len(self.c) - 1) + 1):
                s += self.c[k] * b[n - k]
            b[n] = -s % P * a0i % P
        return LS(-self.val, b, -self.val + L)

    def pow(self, n):
        if n < 0:
            return self.inverse().pow(-n)
        r = LS(0, [1], 10 ** 6)
        for _ in range(n):
            r = r * self
        return r

    def deriv(self):
        c = [(self.val + i) * x for i, x in enumerate(self.c)]
        return LS(self.val - 1, c, self.aprec - 1)

    def polar(self):
        return {self.val + i: x for i, x in enumerate(self.c) if self.val + i < 0 and x != 0}

    def isknownzero_polar(self):
        return self.aprec >= 0 and not self.polar()


class Base:
    """Series on the real branch t = -(1-xi)^{1/3}/3 in xi, to absolute precision N."""

    def __init__(self, N):
        self.N = N
        # (1-xi)^{1/3}
        c = [1]
        for k in range(1, N + 3):
            c.append(c[-1] * fp(Fr(1, 3) - (k - 1)) % P * inv(k) % P * (P - 1) % P)
        cube = LS(0, c, N + 3)
        self.t = cube * fp(Fr(-1, 3))
        self.xi = LS(1, [1], N + 3)
        # wc = 2F1(1/3,2/3;1;xi)
        w = [1]
        for k in range(0, N + 3):
            w.append(w[-1] * fp((Fr(1, 3) + k) * (Fr(2, 3) + k) / ((k + 1) ** 2)) % P)
        self.wc = LS(0, w, N + 3)
        # t_c: dt/dxi = wc/(1-xi)
        onemxi_inv = LS(0, [1] * (N + 3), N + 3)
        d = self.wc * onemxi_inv
        tc = [0] + [d.coef(k) * inv(k + 1) % P for k in range(0, N + 2)]
        self.tc = LS(0, tc, N + 3)
        # XX = theta log wc, theta = (xi-1) d/dxi
        self.XX = (self.xi + (P - 1)) * self.wc.deriv() * self.wc.inverse()

    def theta(self, f):
        return (self.xi + (P - 1)) * f.deriv()


class Rows:
    """Rows of draft Def 1.3 in normalised form; eps as Fractions."""

    def __init__(self, e1, e2, N, zmax, base=None, FF=None):
        self.e1, self.e2 = Fr(e1), Fr(e2)
        self.E = self.e1 + self.e2
        self.b = base or Base(N)
        self.N = N
        self.zmax = zmax
        self.E_ = fp(self.E)
        self.pp = fp(self.e1 * self.e2)
        self.Einv = inv(self.E_)
        self.FF = FF if FF is not None else self.b.XX
        # r_l
        r = [Fr(1)]
        rm = Fr(0)
        for m in range(0, zmax + 3):
            nxt = (3 * (3 * m * m + 3 * m + 1) * r[m] - m ** 3 * (r[m - 1] if m >= 1 else 0)) / (27 * (m + 1))
            r.append(nxt)
        self.r = [fp(x) for x in r]
        # S_n as dict pow->coef
        S = [dict() for _ in range(zmax + 2)]
        for l in range(0, zmax + 2):
            # series of 1/prod_{j<=l+1}(1+jz/E) in z
            ser = [0] * (zmax + 2)
            ser[0] = 1
            for j in range(1, l + 2):
                a = j * self.Einv % P
                for n in range(1, zmax + 2):
                    ser[n] = (ser[n] - a * ser[n - 1]) % P
            for n in range(l, zmax + 2):
                c = self.r[l] * ser[n - l] % P
                if c:
                    S[n][l] = (S[n].get(l, 0) + c) % P
        self.S = S
        self.Q = self.T(S)
        self.TQ = self.T(self.Q)
        self.cache_v = {}
        self.colour_data = {}
        for i in range(3):
            self.colour_data[i] = self._colour(i)
        self.ecache = {}

    def T(self, f):
        """(T f)_n = Y f_n + (Y/E) f_{n-1} + theta f_{n-1}, Y = -E/3 + 1/v, theta v = 1/E - v/3"""
        third = fp(Fr(1, 3))
        out = [dict() for _ in f]
        for n in range(len(f)):
            d = out[n]

            def add(pw, c):
                d[pw] = (d.get(pw, 0) + c) % P
            for k, c in f[n].items():
                add(k, -self.E_ * third % P * c)
                add(k - 1, c)
            if n >= 1:
                for k, c in f[n - 1].items():
                    # (Y/E) f_{n-1}
                    add(k, -third * c)
                    add(k - 1, self.Einv * c)
                    # theta v^k
                    add(k - 1, k * self.Einv % P * c)
                    add(k, -k * third % P * c)
        return out

    def _colour(self, i):
        b = self.b
        w = pow(OMEGA, i, P)
        ti = b.t * w
        vi = ti.inverse() * self.Einv + LS(0, [3 * self.Einv % P], b.N + 3)
        vinv = vi.inverse()
        pows = {}

        def vp(k):
            if k not in pows:
                if k == 0:
                    pows[k] = LS(0, [1], b.N + 3)
                elif k > 0:
                    pows[k] = vp(k - 1) * vi
                else:
                    pows[k] = vp(k + 1) * vinv
            return pows[k]

        def ev(lst):
            res = []
            for d in lst:
                s = None
                for k, c in d.items():
                    if c == 0:
                        continue
                    term = vp(k) * c
                    s = term if s is None else s + term
                if s is None:
                    s = LS(b.N + 3, [], b.N + 3)
                res.append(s)
            return res
        data = {}
        data['v'] = vi
        data['S'] = ev(self.S)
        data['Q'] = ev(self.Q)
        data['TQ'] = ev(self.TQ)
        data['Delta'] = vinv * vinv * (3 * self.pp % P * self.E_ % P)
        return data

    def Dil(self, i, d):
        # (-1)^d [z^{d-1}] S
        x = self.colour_data[i]['S'][d - 1]
        return x if d % 2 == 0 else -x

    def pair(self, i, j, M, b):
        """[z^M w^b] <P^i(z),P^j(w)>"""
        ci, cj = self.colour_data[i], self.colour_data[j]
        pp, E = self.pp, self.E_
        eta00 = pp * pow(E, 3, P) % P * inv(9) % P
        xi = self.b.xi
        r = ci['S'][M] * cj['S'][b] * eta00
        r = r + ci['Q'][M] * cj['Q'][b] * (-3 * pp * E % P)
        tmp = ci['Q'][M] * cj['TQ'][b] + ci['TQ'][M] * cj['Q'][b]
        r = r + tmp * xi * (-3 * pp % P)
        f = None
        if M >= 1:
            f = ci['Q'][M - 1] * cj['Q'][b]
        if b >= 1:
            g = ci['Q'][M] * cj['Q'][b - 1]
            f = g if f is None else f + g
        if f is not None:
            r = r + f * xi * self.FF * (3 * pp % P)
        return r

    def edge(self, i, j, a, b):
        key = (i, j, a, b)
        if key in self.ecache:
            return self.ecache[key]
        s = None
        for k in range(0, b + 1):
            M, bb = a + k + 1, b - k
            n = self.pair(i, j, M, bb)
            if (M + bb) % 2:
                n = -n
            term = n if k % 2 == 0 else -n
            s = term if s is None else s + term
        s = -s
        self.ecache[key] = s
        return s


class Vert:
    """Vertex integrals with Omega_w, w = (E,-e1,-e2) (default mode B)."""

    def __init__(self, e1, e2, w=None):
        e1, e2 = Fr(e1), Fr(e2)
        self.w = w or (e1 + e2, -e1, -e2)
        self.s = {}
        for c in range(1, 80, 2):
            pc = sum(Fr(1) / x ** c for x in self.w)
            self.s[c] = fp(-factorial(c - 1) * pc)
        self.cache = {}

    def VI(self, g, ds):
        ds = tuple(sorted(ds))
        key = (g, ds)
        if key in self.cache:
            return self.cache[key]
        n = len(ds)
        deg = 3 * g - 3 + n - sum(ds)
        if deg < 0:
            self.cache[key] = 0
            return 0
        tot = 0
        for chs in hodge.ch_partitions(deg):
            if g == 0 and chs:
                continue
            coef = 1
            for c in set(chs):
                m = chs.count(c)
                coef = coef * pow(self.s[c], m, P) % P * inv(factorial(m)) % P
            h = hodge.hodge(g, chs, ds)
            tot = (tot + coef * fp(h)) % P
        self.cache[key] = tot
        return tot


def partitions_ge2(budget):
    """multisets of d>=2 with sum(d-1) <= budget. yields tuple (sorted) and 1/prod mult!"""
    def rec(rem, mn):
        yield ()
        for d in range(mn, rem + 2):
            for tail in rec(rem - (d - 1), d):
                yield (d,) + tail
    seen = set()
    for x in rec(budget, 2):
        if x not in seen:
            seen.add(x)
            yield x


class Evaluator:
    def __init__(self, rows, vert):
        self.R = rows
        self.V = vert
        self.vc = {}

    def vertex(self, g, col, a_tuple):
        a_tuple = tuple(sorted(a_tuple))
        key = (g, col, a_tuple)
        if key in self.vc:
            return self.vc[key]
        n = len(a_tuple)
        budget = 3 * g - 3 + n - sum(a_tuple)
        s = None
        if budget >= 0:
            for D in partitions_ge2(budget):
                val = self.V.VI(g, a_tuple + D)
                if val == 0:
                    continue
                w = val
                for d in set(D):
                    w = w * inv(factorial(D.count(d))) % P
                term = LS(0, [w], 10 ** 6)
                for d in D:
                    term = term * self.R.Dil(col, d)
                s = term if s is None else s + term
        if s is None:
            s = LS(10 ** 6, [], 10 ** 6)
        Delta = self.R.colour_data[col]['Delta']
        s = s * Delta.pow(g - 1) if g != 1 else s
        self.vc[key] = s
        return s

    def graph(self, genera, edges, col):
        """decoration sum T(Gamma,c), no Aut factor. edges = list of (u,v). returns LS."""
        nv = len(genera)
        hes = [[] for _ in range(nv)]  # half-edge ids per vertex
        for k, (u, v) in enumerate(edges):
            hes[u].append((k, 0))
            hes[v].append((k, 1))
        maxa = [3 * genera[i] - 3 + len(hes[i]) for i in range(nv)]
        ne = len(edges)
        total = [None]

        def rec(k, cur, assign, used):
            if k == ne:
                # vertex factors
                res = cur
                for i in range(nv):
                    a_t = tuple(assign[h] for h in hes[i])
                    vv = self.vertex(genera[i], col[i], a_t)
                    if not vv.c and vv.aprec >= 10 ** 6:
                        return
                    res = res * vv
                total[0] = res if total[0] is None else total[0] + res
                return
            u, v = edges[k]
            for a in range(0, maxa[u] - used[u] + 1):
                for b in range(0, maxa[v] - used[v] + 1):
                    if u == v and a + b > maxa[u] - used[u]:
                        continue
                    e = self.R.edge(col[u], col[v], a, b)
                    assign[(k, 0)] = a
                    assign[(k, 1)] = b
                    used[u] += a
                    used[v] += b
                    rec(k + 1, cur * e, assign, used)
                    used[u] -= a
                    used[v] -= b
        rec(0, LS(0, [1], 10 ** 6), {}, [0] * nv)
        if total[0] is None:
            return LS(10 ** 6, [], 10 ** 6)
        return total[0]


def beta(n, e1, e2):
    """beta_n = (n-1)! [t^n] 1/(4 sinh sinh), exact"""
    e1, e2 = Fr(e1), Fr(e2)
    B = hodge.BERN

    def coef(k):  # [x^k] (x/2)/sinh(x/2)
        if k % 2:
            return Fr(0)
        return (Fr(2) ** (1 - k) - 1) * B[k] / factorial(k)
    s = Fr(0)
    for a in range(0, n + 3):
        b = n + 2 - a
        s += coef(a) * coef(b) * e1 ** a * e2 ** b
    return factorial(n - 1) * s / (e1 * e2)


def stable_graphs(g):
    """All connected stable graphs of genus g, no legs: list of (genera, edges, |Aut|)."""
    out = {}
    for nv in range(1, 2 * g - 1):
        # h1 = E - V + 1, E = g - sum g_v + V - 1
        for genera in product(range(g + 1), repeat=nv):
            if list(genera) != sorted(genera, reverse=True):
                continue
            h1 = g - sum(genera)
            if h1 < 0:
                continue
            ne = h1 + nv - 1
            pairs = [(i, j) for i in range(nv) for j in range(i, nv)]

            def rec(idx, rem, cur):
                if idx == len(pairs):
                    if rem == 0:
                        yield dict(cur)
                    return
                for m in range(rem + 1):
                    cur[pairs[idx]] = m
                    yield from rec(idx + 1, rem - m, cur)
                cur[pairs[idx]] = 0
            for mult in rec(0, ne, {}):
                # connected?
                adj = [set() for _ in range(nv)]
                val = [0] * nv
                for (i, j), m in mult.items():
                    if m:
                        adj[i].add(j)
                        adj[j].add(i)
                        val[i] += m
                        val[j] += m
                seen = {0}
                st = [0]
                while st:
                    x = st.pop()
                    for y in adj[x]:
                        if y not in seen:
                            seen.add(y)
                            st.append(y)
                if len(seen) != nv:
                    continue
                if any(2 * genera[i] - 2 + val[i] <= 0 for i in range(nv)):
                    continue
                # canonical form
                best = None
                naut = 0
                for perm in permutations(range(nv)):
                    gen2 = tuple(genera[perm[i]] for i in range(nv))
                    mm = tuple(mult[(min(perm[i], perm[j]), max(perm[i], perm[j]))] for (i, j) in pairs)
                    key = (gen2, mm)
                    if best is None or key > best:
                        best = key
                for perm in permutations(range(nv)):
                    if tuple(genera[perm[i]] for i in range(nv)) != tuple(genera):
                        continue
                    if all(mult[(min(perm[i], perm[j]), max(perm[i], perm[j]))] == mult[(i, j)] for (i, j) in pairs):
                        naut += 1
                autv = naut
                for (i, j), m in mult.items():
                    if i == j:
                        autv *= (2 ** m) * factorial(m)
                    else:
                        autv *= factorial(m)
                if best not in out:
                    edges = []
                    for (i, j), m in mult.items():
                        edges += [(i, j)] * m
                    out[best] = (tuple(genera), edges, autv)
    return list(out.values())
