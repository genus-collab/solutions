# -*- coding: utf-8 -*-
"""xy_trgen.py -- direct TR mod p for a GENERIC point of the Hurwitz family (no Z/5 symmetry used):
   x = z^4 + c3 z^3 + c2 z^2 + c1 z + V/z,  y = z,  B = dz dz'/(z - z')^2,
with the five branch points r_j chosen in F_p first (z^2 x' = 4 prod (z - r_j) forces Sum 1/r_j = 0).
Same kernel convention as xy_tr.py.  Checks F_g = B_2g/(2g(2g-2)) V^(2-2g) off the symmetric line.
Usage: py -I -B xy_trgen.py [p89|p61] [seed] [gmax]
"""
import os, sys, random, itertools, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from math import comb
from xy_closed import bernoulli

P = {"p89": 2 ** 89 - 1, "p61": 2 ** 61 - 1}[sys.argv[1] if len(sys.argv) > 1 else "p89"]
inv = lambda a: pow(a % P, P - 2, P)


def smul(a, b, n):
    out = [0] * (n + 1)
    for i, x in enumerate(a[:n + 1]):
        if x:
            for j, y in enumerate(b[:n + 1 - i]): out[i + j] += x * y
    return [c % P for c in out]


def sinv(a, n):
    a0 = inv(a[0]); out = [a0] + [0] * n
    for k in range(1, n + 1):
        out[k] = -sum(a[i] * out[k - i] for i in range(1, min(k, len(a) - 1) + 1)) * a0 % P
    return out


def ssqrt1(a, n):
    b = [1] + [0] * n; h = inv(2)
    for k in range(1, n + 1): b[k] = (a[k] - sum(b[i] * b[k - i] for i in range(1, k))) * h % P
    return b


def scomp(f, g, n):
    out = [0] * (n + 1)
    for c in reversed(f[:n + 1]):
        out = smul(out, g, n); out[0] = (out[0] + c) % P
    return out


def srev(f, n):
    g = [0, inv(f[1])] + [0] * (n - 1); i1 = inv(f[1])
    for k in range(2, n + 1): g[k] = -scomp(f, g, k)[k] * i1 % P
    return g


def lmul(A, B, top):
    va, ca = A; vb, cb = B; v = va + vb; n = top - v + 1
    if n <= 0: return (v, [])
    out = [0] * n
    for i, a in enumerate(ca[:n]):
        if a:
            for j, b in enumerate(cb[:n - i]): out[i + j] += a * b
    return (v, [c % P for c in out])


def ladd(acc, A, cf=1):
    v, c = A
    for j, a in enumerate(c): acc[v + j] = (acc.get(v + j, 0) + a * cf) % P


class Curve:
    def __init__(self, seed):
        rng = random.Random(seed)
        r = [rng.randrange(1, P) for _ in range(4)]
        r.append(inv(-sum(inv(x) for x in r)))
        e = [1, 0, 0, 0, 0, 0]
        for x in r:
            for k in range(5, 0, -1): e[k] = (e[k] + e[k - 1] * x) % P
        assert e[4] == 0
        self.r = r; self.c3 = -4 * e[1] * inv(3) % P; self.c2 = 2 * e[2] % P; self.c1 = -4 * e[3] % P; self.V = 4 * e[5] % P

    def taylor_x(self, r, n):
        out = []
        for k in range(n + 1):
            s = comb(4, k) * pow(r, 4 - k, P) if k <= 4 else 0
            if k <= 3: s += self.c3 * comb(3, k) * pow(r, 3 - k, P)
            if k <= 2: s += self.c2 * comb(2, k) * pow(r, 2 - k, P)
            if k <= 1: s += self.c1 * pow(r, 1 - k, P)
            s += self.V * (-1) ** k * pow(inv(r), k + 1, P)
            out.append(s % P)
        return out


class TRgen:
    def __init__(self, cur, gmax, nmax):
        self.C = cur; kmax = max(6 * g - 5 + 2 * n for g in range(gmax + 1) for n in range(1, nmax + 1) if 2 * g - 2 + n > 0)
        self.kmax = kmax; self.Pt = kmax + 3; N = 2 * kmax + 12; self.N = N; self.om = {}; self.loc = []
        for j, r in enumerate(cur.r):
            c = cur.taylor_x(r, N + 3); assert c[1] == 0
            u = [1] + [c[k + 2] * inv(c[2]) % P for k in range(1, N + 1)]
            tau = [0] + ssqrt1(u, N)[:N]; sig = scomp(srev(tau, N), [-a % P for a in tau], N)
            dsig = [k * sig[k] % P for k in range(1, N + 1)] + [0]; su = sig[1:] + [0]
            dx = [(k + 1) * c[k + 1] % P for k in range(1, N + 2)]          # x'(r+t)/t
            i2 = sinv([2 * a % P for a in dx], N)
            kap = {}
            for m in range(1, kmax + 4):
                q = [0] * (N + 1); sj = [1] + [0] * N
                for jj in range(m):
                    for i in range(N + 1 - (m - 1 - jj)): q[i + m - 1 - jj] = (q[i + m - 1 - jj] + sj[i]) % P
                    sj = smul(sj, sig, N)
                kap[m] = (-1, smul(q, i2, N))
            w = sinv([(1 if k == 1 else 0) - sig[k] for k in range(1, N + 1)] + [0], N)
            Bss = (-2, smul(smul(w, w, N), dsig, N))
            sp = [dsig[:N + 1]]
            for n in range(1, N + 1): sp.append(smul(sp[-1], sig, N))
            sun = {1: sinv(su, N)}
            for k in range(2, kmax + 3): sun[k] = smul(sun[k - 1], sun[1], N)
            self.loc.append(dict(kap=kap, Bss=Bss, sp=sp, sun=sun, dsig=dsig))
        self._ep = {}; self._es = {}

    def e_p(self, j, i, k):
        key = (j, i, k)
        if key not in self._ep:
            if i == j: s = (-(k + 1), [1])
            else:
                w = inv(self.C.r[j] - self.C.r[i])
                s = (0, [(-1) ** n * comb(k + n, n) * pow(w, k + 1 + n, P) % P for n in range(self.Pt + 1)])
            self._ep[key] = s
        return self._ep[key]

    def e_s(self, j, i, k):
        key = (j, i, k)
        if key not in self._es:
            L = self.loc[j]
            if i == j:
                s = (-(k + 1), smul(L["dsig"], L["sun"][k + 1], self.N)[: self.Pt + k + 2])
            else:
                w = inv(self.C.r[j] - self.C.r[i]); acc = [0] * (self.Pt + 1)
                for n in range(self.Pt + 1):
                    cf = (-1) ** n * comb(k + n, n) * pow(w, k + 1 + n, P) % P
                    for t in range(n, self.Pt + 1): acc[t] = (acc[t] + cf * L["sp"][n][t]) % P
                s = (0, acc)
            self._es[key] = s
        return self._es[key]

    def expand(self, j, g, n, sg):
        L = self.loc[j]; out = {}
        if (g, n) == (0, 2):
            for m in range(self.Pt + 1):
                out[((j, m + 1),)] = (0, [(m + 1) * a % P for a in L["sp"][m][: self.Pt + 1]]) if sg else (m, [m + 1])
            return out
        acc = {}
        for key, cf in self.om[(g, n)].items():
            ser = self.e_s(j, *key[0]) if sg else self.e_p(j, *key[0])
            ladd(acc.setdefault(key[1:], {}), ser, cf)
        for rr, d in acc.items():
            v = min(d); out[rr] = (v, [d.get(e, 0) for e in range(v, self.Pt + 1)])
        return out

    def step(self, g, n1):
        n = n1 - 1; S = tuple(range(n)); out = {}
        for j in range(5):
            L = self.loc[j]; R = {}
            if g >= 1:
                if (g - 1, n + 2) == (0, 2): ladd(R.setdefault((), {}), (L["Bss"][0], L["Bss"][1][:3]))
                else:
                    for key, cf in self.om[(g - 1, n + 2)].items():
                        ladd(R.setdefault(key[2:], {}), lmul(self.e_p(j, *key[0]), self.e_s(j, *key[1]), 0), cf)
            for g1 in range(g + 1):
                for sz in range(n + 1):
                    for I in itertools.combinations(S, sz):
                        J = tuple(s for s in S if s not in I)
                        if (g1, len(I)) == (0, 0) or (g - g1, len(J)) == (0, 0): continue
                        E1 = self.expand(j, g1, len(I) + 1, False); E2 = self.expand(j, g - g1, len(J) + 1, True)
                        for r1, s1 in E1.items():
                            for r2, s2 in E2.items():
                                rest = [None] * n
                                for pos, var in enumerate(I): rest[var] = r1[pos]
                                for pos, var in enumerate(J): rest[var] = r2[pos]
                                ladd(R.setdefault(tuple(rest), {}), lmul(s1, s2, 0))
            for rest, d in R.items():
                for m in range(1, 6 * g - 3 + 2 * n1):
                    kv, kc = L["kap"][m]; v = sum(a * kc[-1 - e - kv] for e, a in d.items() if 0 <= -1 - e - kv < len(kc)) % P
                    if v: out[((j, m),) + rest] = v
        assert all(k[0][1] <= 6 * g - 5 + 2 * n1 for k in out)
        self.om[(g, n1)] = out

    def free_energy(self, g):
        C = self.C; tot = 0
        for ((j, k),), cf in self.om[(g, 1)].items():
            r = C.r[j]
            poly = [4 * inv(5), 3 * C.c3 * inv(4), 2 * C.c2 * inv(3), C.c1 * inv(2)]   # z^5, z^4, z^3, z^2
            ph = sum(cf_ * comb(e, k) * pow(r, e - k, P) for cf_, e in zip(poly, (5, 4, 3, 2)) if k <= e)
            ph -= C.V * (-1) ** (k + 1) * inv(k) * pow(inv(r), k, P)
            tot += cf * ph
        return tot * inv(2 - 2 * g) % P


def main():
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    gmax = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    cur = Curve(seed); t0 = time.time()
    print("p =", P, "seed", seed, "branch points", cur.r[:2], "...; c3,c2,c1,V nonzero:", all([cur.c3, cur.c2, cur.c1, cur.V]))
    tr = TRgen(cur, gmax, 2 if gmax == 2 else 3)
    order = [(0, 3), (1, 1), (0, 4), (1, 2), (2, 1)] + ([(0, 5), (1, 3), (2, 2), (3, 1)] if gmax >= 3 else [])
    for (g, n) in order: tr.step(g, n)
    B = bernoulli(2 * gmax + 2)
    for g in range(2, gmax + 1):
        Fg = tr.free_energy(g); c = B[2 * g] / (2 * g * (2 * g - 2))
        exp = c.numerator * inv(c.denominator) * pow(inv(cur.V), 2 * g - 2, P) % P
        print("CHECK F_%d = B_2g/(2g(2g-2)) V^(2-2g) at a generic point of the Hurwitz family: %s  (%.1fs)" % (g, Fg == exp, time.time() - t0))


if __name__ == "__main__":
    main()
