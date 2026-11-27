"""Referee's TR recursion in the odd dxi-basis (see rf_tr.py docstring). Independent of the worker's code."""
from collections import defaultdict
from math import comb
from rf_tr import inv, zeta5, smul, sinv, spow, sderiv, ssub_neg, Local

class Curve:
    def __init__(self, P, KM):
        self.P, self.KM = P, KM
        z = zeta5(P)
        self.c = [pow(z, b, P) for b in range(5)]
        N = 4 * KM + 12
        self.N = N
        self.V0 = -(2 * KM + 2)
        self.E = 2 * KM + 6            # exponents V0..E stored
        self.loc = [Local(c, N, P) for c in self.c]
        self.nidx = 5 * (KM + 1)
        self._series()
        self._Ccache, self._Bcache, self._Acache = {}, {}, {}

    def idx(self, b, k): return b * (self.KM + 1) + k
    def bk(self, i): return divmod(i, self.KM + 1)

    def _series(self):
        P, N, KM, V0, E = self.P, self.N, self.KM, self.V0, self.E
        L = E - V0 + 1
        self.T = {}   # T[b',k][m] = [zeta'^(2k+1)] t_{b'}^m
        for bp in range(5):
            t = self.loc[bp].t
            pw = [1] + [0] * (2 * KM + 3)
            pows = [pw]
            for m in range(1, 2 * KM + 2):
                pw = smul(pw, t, 2 * KM + 4, P); pows.append(pw)
            for k in range(KM + 1):
                self.T[bp, k] = {m: pows[m][2 * k + 1] for m in range(1, 2 * k + 2)}
        self.s, self.ts, self.khat, self.tser = {}, {}, {}, {}
        for b in range(5):
            lo = self.loc[b]; t = lo.t + [0] * N; t = t[:N]; td = sderiv(t, P)
            self.tser[b] = t
            w = t[1:] + [0]                       # t = zeta w
            winv = sinv(w, N, P)
            D = [(2 * t[i + 1]) % P if i % 2 == 0 else 0 for i in range(N - 1)] + [0]   # (t(z)-t(-z))/z
            self.khat[b] = [(x * inv(2 * lo.alpha, P)) % P for x in sinv(D, N, P)]
            # e_{b',m} at c_b as arrays over exponents V0..E
            e = {}
            for bp in range(5):
                if bp == b:
                    wp = winv[:]                          # after the first product: w^-(m+1)
                    for m in range(1, 2 * KM + 2):
                        wp = smul(wp, winv, N, P)
                        ser = smul(td, wp, N, P)          # times zeta^-(m+1)
                        arr = [0] * L
                        for n_ in range(N):
                            ex = n_ - (m + 1)
                            if V0 <= ex <= E: arr[ex - V0] = ser[n_]
                        e[bp, m] = arr
                else:
                    d = (self.c[b] - self.c[bp]) % P
                    dt = t[:]; dt[0] = (dt[0] + d) % P
                    dinv = sinv(dt, N, P)
                    wp = dinv[:]
                    for m in range(1, 2 * KM + 2):
                        wp = smul(wp, dinv, N, P)
                        ser = smul(td, wp, N, P)
                        arr = [0] * L
                        for n_ in range(N):
                            if n_ <= E: arr[n_ - V0] = ser[n_]
                        e[bp, m] = arr
            for bp in range(5):
                for k in range(KM + 1):
                    arr = [0] * L
                    for m, cf in self.T[bp, k].items():
                        em = e[bp, m]
                        for q in range(L): arr[q] = (arr[q] + cf * em[q]) % P
                    i = self.idx(bp, k)
                    if bp == b:   # principal part must be exactly zeta^-(2k+2)
                        for ex in range(V0, 0):
                            want = 1 if ex == -(2 * k + 2) else 0
                            assert arr[ex - V0] == want, ("dxi principal part", b, k, ex)
                    self.s[b, i] = arr
                    self.ts[b, i] = [(-arr[q]) % P if (q + V0) % 2 == 0 else arr[q] for q in range(L)]  # -f(-z)

    def _A(self, b, k0, i):
        key = (b, k0, i)
        if key not in self._Acache:
            P, V0, E = self.P, self.V0, self.E
            kh, s = self.khat[b], self.s[b, i]
            # A = zeta^(2k0-1) khat * s, exponents V0+2k0-1 .. E ; store as dict ex->coef only for ex <= -1-V0
            A = {}
            top = -1 - V0
            for q, sv in enumerate(s):
                if sv:
                    ex0 = q + V0 + 2 * k0 - 1
                    for a in range(0, top - ex0 + 1):
                        if a < len(kh) and kh[a]:
                            A[ex0 + a] = (A.get(ex0 + a, 0) + sv * kh[a]) % P
            self._Acache[key] = A
        return self._Acache[key]

    def Clist(self, i, j):
        """[(a0, Res kappa_k0 s_i t_j)] over all branch points b and k0 <= KM (nonzero only)."""
        key = (i, j)
        if key in self._Ccache: return self._Ccache[key]
        P, V0, E = self.P, self.V0, self.E
        out = []
        for b in range(5):
            tj = self.ts[b, j]
            for k0 in range(self.KM + 1):
                A = self._A(b, k0, i)
                r = 0
                for ex, av in A.items():
                    ex2 = -1 - ex
                    if ex2 < V0: continue
                    assert ex2 <= E, "truncation too short"
                    r += av * tj[ex2 - V0]
                r %= P
                if r: out.append((self.idx(b, k0), r))
        self._Ccache[key] = out
        return out

    def _resk(self, b, k0, f, shift):
        """Res kappa_k0 zeta^shift f  (f an array over exponents V0..E)."""
        kh, V0, E = self.khat[b], self.V0, self.E
        r = 0
        for a, kv in enumerate(kh):
            ex = -1 - (2 * k0 - 1) - shift - a
            if ex < V0: break
            assert ex <= E, "truncation too short"
            r += kv * f[ex - V0]
        return r % self.P

    def Blist(self, j):
        """[(a0, u, Res kappa_k0 [sB_u t_j + s_j tB_u])], u = (b,ku): B(q,p_l) projected on its odd part."""
        if j in self._Bcache: return self._Bcache[j]
        P, out = self.P, []
        for b in range(5):
            diff = [(x - y) % P for x, y in zip(self.ts[b, j], self.s[b, j])]
            for k0 in range(self.KM + 1):
                for ku in range(self.KM + 1):
                    r = (2 * ku + 1) * self._resk(b, k0, diff, 2 * ku) % P
                    if r: out.append((self.idx(b, k0), self.idx(b, ku), r))
        self._Bcache[j] = out
        return out
