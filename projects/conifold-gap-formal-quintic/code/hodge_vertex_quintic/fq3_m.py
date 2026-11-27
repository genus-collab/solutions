# -*- coding: utf-8 -*-
"""fq3_m.py -- symplectic matrix series to zhat-order 6 for the genus-3 programme: the Lie algebra sp (75 coordinates), the classes
of N, the tautological M_0 = R~ mod zhat^7 from the conifold symbol, the u-dependent factorisation M(zhat; u) = M_0 exp(Sum_k u^k N_k)
(RIGHT placement; LEFT: exp(Sum u^k N_k) M_0), and the order-6 row factorisation R' = R M^{-1}.

sp mod zhat^7: N = Sum_{k=1}^{6} N_k zhat^k with N_k^T = (-1)^{k+1} N_k (symmetric for odd k, antisymmetric for even k): coordinates
  v = 0..74 in blocks  k=1: 0..14 (PAIRS1, i <= j), k=2: 15..24 (PAIRS2, i < j), k=3: 25..39, k=4: 40..49, k=5: 50..64, k=6: 65..74
  (fq_mfit.PAIRS1/PAIRS2 conventions; the first 40 coordinates are those of fq_leadsolve.mu_params for k <= 3).
Classes (linear parametrisations N = A t, A = list over the 75 coordinates of {t_index: coeff}):
  cov   t = (b_1, b_3, b_5): N_k = b_k diag(zeta^{-ki}) for k = 1, 3, 5 (Z/5-covariant scalar; N_5 = b_5 . 1);
  diag  t = (X_{k,a}), k = 1, 3, 5, a = 0..4: (N_k)_{ii} = (1/5) Sum_a zeta^{-ai} X_{k,a} (Fourier modes of the diagonal class);
  circ  Z/5-covariant: (N_k)_{ij} = zeta^{-ki} c_k(j-i), c_k(5-d) = (-1)^{k+1} zeta^{kd} c_k(d) (d = 1, 2), c_k(0) = 0 for even k:
        t = (c_1(0), c_1(1), c_1(2), c_2(1), c_2(2), c_3(0), c_3(1), c_3(2), c_4(1), c_4(2), c_5(0), c_5(1), c_5(2), c_6(1), c_6(2)) (15);
  full  the 75 coordinates.
MConst (fq3_mvertex): constant M, lists .M, .Minv of 7 matrices.  M0_from_symbol(sym, R): M_0 = R~ mod zhat^7 with M_k = C_k
(fq_symbol.Symbol.colour_matrix) and M^{-1}_k = (-1)^k C_k^T (M(zhat)^T M(-zhat) = 1); returns (MConst, mu = log, checks).
MUprod(R, M0, Ls, Rs, K): two-sided exp(Sum u^k L_k) M_0 exp(Sum u^k R_k).  MU(R, M0, Ns, K, placement): M(zhat; u) mod (zhat^7, u^K), Ns = [N_1, N_2, ...] (sp series), E(u) = exp(Sum_k u^k N_k) computed as a
bivariate truncated series; .M[k][i][j] and .Minv[k][i][j] are u-polynomials (lists of K F_p coefficients), .at(eps) -> MConst at u = eps,
.M_ser / .Minv_ser (Ser in u over R), .is_identity = False, .kmax = 6, .K.
rows_prime6(S, M): Q'^{(i)}_{a,n} = Sum_{k <= min(n,6)} Sum_j Q^{(j)}_{a,n-k} (M^{-1}_k)_{ji}, M^{-1} entries F_p (MConst) or Ser (MU.Minv_ser).
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fq_series import Ser
from fq_mfit import PAIRS1, PAIRS2
from fq_symbol import mat_series_mul, mat_series_exp, mat_series_log, _mat_add, _mat_zero
from fq_sser import mat_T, mat_eq, mat_mul, mat_id
from fq3_mvertex import MConst, KMAX

BLOCKS = []            # (k, pairs, offset)
_off = 0
for _k in range(1, KMAX + 1):
    _pairs = PAIRS1 if _k % 2 == 1 else PAIRS2
    BLOCKS.append((_k, _pairs, _off)); _off += len(_pairs)
NSP = _off             # 75
SP_NAMES = ["n%d_%d%d" % (k, i, j) for (k, pairs, off) in BLOCKS for (i, j) in pairs]


def Z5(): return [[0] * 5 for _ in range(5)]


def sp_mats(par, p):
    """75 coordinates -> [0, N_1, .., N_6]"""
    N = [Z5() for _ in range(KMAX + 1)]
    for k, pairs, off in BLOCKS:
        for v, (i, j) in enumerate(pairs):
            x = par[off + v] % p
            N[k][i][j] = x; N[k][j][i] = x if k % 2 == 1 else (-x) % p
    return N


def sp_params(N, p):
    assert in_sp(N, p), "not in sp"
    return [N[k][i][j] % p for k, pairs, off in BLOCKS for (i, j) in pairs]


def in_sp(N, p):
    if any(x % p for r in N[0] for x in r): return False
    return all((N[k][i][j] - (-1) ** (k + 1) * N[k][j][i]) % p == 0 for k in range(1, KMAX + 1) for i in range(5) for j in range(5))


def sp_direction(v, p):
    par = [0] * NSP; par[v] = 1; return sp_mats(par, p)


def sp_add(X, Y, p, c=1): return [_mat_add(X[k], Y[k], p, c) for k in range(KMAX + 1)]
def sp_scale(X, c, p): return [[[x * c % p for x in r] for r in Xk] for Xk in X]


def bracket(A, B, p):
    return [_mat_add(mat_series_mul(A, B, p, KMAX)[k], mat_series_mul(B, A, p, KMAX)[k], p, p - 1) for k in range(KMAX + 1)]


def conj(M0, N, p):
    """M_0^{-1} N M_0 mod zhat^7 (LEFT-placed N' <-> RIGHT-placed N = M_0^{-1} N' M_0)"""
    return mat_series_mul(mat_series_mul(M0.Minv, N, p, KMAX), M0.M, p, KMAX)


# ------------------------------------------------------------------ classes
def class_cov(R):
    A = [dict() for _ in range(NSP)]; t = 0; names = []
    for k, pairs, off in BLOCKS:
        if k % 2 == 0: continue
        for i in range(5): A[off + pairs.index((i, i))] = {t: R.z(-k * i)}
        names.append("b%d" % k); t += 1
    return A, names


def class_diag(R):
    p = R.p; inv5 = R.inv(5); A = [dict() for _ in range(NSP)]; t = 0; names = []
    for k, pairs, off in BLOCKS:
        if k % 2 == 0: continue
        for i in range(5): A[off + pairs.index((i, i))] = {t + a: inv5 * R.z(-a * i) % p for a in range(5)}
        names += ["X%d_%d" % (k, a) for a in range(5)]; t += 5
    return A, names


def class_circ(R):
    p = R.p; A = [dict() for _ in range(NSP)]; names = []; t = 0
    for k, pairs, off in BLOCKS:
        sgn = 1 if k % 2 == 1 else p - 1
        # free parameters of this k: c_k(0) (odd k only), c_k(1), c_k(2)
        red = {1: (1, 1), 2: (2, 1), 3: (2, sgn * R.z(2 * k) % p), 4: (1, sgn * R.z(k) % p)}
        if k % 2 == 1: red[0] = (0, 1); loc = {0: t, 1: t + 1, 2: t + 2}; names += ["c%d_0" % k, "c%d_1" % k, "c%d_2" % k]; t += 3
        else: loc = {1: t, 2: t + 1}; names += ["c%d_1" % k, "c%d_2" % k]; t += 2
        for v, (i, j) in enumerate(pairs):
            d = (j - i) % 5
            if d not in red: continue
            dd, c = red[d]; A[off + v] = {loc[dd]: c * R.z(-k * i) % p}
    return A, names


def class_full(R):
    return [{v: 1} for v in range(NSP)], list(SP_NAMES)


CLASSES = {'cov': class_cov, 'diag': class_diag, 'circ': class_circ, 'full': class_full}


def class_columns(R, cls, M0=None, placement='right'):
    """list over the class parameters of 75-vectors (right-placed coordinates); left: conjugated by M_0"""
    p = R.p; A, names = CLASSES[cls](R); cols = []
    for j in range(len(names)):
        col = [A[v].get(j, 0) % p for v in range(NSP)]
        if placement == 'left': col = sp_params(conj(M0, sp_mats(col, p), p), p)
        cols.append(col)
    return cols, names


def expand(cols, t, p):
    return [sum(c[v] * t[j] for j, c in enumerate(cols)) % p for v in range(NSP)]


# ------------------------------------------------------------------ M_0 from the symbol
def M0_from_symbol(sym, R, kmax=KMAX):
    """M_0 = R~ mod zhat^{kmax+1}: (MConst, mu = log M_0, dict of checks)"""
    p = R.p; C = sym.colour_matrix()
    assert len(C) >= kmax + 1, "symbol has only %d zhat-orders" % (len(C) - 1)
    M = [C[k] for k in range(kmax + 1)]
    Minv = [[[(C[k][j][i] * (1 if k % 2 == 0 else p - 1)) % p for j in range(5)] for i in range(5)] for k in range(kmax + 1)]
    M0 = MConst(R, M, Minv, check=True)                       # asserts symplecticity and M M^-1 = 1 to order kmax
    mu = mat_series_log(C, R, kmax)
    checks = dict(exp_log=all(mat_eq(mat_series_exp(mu, R, kmax)[k], C[k], p) for k in range(kmax + 1)),
                  mu_in_sp=in_sp(mu, p), symbol_sympl=not sym.symplectic_check(kmax))
    return M0, mu, checks


# ------------------------------------------------------------------ u-dependent M
def _biv_mul(X, Y, p, K):
    """product of bivariate truncated series {(k, j): matrix}, zhat-order k <= KMAX, u-order j < K"""
    out = {}
    for (k1, j1), A in X.items():
        for (k2, j2), B in Y.items():
            k = k1 + k2; j = j1 + j2
            if k > KMAX or j >= K: continue
            P = mat_mul(A, B, p)
            key = (k, j)
            out[key] = _mat_add(out[key], P, p) if key in out else P
    return out


def _biv_exp(X, p, K, R):
    """exp of a bivariate series with X[(0, 0)] absent and every term of zhat-order >= 1 and u-order >= 1"""
    assert all(k >= 1 and j >= 1 for (k, j) in X)
    out = {(0, 0): mat_id(5, p)}; power = {(0, 0): mat_id(5, p)}; fact = 1
    for m in range(1, min(KMAX, K - 1) + 1):
        power = _biv_mul(power, X, p, K); fact *= m; c = R.inv(fact)
        for key, A in power.items():
            out[key] = _mat_add(out.get(key, Z5()), A, p, c)
    return out


class MU(object):
    """M(zhat; u) = M_0 exp(Sum_k u^k N_k) (right) or exp(Sum_k u^k N_k) M_0 (left), mod (zhat^7, u^K)"""
    def __init__(self, R, M0, Ns, K, placement='right'):
        p = R.p; self.R = R; self.p = p; self.K = K; self.kmax = KMAX; self.placement = placement; self.M0 = M0; self.Ns = Ns
        self.is_identity = False
        X = {}
        for jj, N in enumerate(Ns):
            j = jj + 1
            if j >= K: break
            for k in range(1, KMAX + 1):
                if any(x % p for r in N[k] for x in r): X[(k, j)] = [[x % p for x in r] for r in N[k]]
        E = _biv_exp(X, p, K, R); Ei = _biv_exp({key: [[(-x) % p for x in r] for r in A] for key, A in X.items()}, p, K, R)
        B0 = {(k, 0): M0.M[k] for k in range(KMAX + 1)}; B0i = {(k, 0): M0.Minv[k] for k in range(KMAX + 1)}
        if placement == 'right': Mb = _biv_mul(B0, E, p, K); Mbi = _biv_mul(Ei, B0i, p, K)
        else: Mb = _biv_mul(E, B0, p, K); Mbi = _biv_mul(B0i, Ei, p, K)
        # u-polynomial entries: self.M[k][i][j] = [c_0, .., c_{K-1}]
        self.M = [[[[Mb.get((k, j), Z5())[a][b] for j in range(K)] for b in range(5)] for a in range(5)] for k in range(KMAX + 1)]
        self.Minv = [[[[Mbi.get((k, j), Z5())[a][b] for j in range(K)] for b in range(5)] for a in range(5)] for k in range(KMAX + 1)]
        self.M_ser = [[[Ser(R, 0, list(self.M[k][a][b])) for b in range(5)] for a in range(5)] for k in range(KMAX + 1)]
        self.Minv_ser = [[[Ser(R, 0, list(self.Minv[k][a][b])) for b in range(5)] for a in range(5)] for k in range(KMAX + 1)]

    def at(self, eps, check=False):
        """MConst at u = eps (F_p)"""
        p = self.p; e = eps % p
        def ev(c):
            r = 0
            for x in reversed(c): r = (r * e + x) % p
            return r
        M = [[[ev(self.M[k][a][b]) for b in range(5)] for a in range(5)] for k in range(KMAX + 1)]
        Mi = [[[ev(self.Minv[k][a][b]) for b in range(5)] for a in range(5)] for k in range(KMAX + 1)]
        return MConst(self.R, M, Mi, check)

    def check(self):
        """M M^{-1} = 1 and M(zhat)^T M(-zhat) = 1 mod (zhat^7, u^K), and at(0) = M_0"""
        p = self.p; K = self.K; ok = True
        Mb = {(k, j): [[self.M[k][a][b][j] for b in range(5)] for a in range(5)] for k in range(KMAX + 1) for j in range(K)}
        Mbi = {(k, j): [[self.Minv[k][a][b][j] for b in range(5)] for a in range(5)] for k in range(KMAX + 1) for j in range(K)}
        P = _biv_mul(Mb, Mbi, p, K)
        MT = {key: mat_T(A) for key, A in Mb.items()}
        Mm = {(k, j): [[(x if k % 2 == 0 else (-x) % p) for x in r] for r in A] for (k, j), A in Mb.items()}
        S = _biv_mul(MT, Mm, p, K)
        for key in set(P) | set(S):
            tgt = mat_id(5, p) if key == (0, 0) else Z5()
            if not mat_eq(P.get(key, Z5()), tgt, p) or not mat_eq(S.get(key, Z5()), tgt, p): ok = False
        M00 = self.at(0)
        ok = ok and all(mat_eq(M00.M[k], self.M0.M[k], p) and mat_eq(M00.Minv[k], self.M0.Minv[k], p) for k in range(KMAX + 1))
        return ok


class MUprod(MU):
    """two-sided dressing M(zhat; u) = exp(Sum_k u^k L_k) M_0 exp(Sum_k u^k R_k) mod (zhat^7, u^K); Ls, Rs lists of sp series"""
    def __init__(self, R, M0, Ls, Rs, K):
        p = R.p; self.R = R; self.p = p; self.K = K; self.kmax = KMAX; self.placement = 'both'; self.M0 = M0; self.Ns = (Ls, Rs)
        self.is_identity = False
        def biv(Ns):
            X = {}
            for jj, N in enumerate(Ns):
                j = jj + 1
                if j >= K: break
                for k in range(1, KMAX + 1):
                    if any(x % p for r in N[k] for x in r): X[(k, j)] = [[x % p for x in r] for r in N[k]]
            return X
        XL = biv(Ls); XR = biv(Rs)
        EL = _biv_exp(XL, p, K, R); ELi = _biv_exp({key: [[(-x) % p for x in r] for r in A] for key, A in XL.items()}, p, K, R)
        ER = _biv_exp(XR, p, K, R); ERi = _biv_exp({key: [[(-x) % p for x in r] for r in A] for key, A in XR.items()}, p, K, R)
        B0 = {(k, 0): M0.M[k] for k in range(KMAX + 1)}; B0i = {(k, 0): M0.Minv[k] for k in range(KMAX + 1)}
        Mb = _biv_mul(_biv_mul(EL, B0, p, K), ER, p, K); Mbi = _biv_mul(_biv_mul(ERi, B0i, p, K), ELi, p, K)
        self.M = [[[[Mb.get((k, j), Z5())[a][b] for j in range(K)] for b in range(5)] for a in range(5)] for k in range(KMAX + 1)]
        self.Minv = [[[[Mbi.get((k, j), Z5())[a][b] for j in range(K)] for b in range(5)] for a in range(5)] for k in range(KMAX + 1)]
        self.M_ser = [[[Ser(R, 0, list(self.M[k][a][b])) for b in range(5)] for a in range(5)] for k in range(KMAX + 1)]
        self.Minv_ser = [[[Ser(R, 0, list(self.Minv[k][a][b])) for b in range(5)] for a in range(5)] for k in range(KMAX + 1)]


# ------------------------------------------------------------------ rows
def rows_prime6(S, M):
    """Q'^{(i)}_{a,n} = Sum_k Sum_j Q^{(j)}_{a,n-k} (M^{-1}_k)_{ji}; M^{-1} entries F_p ints (MConst) or Ser (MU.Minv_ser)"""
    R = S.R; p = R.p
    if getattr(M, 'is_identity', False): return S.Q
    Minv = M.Minv_ser if hasattr(M, 'Minv_ser') else M.Minv
    kmax = len(Minv) - 1
    Qp = [[[None] * (S.Nz + 1) for i in range(5)] for a in range(5)]
    for a in range(5):
        for i in range(5):
            for n in range(S.Nz + 1):
                acc = None
                for k in range(0, min(n, kmax) + 1):
                    for j in range(5):
                        c = Minv[k][j][i]
                        if isinstance(c, Ser):
                            if not c.c: continue
                            t = S.Q[a][j][n - k].mul(c)
                        else:
                            if not c % p: continue
                            t = S.Q[a][j][n - k].scale(c)
                        acc = t if acc is None else acc.add(t)
                Qp[a][i][n] = acc if acc is not None else Ser(R, S.Q[a][i][n].prec, [])
    return Qp


if __name__ == "__main__":
    import random
    from fq_fp import FpRing, P2
    R = FpRing(P2); p = R.p; rng = random.Random(11)
    # sp coordinates round trip, classes in sp
    par = [rng.randrange(p) for _ in range(NSP)]; N = sp_mats(par, p); assert in_sp(N, p) and sp_params(N, p) == par
    for cls in CLASSES:
        A, names = CLASSES[cls](R)
        for j in range(len(names)):
            col = [A[v].get(j, 0) % p for v in range(NSP)]; assert in_sp(sp_mats(col, p), p), (cls, names[j])
        print("class %s: %d parameters %s" % (cls, len(names), names if len(names) <= 15 else "..."))
    # circulant class: covariance (N_k)_{i+1,j+1} = zeta^{-k} (N_k)_{ij}
    A, names = class_circ(R); t = [rng.randrange(p) for _ in names]
    Nc = sp_mats([sum(A[v].get(j, 0) * t[j] for j in range(len(t))) % p for v in range(NSP)], p)
    assert all((Nc[k][(i + 1) % 5][(j + 1) % 5] - R.z(-k) * Nc[k][i][j]) % p == 0 for k in range(1, 7) for i in range(5) for j in range(5))
    # covariant class is inside the circulant class and the diagonal class
    Ac, nc = class_cov(R); b = [rng.randrange(p) for _ in nc]
    Ncov = sp_mats([sum(Ac[v].get(j, 0) * b[j] for j in range(3)) % p for v in range(NSP)], p)
    assert all((Ncov[k][(i + 1) % 5][(j + 1) % 5] - R.z(-k) * Ncov[k][i][j]) % p == 0 for k in range(1, 7) for i in range(5) for j in range(5))
    assert all(Ncov[k][i][j] == 0 for k in range(1, 7) for i in range(5) for j in range(5) if i != j or k % 2 == 0)
    print("classes ok (circulant covariance, cov inside circ)")
    # MU: random N_1, N_2, N_3 and M_0 = exp(random sp), K = 6; checks; at(eps) equals direct exponential at eps
    mu0 = sp_mats([rng.randrange(p) for _ in range(NSP)], p); M0 = MConst.from_log(R, mu0)
    Ns = [sp_mats([rng.randrange(p) for _ in range(NSP)], p) for _ in range(3)]
    for placement in ('right', 'left'):
        MUo = MU(R, M0, Ns, 6, placement); assert MUo.check(), placement
        eps = rng.randrange(p)
        X = sp_add(sp_add(sp_scale(Ns[0], eps, p), sp_scale(Ns[1], eps * eps % p, p), p), sp_scale(Ns[2], pow(eps, 3, p), p), p)
        E = mat_series_exp(X, R, KMAX)
        direct = mat_series_mul(M0.M, E, p, KMAX) if placement == 'right' else mat_series_mul(E, M0.M, p, KMAX)
        # NOTE: MU truncates in u at K = 6: with eps generic the comparison needs the full u-dependence -> compare with K large enough
        MUf = MU(R, M0, Ns, 19, placement); Me = MUf.at(eps, check=True)
        assert all(mat_eq(Me.M[k], direct[k], p) for k in range(KMAX + 1)), placement
        print("MU %s: bivariate exponential, inverse and symplecticity ok; at(eps) = M_0 exp(eps N_1 + eps^2 N_2 + eps^3 N_3) exactly (K = 19)" % placement)
    # conj: right-placed equivalent of a left-placed N'
    Np = Ns[0]; Nr = conj(M0, Np, p)
    L = MU(R, M0, [Np], 4, 'left'); Rr = MU(R, M0, [Nr], 4, 'right')
    assert all(L.M[k][a][b] == Rr.M[k][a][b] for k in range(KMAX + 1) for a in range(5) for b in range(5))
    print("conj: exp(uN') M_0 = M_0 exp(u M_0^-1 N' M_0) ok")
    print("fq3_m ok")
