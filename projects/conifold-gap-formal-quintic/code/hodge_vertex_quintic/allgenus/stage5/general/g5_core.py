"""g5_core.py -- Stage 5 (general): Laurent expansion in eps of a normalised formal stationary-phase (SP) series
with ONE factor singular at the colour point.  New code, stdlib only; imports no project or agent module.

Setting (sigma = z s).  The 1-form  Gamma_St(e* + c*(sigma - s*)/z)^M  B(sigma; z)  q^(sigma/z) dsigma, with
  B(sigma; z) = exp(Phi(sigma)/z + A_0(sigma) + Sum_(n>=1) z^n A_n(sigma))   regular at sigma = s*  (Taylor data given),
  Gamma_St(e + Y) = exp((Y + e - 1/2) log Y - Y + (1/2) log 2pi + Sum_n beta_n Y^-n),  beta_n = beta_n(e) (Bernoulli).
Critical point sigma_c = s* + eps; q is then fixed by the critical equation, and no derivative of order >= 2 of the
phase involves q, so the normalised SP series depends on eps only.  Scaled variables sigma = s* + eps (1 + T),
nu = z/eps.  The normalised SP is a double series ct(nu, eps) = Sum_n nu^n ct_n(eps), and the z^n coefficient is
   c_n(sigma_c) = eps^-n ct_n(eps).
Removable singularity at eps = 0  <=>  ord_eps ct_n >= n for every n;  value c_n(s*) = [eps^n] ct_n.
Method (Lagrange-Buermann, exact):  phase = F(T)/nu with F = -(a/2) T^2 Theta(T); then
   normalised SP of e^(F/nu) g(T) dT = Sum_m (2m-1)!! (nu/a)^m [T^(2m)] g(T) Theta(T)^(-(2m+1)/2),
and g is expanded in powers of U = 1 + T:
   g = C(eps) U^(M(e*-1/2)) Sum_p S_p nu^p U^-p  Sum_(k,r) Bk[k][r] nu^k eps^(k+r) U^r,
   S_p = [y^p] exp(M Sum_n beta*_n c*^-n y^n),   Bk = [z^k x^r] exp(A_0(s*+x) - A_0(s*) + Sum_n z^n A_n(s*+x)),
   C = 1/Bk[0](eps).   Field: exact Q (P = None) or F_p.
"""
from fractions import Fraction as Fr
from math import comb

P = None


def setfield(p):
    global P
    P = p


def cv(x):
    x = Fr(x)
    if P is None:
        return x
    return x.numerator % P * pow(x.denominator % P, -1, P) % P


def finv(x):
    return 1 / x if P is None else pow(x, -1, P)


def red(x):
    return x if P is None else x % P


def iszero(x):
    return (x == 0) if P is None else (x % P == 0)


def smul(A, B, K):
    """product of two eps-series, truncated at eps^K"""
    C = [0] * (K + 1)
    for i in range(min(len(A), K + 1)):
        x = A[i]
        if x:
            for j in range(min(len(B), K + 1 - i)):
                y = B[j]
                if y:
                    C[i + j] += x * y
    return [red(c) for c in C]


def sinv(A, K):
    """1/A for an eps-series with A[0] != 0"""
    i0 = finv(A[0]); C = [0] * (K + 1); C[0] = i0
    for m in range(1, K + 1):
        s = 0
        for i in range(1, min(m, len(A) - 1) + 1):
            if A[i]:
                s += A[i] * C[m - i]
        C[m] = red(-s * i0)
    return C


def sexp(A, K):
    """exp(A) for a scalar series with A[0] = 0"""
    E = [0] * (K + 1); E[0] = cv(1)
    for m in range(1, K + 1):
        s = 0
        for i in range(1, min(m, len(A) - 1) + 1):
            if A[i]:
                s += i * A[i] * E[m - i]
        E[m] = red(s * finv(cv(m)))
    return E


def spow_T(Th, alpha, J, K):
    """Theta(T)^alpha through T^J, coefficients eps-series (length K+1); Th[0] = 1 (J.C.P. Miller recurrence)"""
    B = [None] * (J + 1); B[0] = [cv(1)] + [0] * K
    for j in range(1, J + 1):
        acc = [0] * (K + 1)
        for k in range(1, j + 1):
            cf = (alpha + 1) * k - j
            if cf == 0:
                continue
            cfv = cv(cf); t = Th[k]; b = B[j - k]
            for i1 in range(K + 1):
                x = t[i1]
                if x:
                    xx = cfv * x
                    for i2 in range(K + 1 - i1):
                        y = b[i2]
                        if y:
                            acc[i1 + i2] += xx * y
        ij = finv(cv(j))
        B[j] = [red(a * ij) for a in acc]
    return B


def gbinom(e, j):
    """generalised binomial coefficient C(e, j), e rational"""
    r = Fr(1)
    for i in range(j):
        r = r * (Fr(e) - i) / (i + 1)
    return r


def engine(loc, N, K=None, verbose=False):
    """loc: M, cs, es (singular factor); bstar[n] (n = 0..N; beta*_n, index 0 unused); Phi[m] (Taylor coefficients of
    the regular phase at s*, m = 0..K+1, only m >= 2 used); Bk[k][r] (k <= N, r <= K); optional 'extra' = list of
    further singular factors {M, cs, es, bstar} at the same point (degenerate tests).
    Returns ct[n] (eps-series of length K+1), n = 0..N."""
    K = N if K is None else K
    M, cs, es = loc['M'], Fr(loc['cs']), Fr(loc['es'])
    Phi, Bk, bstar = loc['Phi'], loc['Bk'], loc['bstar']
    sings = [(Fr(M), cs, es, bstar)] + [(Fr(x['M']), Fr(x['cs']), Fr(x['es']), x['bstar']) for x in loc.get('extra', [])]
    Mc = sum(m_ * c_ for m_, c_, _, _ in sings)          # several singular factors at s* (tests only)
    f = {}
    for k in range(2, 2 * N + 3):
        s = [0] * (K + 1)
        s[0] = cv(Mc * (-1) ** k / (k * (k - 1)))
        for m in range(k, min(K + 1, len(Phi) - 1) + 1):
            s[m - 1] = red(s[m - 1] + cv(comb(m, k)) * cv(Phi[m]))
        f[k] = s
    a = [red(-2 * x) for x in f[2]]
    ainv = sinv(a, K)
    Th = [[cv(1)] + [0] * K] + [smul([red(-2 * x) for x in f[k + 2]], ainv, K) for k in range(1, 2 * N + 1)]
    e0 = sum(m_ * (e_ - Fr(1, 2)) for m_, _, e_, _ in sings)
    dlo, dhi = -N, K
    Lm = []
    for m in range(N + 1):
        Psi = spow_T(Th, Fr(-(2 * m + 1), 2), 2 * m, K)
        row = {}
        for d in range(dlo, dhi + 1):
            acc = [0] * (K + 1)
            for j in range(2 * m + 1):
                bj = cv(gbinom(e0 + d, j))
                if not iszero(bj):
                    ps = Psi[2 * m - j]
                    for i in range(K + 1):
                        if ps[i]:
                            acc[i] += bj * ps[i]
            row[d] = [red(x) for x in acc]
        Lm.append(row)
        if verbose:
            print('   engine: m = %d done' % m, flush=True)
    sl = [0] * (N + 1)
    for n in range(1, N + 1):
        sl[n] = red(sum(cv(m_) * cv(b_[n]) * cv(c_ ** (-n)) for m_, c_, _, b_ in sings))
    S = sexp(sl, N)
    C = sinv([Bk[0][r] for r in range(K + 1)], K)
    w = {}
    for n1 in range(N + 1):
        for d in range(-n1, K - n1 + 1):
            s = 0
            for k in range(0, min(n1, n1 + d) + 1):
                r = d + n1 - k
                if k < len(Bk) and r < len(Bk[k]) and Bk[k][r]:
                    s += S[n1 - k] * Bk[k][r]
            w[(n1, d)] = red(s)
    ct = [[0] * (K + 1) for _ in range(N + 1)]
    apow = [cv(1)] + [0] * K; dfm = 1
    for m in range(N + 1):
        if m > 0:
            apow = smul(apow, ainv, K); dfm *= (2 * m - 1)
        pref = [red(cv(dfm) * x) for x in smul(C, apow, K)]
        for n1 in range(N - m + 1):
            inner = [0] * (K + 1)
            for d in range(-n1, K - n1 + 1):
                wv = w[(n1, d)]
                if wv:
                    sh = d + n1; Ld = Lm[m][d]
                    for i in range(K + 1 - sh):
                        if Ld[i]:
                            inner[sh + i] += wv * Ld[i]
            inner = [red(x) for x in inner]
            add = smul(pref, inner, K)
            ct[m + n1] = [red(x + y) for x, y in zip(ct[m + n1], add)]
    return ct


def laurent_report(ct, N, K=None):
    """returns (first n with a pole or None, list of values c_n(s*) = [eps^n] ct_n, pole table)"""
    K = N if K is None else K
    first = None; vals = []; poles = {}
    for n in range(N + 1):
        pl = [r for r in range(min(n, K + 1)) if not iszero(ct[n][r])]
        if pl:
            poles[n] = pl
            if first is None:
                first = n
        vals.append(ct[n][n] if n <= K else None)
    return first, vals, poles
