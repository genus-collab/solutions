# Check 5: dilaton calculus, uniform Lemma DC, by exact series in u, independent of the MV code path.
# (a) operator identity on e^{Ju} w^n in Q[[u]]; (b) genus 0: (J+D)^k G_U = e^{uJ} w^{r-2+k} (falling(delta,k) X_0(U) + O(w)) as Laurent polynomial in w, verified as series in u;
# (c) empty-profile table = uniform lemma at U=1^v, recomputed directly from G^(0)_empty = Sum_{a>=1} a^(a-3) y^a/a!, and genus 1 from y d_y G^(1)_empty;
# (d) controls: k=0 v=0 functions are more singular, as the referees reported (D G_0 etc.), so v>=1 is essential.
# All series are in u to order NN. y = u e^{-u}; substitution by composing power series.
from fractions import Fraction as Fr
from math import factorial
from tl_lib import *
NN = 14

def mul(a, b): return ser_mul(a, b, NN)
def exp_mu(m):  # e^{m u}
    return [Fr(m)**i / factorial(i) if (m != 0 or i == 0) else Fr(0) for i in range(NN + 1)]
one = [Fr(1)] + [Fr(0)] * NN
def inv_1mu_pow(n):  # w^n = (1-u)^{-n}, n integer
    r = list(one)
    if n >= 0:
        f = [Fr(1)] * (NN + 1)
    else:
        f = [Fr(1), Fr(-1)] + [Fr(0)] * (NN - 1)
    for _ in range(abs(n)): r = mul(r, f)
    return r
yu = mul([Fr(0), Fr(1)] + [Fr(0)] * (NN - 1), exp_mu(-1))
ypow = [one]
for a in range(1, NN + 1): ypow.append(mul(ypow[-1], yu))
def comp(coefs):
    """sum_a coefs[a] y^a with y=u e^{-u}, as series in u."""
    r = [Fr(0)] * (NN + 1)
    for a, c in enumerate(coefs):
        if a > NN: break
        for i in range(NN + 1): r[i] += c * ypow[a][i]
    return r
def lp(terms, m):
    """e^{m u} * Sum_n c_n w^n, terms dict n->c, as a series in u."""
    r = [Fr(0)] * (NN + 1)
    for n, c in terms.items():
        t = inv_1mu_pow(n)
        for i in range(NN + 1): r[i] += c * t[i]
    return mul(r, exp_mu(m))
def opJD(terms, J):
    """Lemma 2.2: e^{-uJ}(J+D)e^{uJ} w^n = (J-n) w^{n+1} + n w^{n+2}."""
    out = {}
    for n, c in terms.items():
        out[n + 1] = out.get(n + 1, 0) + c * (J - n)
        out[n + 2] = out.get(n + 2, 0) + c * n
    return out
def opdy(terms, J):
    """e^{-u(J+1)} d_y e^{uJ} w^n = J w^{n+1} + n w^{n+2}."""
    out = {}
    for n, c in terms.items():
        out[n + 1] = out.get(n + 1, 0) + c * J
        out[n + 2] = out.get(n + 2, 0) + c * n
    return out

def base_g0(U):
    """e^{-uJ} G^(0)_U as Laurent polynomial in w, from the A_m recursion."""
    r = len(U); J = sum(U)
    if r == 1:
        return {-1: Fr(1, J * (J + 1)), 0: Fr(1, J * J * (J + 1))}
    if r == 2:
        return {0: Fr(1, J)}
    t = {0: Fr(1, J)}
    for _ in range(r - 2):
        t = opJD(t, J)
    return {n: c for n, c in t.items() if c != 0}

def series_g0(U, k):
    """direct: Sum_a (J+a)^(r+a-3+k) y^a / a!"""
    r = len(U); J = sum(U)
    coefs = []
    for a in range(NN + 1):
        e = r + a - 3 + k
        coefs.append((Fr(J + a)**e) / factorial(a))
    return comp(coefs)

cnt = bad = 0
for U in [(1,), (2,), (3,), (5,), (1, 1), (2, 1), (3, 2), (1, 1, 1), (2, 1, 1), (3, 3, 2), (1, 1, 1, 1), (4, 2, 1, 1), (2, 2, 2, 2, 2)]:
    J = sum(U); r = len(U)
    t = base_g0(U)
    d0 = (J - r) + 2
    X0num = t.get(r - 2)
    # Gamma(J)/Gamma(J-r+3)
    from mvx import x0, falling
    assert X0num == x0(U), (U, X0num, x0(U))
    assert min(t) == r - 2
    for k in range(0, 6):
        s = series_g0(U, k)
        l = lp(t, J)
        cnt += 1
        if s != l:
            bad += 1; print("DC g0 fail", U, k)
        lead = t.get(r - 2 + k, 0)
        if lead != falling(d0, k) * x0(U) or (t and min(t) < r - 2 + k):
            bad += 1; print("lead fail", U, k, lead)
        t = opJD(t, J)
print(f"(a,b) genus-0 lemma: {cnt} series comparisons in u through u^{NN}, {bad} failures")

# (c) table
sigma = Fr(18, 49)
print("(c) empty-profile table, recomputed from G^(0)_empty directly and from y d_y G^(1)_empty = (w-1)(w-1-1/sigma)/24:")
# genus 0: series of d_y^v D^k G^(0)_empty: coefficient of y^(a-v) with a>=1 : a^(a-3+k) * a(a-1)..(a-v+1)/a!
def g0_empty_series(v, k):
    coefs = [Fr(0)] * (NN + 1)
    for a in range(1, NN + v + 1):
        if a - v < 0: continue
        if a - v > NN: break
        c = Fr(a)**(a - 3 + k) / factorial(a)
        for i in range(v): c *= (a - i)
        coefs[a - v] = c
    return comp(coefs)
okc = 0; tot = 0
for v in range(1, 6):
    # uniform lemma: U=1^v
    U = (1,) * v
    t = base_g0(U)
    for k in range(0, 5):
        tot += 1
        s = g0_empty_series(v, k)
        l = lp(t, v)
        lead = t.get(v - 2 + k, Fr(0))
        table = factorial(v - 1) * (Fr(1, 2) if k == 0 else (Fr(1) if k in (1, 2) else Fr(0)))
        if s == l and lead == table: okc += 1
        else: print("   table mismatch g=0", v, k, lead, table)
        t = opJD(t, v)
print(f"   genus 0: {okc}/{tot} entries (v<=5,k<=4) agree with (v-1)! * (1/2,1,1,0,0)")
# genus 1
okc = 0; tot = 0
X11 = -(1 + 1 / sigma) / 24
for v in range(1, 6):
    # d_y G^(1)_empty = e^u w (w-1-1/sigma)/24
    t = {1: -(1 + 1 / sigma) / 24, 2: Fr(1, 24)}
    J = 1
    for i in range(v - 1):
        t = opdy(t, J); J += 1
    assert J == v
    # k-loop: uniform: order 2g-2+r+k = v+k
    for k in range(0, 5):
        tot += 1
        lead = t.get(v + k, Fr(0))
        table = factorial(v - 1) * X11 if k == 0 else Fr(0)
        if lead == table and min(n for n, c in t.items() if c != 0) >= v + k: okc += 1
        else: print("   table mismatch g=1", v, k, lead, table)
        t = opJD(t, v)
print(f"   genus 1: {okc}/{tot} entries agree with ((v-1)! X_1(1), 0, 0, ...)")
# (d) v=0 functions
print("(d) v=0 is excluded: D G^(0)_empty =? 1/2 - w^-2/2:", comp([Fr(0)] + [Fr(a)**(a - 2) / factorial(a) for a in range(1, NN + 1)]) == lp({0: Fr(1, 2), -2: Fr(-1, 2)}, 0),
      " (order -2, while the nominal order 2g-2+k+v at (g,k,v)=(0,1,0) would be -1)")
