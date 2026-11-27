# Check 1: the fermionic bilinear identity and its members P_k, generic content-product tau functions.
# (a) full bilinear identity [z^0] e^{xi(t-t2,z)} tau_0(t-[1/z]) tau_1(t2+[1/z]) = tau_1(t) tau_0(t2), level by level, random t,t2,p',B
# (b) P_k tau_0.tau_1 = 0, k=2..kmax, all monomials
# (c) controls: swapped pair, tau_{-1}, sign flip of the D_k coefficient, printed exponent sign of Harnad-Orlov (2.2)
# (d) P_k in the span of the Jimbo-Miwa Appendix 1 table entries (degree 2..5)
import sys, random
from fractions import Fraction as Fr
from tl_lib import *

d = int(sys.argv[1]) if len(sys.argv) > 1 else 8
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 7
random.seed(seed)
def rr(): return Fr(random.randint(-9, 9) or 1, random.randint(1, 7))
B = rr() if len(sys.argv) < 4 else Fr(sys.argv[3])
if B in (0, 1, -1): B = Fr(3, 2)
pprime = {k: rr() for k in range(1, d + 1)}
print("d =", d, "B =", B)

allnu = [nu for n in range(d + 1) for nu in partitions(n)]
snu = {nu: schur_in_p(nu) for nu in allnu}
snuprime = {nu: sum((c * eval_p(mu, pprime) for mu, c in snu[nu].items()), Fr(0)) if nu else Fr(1) for nu in allnu}
snu[()] = {(): Fr(1)}

def tau(n):
    out = {}
    for nu in allnu:
        w = B**(f2(nu) + n * sum(nu)) * snuprime[nu]
        if w == 0: continue
        for mu, c in snu[nu].items():
            out[mu] = out.get(mu, 0) + w * c
    return out

t0, t1, tm = tau(0), tau(1), tau(-1)

# (b) members
def apply_member(P, f, g, deg):
    tot = {}
    for alpha, c in P.items():
        h = hirota(f, g, list(alpha), deg)
        for m, v in h.items():
            tot[m] = tot.get(m, 0) + c * v
    return tot

kmax = min(d - 1, 7)
for k in range(2, kmax + 1):
    P = member_Pk(k)
    deg = d - k
    res = apply_member(P, t0, t1, d)
    # only monomials of degree <= d-k are exact under truncation at d
    bad = sum(1 for m, v in res.items() if sum(m) <= deg and v != 0)
    tot = sum(1 for m, v in res.items() if sum(m) <= deg)
    cs = apply_member(P, t1, t0, d)
    csw = sum(1 for m, v in cs.items() if sum(m) <= deg and v != 0)
    cm = apply_member(P, t0, tm, d)
    cmm = sum(1 for m, v in cm.items() if sum(m) <= deg and v != 0)
    Pw = dict(P); Pw[(k,)] = Pw[(k,)] - 2
    cw = apply_member(Pw, t0, t1, d)
    cwn = sum(1 for m, v in cw.items() if sum(m) <= deg and v != 0)
    print(f"P_{k}: monomials {tot}, nonzero {bad} | controls nonzero: swap {csw}, tau_-1 {cmm}, wrong D_k coeff {cwn}")

# (a) full bilinear identity level by level. Variables: t, t2 random rationals; zeta = 1/z.
tv = {k: rr() for k in range(1, d + 1)}
t2v = {k: rr() for k in range(1, d + 1)}

def elem_schur(tt, L):
    # p_j(t): exp(sum t_k z^k)
    a = [Fr(0)] + [tt.get(k, Fr(0)) for k in range(1, L + 1)]
    return ser_exp(a, L)

def poly_shift(nu, tt, sign, L):
    """s_nu(t + sign*[zeta]) as polynomial in zeta, list up to degree L; p_k = k t_k + sign*zeta^k."""
    r = [Fr(0)] * (L + 1)
    for mu, c in snu[nu].items():
        cur = [Fr(1)] + [Fr(0)] * L
        for k in mu:
            fac = [Fr(0)] * (L + 1)
            fac[0] = k * tt[k]
            if k <= L: fac[k] = Fr(sign)
            cur = ser_mul(cur, fac, L)
        for i in range(L + 1):
            r[i] += c * cur[i]
    return r

def wt(nu, n): return B**(f2(nu) + n * sum(nu)) * snuprime[nu]

def check_full(sign_exp, n_a, n_b, n_ra, n_rb, sign_a=-1, sign_b=+1):
    """level m: sum_j p_j(sign_exp*(t-t2)) * [zeta^j] sum_{|nu|+|nu'|=m} w_nu w_nu' s_nu(t+sign_a zeta) s_nu'(t2+sign_b zeta)
       against sum_{|nu|+|nu'|=m} w s_nu(t) s_nu'(t2) (charges n_ra at t, n_rb at t2)."""
    diff = {k: sign_exp * (tv[k] - t2v[k]) for k in tv}
    out = []
    for m in range(0, d + 1):
        pj = elem_schur(diff, m)
        lhs = Fr(0)
        for n1 in range(m + 1):
            n2 = m - n1
            for nu in partitions(n1):
                a = poly_shift(nu, tv, sign_a, m)
                wa = wt(nu, n_a)
                for nu2 in partitions(n2):
                    b = poly_shift(nu2, t2v, sign_b, m)
                    pr = ser_mul(a, b, m)
                    w = wa * wt(nu2, n_b)
                    lhs += w * sum(pj[j] * pr[j] for j in range(m + 1))
        rhs = Fr(0)
        for n1 in range(m + 1):
            n2 = m - n1
            for nu in partitions(n1):
                for nu2 in partitions(n2):
                    rhs += wt(nu, n_ra) * wt(nu2, n_rb) * sum((c * eval_p(mu, {k: k * tv[k] for k in tv}) for mu, c in snu[nu].items()), Fr(0)) * \
                           sum((c * eval_p(mu, {k: k * t2v[k] for k in t2v}) for mu, c in snu[nu2].items()), Fr(0))
        out.append(lhs - rhs)
    return out

L = min(d, 6)
save = d
d = L
r = check_full(+1, 0, 1, 1, 0)
print("full identity [z^0] e^{xi(t-t2,z)} tau_0(t-[1/z]) tau_1(t2+[1/z]) - tau_1(t) tau_0(t2), levels 0..%d:" % L, ["0" if x == 0 else "NZ" for x in r])
r = check_full(+1, 1, 0, 0, 1)
print("control, charges (1,0) with same right side form:", ["0" if x == 0 else "NZ" for x in r])
r = check_full(-1, 0, 1, 1, 0)
print("control, opposite exponent sign:", ["0" if x == 0 else "NZ" for x in r])
# Harnad-Orlov (2.2) as printed, at ds=0, N=1, N'=0: t-slot of N has +[1/z], exponent -xi(dt,z), dt = t' - t, slot N' at t' with -[1/z]
# in our variables t := position of slot N (=1), t2 := t+dt: oint z^{-1} e^{-xi(t2-t,z)} tau_1(t+[1/z]) tau_0(t2-[1/z]) = tau_0(t) tau_1(t2)
r = check_full(+1, 1, 0, 0, 1, sign_a=+1, sign_b=-1)  # exponent e^{+xi(t-t2,z)} = e^{-xi(dt,z)}
print("Harnad-Orlov (2.2) as printed (exponent e^{-xi(dt,z)}, ds=0, N=1, N'=0):", ["0" if x == 0 else "NZ" for x in r])
r = check_full(-1, 1, 0, 0, 1, sign_a=+1, sign_b=-1)  # exponent e^{+xi(dt,z)}
print("Same with exponent e^{+xi(dt,z)} (the form our derivation gives):", ["0" if x == 0 else "NZ" for x in r])
d = save

# (d) Jimbo-Miwa Appendix 1 table, 1st modified KP, P(D) tau_n . tau_{n+1} = 0, degrees 2..5, transcribed from the scan p.996
def mon(*idx): return tuple(sorted(idx, reverse=True))
JM = {
 2: [{mon(1,1):1, mon(2):1}],
 3: [{mon(1,1,1):1, mon(3):-4, mon(2,1):-3}],
 4: [{mon(1,1,1,1):1, mon(3,1):8, mon(2,1,1):3, mon(2,2):-6},
     {mon(2,1,1):-1, mon(2,2):1, mon(4):2}],
 5: [{mon(1,1,1,1,1):1, mon(5):-16, mon(2,2,1):5, mon(4,1):-10},
     {mon(2,1,1,1):1, mon(3,2):-4, mon(2,2,1):3, mon(4,1):6},
     {mon(1,1,1,1,1):1, mon(3,1,1):-4, mon(2,1,1,1):3, mon(2,2,1):6}],
}
def rank_solve(vectors, target, keys):
    # solve target = sum c_i vectors_i over Q; return coefficients or None
    n = len(vectors)
    M = [[Fr(v.get(k, 0)) for v in vectors] + [Fr(target.get(k, 0))] for k in keys]
    piv = []
    row = 0
    for col in range(n):
        p = next((r for r in range(row, len(M)) if M[r][col] != 0), None)
        if p is None: continue
        M[row], M[p] = M[p], M[row]
        M[row] = [x / M[row][col] for x in M[row]]
        for r in range(len(M)):
            if r != row and M[r][col] != 0:
                f = M[r][col]
                M[r] = [a - f * b for a, b in zip(M[r], M[row])]
        piv.append(col)
        row += 1
    if any(all(x == 0 for x in r[:-1]) and r[-1] != 0 for r in M):
        return None
    sol = [Fr(0)] * n
    for i, c in enumerate(piv):
        sol[c] = M[i][-1]
    return sol
for k in range(2, 6):
    P = member_Pk(k)
    keys = list(partitions(k))
    print(f"P_{k} =", {a: str(c) for a, c in P.items()})
    print("   coefficients on the Jimbo-Miwa degree-%d list:" % k, [str(x) if x is not None else None for x in (rank_solve(JM[k], P, keys) or [None])])
