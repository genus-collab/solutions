# Fresh exact checks (Fractions) of the hand-derived generic-weight identities used in LC_CY.  Evidence only.
from fractions import Fraction as Fr
from math import isqrt
import random
rnd = random.Random(20261103)
def rq(): return Fr(rnd.randint(-40, 40), rnd.randint(1, 13))
def wts():
    while True:
        a0, a1, a2, e1 = rq(), rq(), rq(), rq()
        e2 = a0 + a1 + a2 - e1
        E = a0 + a1 + a2; A = a0 + a2; B = a0 + a1
        if a1 * a2 * A * B * e1 * e2 * E != 0 and a0 != 0: return a0, a1, a2, e1, e2
fails = 0
def chk(name, ok):
    global fails
    if not ok: fails += 1
    print(('ok   ' if ok else 'FAIL ') + name)

# (1) amplitude identity (3.16)_CY, pointwise in (s, tau, z)
ok = True
for _ in range(30):
    a0, a1, a2, e1, e2 = wts(); E = a0+a1+a2; A = a0+a2; B = a0+a1
    s, tau, z = rq(), rq(), rq()
    if s in (0, 1) or z == 0 or tau * s == 1: continue
    al = 1 + a2/z; be = 1 + A/z; ga = 2 + E/z
    N = B*tau*s*s - (E + a0*tau)*s + A
    lhs = tau*(1-tau)*al*(al+1)*s*s + (ga - (al+be+1)*tau)*al*s*(1-tau*s) - al*be*(1-tau*s)**2
    rhs = -al*((1-2*s)*(1-tau*s) + 2*tau*s*(1-s) + N/z)
    ok &= (lhs == rhs)
chk('(3.16)_CY amplitude identity times (1-tau s)^3', ok)

# (2) operator identity (4.4)_CY on Laurent monomials tau^n
def lp_add(a, b):
    c = dict(a)
    for k, v in b.items(): c[k] = c.get(k, 0) + v
    return {k: v for k, v in c.items() if v != 0}
def lp_mul(a, b):
    c = {}
    for k1, v1 in a.items():
        for k2, v2 in b.items(): c[k1+k2] = c.get(k1+k2, 0) + v1*v2
    return {k: v for k, v in c.items() if v != 0}
def lp_sc(x, a): return {k: x*v for k, v in a.items() if x*v != 0}
def lp_dtau(a): return {k-1: k*v for k, v in a.items() if k != 0}
ok = True
for _ in range(12):
    a0, a1, a2, e1, e2 = wts(); E = a0+a1+a2; A = a0+a2; B = a0+a1
    z = rq()
    if z == 0: continue
    al = 1 + a2/z; be = 1 + A/z; ga = 2 + E/z
    one = {0: Fr(1)}; tau1 = {1: Fr(1)}; onemt = lp_add(one, lp_sc(-1, tau1))
    omega = lp_add({-1: E + z}, {0: -(A + z)})
    S = lp_add({-1: E}, {0: -a0}); Pp = {-1: a1*a2}
    def Dd(f): return lp_sc(z, lp_mul(onemt, lp_dtau(f)))
    def Op(f): return lp_add(Dd(f), lp_mul(omega, f))
    def Gauss(f):
        return lp_add(lp_add(lp_mul(lp_mul(tau1, onemt), lp_dtau(lp_dtau(f))), lp_mul(lp_add({0: ga}, {1: -(al+be+1)}), lp_dtau(f))), lp_sc(-al*be, f))
    for n in range(-3, 7):
        f = {n: Fr(1)}
        lhs = lp_add(lp_add(Op(Op(f)), lp_sc(-1, lp_mul(S, Op(f)))), lp_mul(Pp, f))
        rhs = lp_mul(lp_sc(z*z, lp_mul(onemt, {-1: Fr(1)})), Gauss(f))
        ok &= (lp_add(lhs, lp_sc(-1, rhs)) == {})
chk('(4.4)_CY operator identity on tau^n, n=-3..6', ok)

# (3) critical-point identities at points where delta^2 is a rational square (both roots)
res = {}
def rec(name, val): res[name] = res.get(name, True) and val
cnt = 0; tries = 0
while cnt < 25 and tries < 50000:
    tries += 1
    a0, a1, a2, e1, e2 = wts(); E = a0+a1+a2; A = a0+a2; B = a0+a1
    tau = rq()
    if tau in (0, 1): continue
    d2 = (E + a0*tau)**2 - 4*A*B*tau
    if d2 <= 0: continue
    pn, pd = d2.numerator, d2.denominator
    rn, rd = isqrt(pn), isqrt(pd)
    if rn*rn != pn or rd*rd != pd: continue
    delta = Fr(rn, rd)
    for sgn in (1, -1):
        dl = sgn*delta
        s = (E + a0*tau - dl)/(2*B*tau)
        N = B*tau*s*s - (E + a0*tau)*s + A
        rec('N(s) = 0', N == 0)
        U = 1 - tau*s
        if U == 0 or s in (0, 1): continue
        Q = 1/(1-tau)
        L1 = A/(tau*s) - a0; L2 = B*s - a0
        rec('L = -a1 Q U/(1-QU) = A/(tau s) - a0', -a1*Q*U/(1-Q*U) == L1)
        rec('L1+L2 = E/tau - a0', L1 + L2 == E/tau - a0)
        rec('L1 L2 = a1 a2/tau', L1*L2 == a1*a2/tau)
        rec('tau(L1-L2) = delta', tau*(L1-L2) == dl)
        rec('a1a2 - E L2 = -B tau s L2', a1*a2 - E*L2 == -B*tau*s*L2)
        rec('(1-tau) a2 s/(1-tau s) = A - B s', (1-tau)*a2*s/(1-tau*s) == A - B*s)
        f2 = -A/s**2 - a1/(1-s)**2 + a2*tau**2/(1-tau*s)**2
        rec("f''(s*) = -delta/(s(1-s)U)", f2 == -dl/(s*(1-s)*U))
        W2 = -a1*Q*U/(1-Q*U)**2 - A*U/(1-U)**2
        rec("W_2(U) = U^2 f''(s)/tau^2", W2 == U*U*f2/tau**2)
        Delta = e1*e2*(a0*a1*a2 + 2*a1*a2*L1 - E*L1*L1)
        rec('Delta = e1e2(L1-L2)(a1a2-E L1)', Delta == e1*e2*(L1-L2)*(a1*a2-E*L1))
        rec('Delta = a1 A e1e2 W_2', Delta == a1*A*e1*e2*W2)
        rec('Delta_1 = -e1e2 E^3 alg^2/tau^2 with alg^2 = a1 A delta U/(E^3 s(1-s))', Delta == -e1*e2*E**3*(a1*A*dl*U/(E**3*s*(1-s)))/tau**2)
        p = (A - B*s)/(1-tau)
        Eik = tau*(1-tau)*p*p + (E - (a2+A)*tau)*p - a2*A
        rec('eikonal (1-tau) Eik(dphi_tau) = B N(s)', (1-tau)*Eik == B*N)
        Ptr = 2*tau*(1-tau)*p + E - (a2+A)*tau
        rec("P_tr = -N'(s) (= delta on the s* branch)", Ptr == -(2*B*tau*s - (E + a0*tau)))
        def pair(a, b): return e1*e2*(a0*a1*a2*a[0]*b[0] + a1*a2*(a[0]*b[1] + a[1]*b[0]) - E*a[1]*b[1])
        rec('<chi1,chi2> = 0', pair((1, L1), (1, L2)) == 0)
        rec('<chi1,chi1> = Delta_1', pair((1, L1), (1, L1)) == Delta)
        Delta2 = e1*e2*(L2-L1)*(a1*a2-E*L2)
        rec('Theta = Delta_2/(L2-L1) = e1e2(a1a2-E L2) = -e1e2 B tau s L2', Delta2/(L2-L1) == e1*e2*(a1*a2-E*L2) == -e1*e2*B*tau*s*L2)
        rec('Delta_2 = e1e2 B delta s* L_2', Delta2 == e1*e2*B*dl*s*L2)
        def Aop(a): return (a[1], -(a1*a2/tau)*a[0] + (E/tau - a0)*a[1])
        aa = (rq(), rq()); bb = (rq(), rq())
        rec('A self-adjoint', pair(Aop(aa), bb) == pair(aa, Aop(bb)))
        Theta = e1*e2*(a1*a2-E*L2)
        lhs = pair((Aop(aa)[0]-L1*aa[0], Aop(aa)[1]-L1*aa[1]), bb)
        rec('(4.8)_CY: <(A-L1)a,b> = pi(a) pi(b)/Theta', lhs == pair(aa, (1, L2))*pair(bb, (1, L2))/Theta)
        rec("L1 - dphi' = E/tau - A", L1 - (A - B*s) == E/tau - A)
        rec('tau p_1 = 1/s, tau p_Q = -(1-tau)/(1-s)', tau/(1-U) == 1/s and tau/(1-Q*U) == -(1-tau)/(1-s))
        s_tau = s*L2/dl
        rec('s_tau = s L2/delta (implicit differentiation)', B*s*s + 2*B*tau*s*s_tau - a0*s - (E + a0*tau)*s_tau == 0)
        p_tau = (-B*s_tau*(1-tau) + (A - B*s))/(1-tau)**2
        Rtr = tau*(1-tau)*p_tau + (2-3*tau)*p - (a2+A)
        rec('R_tr = a0 - 2Bs - B tau s^2 L2/delta', Rtr == a0 - 2*B*s - B*tau*s*s*L2/dl)
        u = (E-2*a1)*(E-2*a2)
        dl_tau = (-(E*E+u) + 2*a0*a0*tau)/(2*dl)
        U_tau = -s - tau*s_tau
        lg = dl_tau/dl + U_tau/U - s_tau/s + s_tau/(1-s)
        rec('(log alg^2)_tau = 2 R_tr/delta', lg == 2*Rtr/dl)
        rec('u + E^2 = 2 a0 E + 4 a1 a2', u + E*E == 2*a0*E + 4*a1*a2)
        rec('delta^2 = E^2 - (E^2+u) tau + a0^2 tau^2', dl*dl == E*E - (E*E+u)*tau + a0*a0*tau*tau)
    cnt += 1
print('critical-point identities tested at', cnt, 'rational points (both roots):')
for k, v in res.items(): chk('  ' + k, v)

# (4) metric
ok = True
for _ in range(10):
    a0, a1, a2, e1, e2 = wts(); E = a0+a1+a2; A = a0+a2; B = a0+a1
    eta = [[E/(a1*a2*A*B*e1*e2), 1/(A*B*e1*e2)], [1/(A*B*e1*e2), -a0/(A*B*e1*e2)]]
    ei = [[e1*e2*a0*a1*a2, e1*e2*a1*a2], [e1*e2*a1*a2, -e1*e2*E]]
    prod = [[sum(eta[i][k]*ei[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    ok &= (prod == [[1, 0], [0, 1]])
chk('eta^-1 = e1 e2 [[a0 a1 a2, a1 a2],[a1 a2, -E]]', ok)

# (5) constant line: power sums
ok = True
for _ in range(10):
    a0, a1, a2, e1, e2 = wts(); E = a0+a1+a2; A = a0+a2
    for k in range(1, 5):
        m = 1 - 2*k
        ok &= (a1**m + A**m + (-e1)**m + (-e2)**m + E**m + (-a1)**m + (-A)**m == E**m + (-e1)**m + (-e2)**m)
        ok &= (a1**m + A**m + (-e1)**m + (-e2)**m - (E**m + (-e1)**m + (-e2)**m) == a1**m + A**m + (-E)**m)
chk('h_nw h_(E,-a1,-A) = h_(E,-e1,-e2) and h_nw h_w^-1 = h_(a1,A,-E) (power sums k<=4)', ok)

# (6) first-order hand values
ok = True
for _ in range(10):
    a0, a1, a2, e1, e2 = wts(); E = a0+a1+a2; A = a0+a2; B = a0+a1; u = (E-2*a1)*(E-2*a2)
    ok &= (2*E*(a2+A) - 4*a2*A - E*E == u)
    l0 = E/2 - a0 - a1*a2/E
    ok &= (u/2 + l0*E == a1*a2)
    ok &= (l0 == E/2 - a0/2 - (u + E*E)/(4*E))
    s1 = A*a1*a2/E**3
    ok &= ((A/E)*(B*A/E - a0)/E == s1)
    dl1 = -(E*E+u)/(2*E)
    ok &= (dl1*(A/E)*(a1*a2/E) + E*s1*(a1*a2/E) + E*(A/E)*B*s1 == 0)
    ok &= (B*A/E - a0 == a1*a2/E)
    ok &= (A*B == a1*a2 + a0*E)
    ok &= (-(E*E+u)/(2*E*E) - (2*a1*a2 + a0*E)/(E*E) == -(E*E+u)/(E*E))
chk("first-order hand values: u, l_0, H-row a1a2, s*'(0), Delta_2'(0)=0, alg_1", ok)
print('FAILURES:', fails)
