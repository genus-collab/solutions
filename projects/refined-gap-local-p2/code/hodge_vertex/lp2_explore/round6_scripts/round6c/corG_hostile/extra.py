from gsum import *
from test_lr import solve_mod
from math import pi, sqrt, log, lgamma, exp
F = FpField()
# (a) string and vanishing, exact
e1, e2 = Fr(3), Fr(7)
E = e1 + e2
print('C_P2: Xi_2(1)', Xi(2, (1,), e1, e2), ' Xi_2(1,1)', Xi(2, (1, 1), e1, e2), ' Xi_3(1)', Xi(3, (1,), e1, e2))
for R in [(3,), (2, 2), (4,)]:
    J = sum(R)
    print('S_P2: Xi_2(1,R)/Xi_2(R) with R=', R, Xi(2, (1,) + R, e1, e2) / Xi(2, R, e1, e2) if Xi(2, R, e1, e2) != 0 else None, ' expected J/E^2 =', Fr(J) / E ** 2)
# (b) edge at t_c=0
D, su = conifold_data(F, e1, e2, 22, 8, 'hodge')
rho_ok = True
vals = []
for i in range(3):
    for j in range(3):
        for a in range(4):
            for b in range(4):
                x = D.edge(i, j, a, b)
                vals.append(x.val if len(x.c) else 99)
print('min t_c valuation of all edges a,b<=3 over colours:', min(vals))
ok = True
for a in range(4):
    for b in range(4):
        x = D.edge(0, 0, a, b)
        exp_ = -(e1 * e2 * E ** 4 / 27) * E ** (-a - 1) * E ** (-b - 1)
        ok = ok and F.is0(F.red(x.coeff(0) - F.fr(exp_)))
print('L3: Ehat_00 at t_c=0 = -(e1 e2 E^4/27) lam1 lam1 coefficientwise a,b<=3:', ok)
# (c) c1-independence at fixed FF
su0 = conifold_setup(F, 22)
t = su0['t']
vs = []
om = F.one
for i in range(3):
    ti = t * om
    one = L.const(F, F.one, 24)
    vs.append((ti * F.fr(3) + one) * (ti * F.fr(E)).inv())
    om = F.red(om * F.omega)
c1b = su0['c1'] * (L.const(F, F.one, 24) + su0['xi'] * F.fr(Fr(2, 7))) * F.fr(5)
D2 = Data(F, e1, e2, vs, su0['xi'], c1b, su0['FF'], 8, 'hodge')
D1 = Data(F, e1, e2, vs, su0['xi'], su0['c1'], su0['FF'], 8, 'hodge')
d = D1.total(2) - D2.total(2)
print('c1-independence of G_2 at fixed F (all coefficients zero):', len(d.c) == 0, 'prec', d.prec)
d = D1.total(3) - D2.total(3)
print('c1-independence of G_3 at fixed F (all coefficients zero):', len(d.c) == 0, 'prec', d.prec)
# (d) coefficient extraction at the conifold
for g in (2, 3):
    samples = [(Fr(3), Fr(7)), (Fr(1), Fr(1)), (Fr(5, 2), Fr(-3, 7)), (Fr(2), Fr(-5, 3))][:g + 1]
    tots = []
    for (x1, x2) in samples:
        Dd, su = conifold_data(F, x1, x2, 26, needed_K(g), 'hodge')
        tots.append((Dd.total(g), su))
    tc = tots[0][1]['tc']
    tci = tc.inv()
    # expand in t_c^(-n) : subtract
    A = []
    for (x1, x2) in samples:
        EE = x1 + x2
        A.append([F.fr(EE ** (2 * g - 2 * j) * (x1 * x2) ** j) for j in range(g + 1)])
    for n in range(2 - 2 * g, 1):
        b = [F.red(F.fr(x1 * x2) * tt.coeff(n)) for (x1, x2), (tt, _) in zip(samples, tots)]
        f = solve_mod(A, b, F.p)
        print('g', g, 'xi^%d coefficients of (F^cf_{g,0},...,F^cf_{0,g}):' % n, [str(rr(x)) for x in f])
    pw = L.const(F, F.one, 24)
    for _ in range(2 * g - 2):
        pw = pw * tci
    print(' 1/t_c^(2g-2) in xi:', [str(rr(pw.coeff(n))) for n in range(2 - 2 * g, 1)])
    bk = []
    for j in range(g + 1):
        # coefficient of E^(2g-2j) (e1 e2)^(j-1) in (-3)^(g-1) beta
        pass
    b = [F.red(F.fr(x1 * x2 * Fr(-3) ** (g - 1) * beta(2 * g - 2, x1, x2))) for (x1, x2) in samples]
    print(' b_{k,h} (k=g..0) from beta:', [str(rr(x)) for x in solve_mod(A, b, F.p)])
# (e) t_c by long series against BS expression
def tc_num(xi, M=4000):
    s = 0.0
    w = 1.0
    pc = 0.0
    p = xi
    for n in range(M):
        pc += w
        s += pc / (n + 1) * p
        p *= xi
        w = w * (1 / 3 + n) * (2 / 3 + n) / (n + 1) ** 2
    return s
def bs(y):
    l = log(-y)
    g = 0.0
    gg = 0.0
    H1 = 0.0
    for d in range(1, 600):
        lg = lgamma(3 * d) - 3 * lgamma(d + 1) + d * log(-y)
        if lg < -50:
            break
        t = exp(lg)
        H1 = sum(1.0 / k for k in range(d + 1, 3 * d))
        g += 3 * t
        gg += 9 * t * H1
    return sqrt(3) / (2 * pi) * (l * l / 2 + g * l + gg - pi * pi / 2)
for y in (-1 / 54, -1 / 81, -1 / 40):
    print('t_c long series', tc_num(1 + 27 * y), ' BS expression, log y = l+i pi', bs(y))
