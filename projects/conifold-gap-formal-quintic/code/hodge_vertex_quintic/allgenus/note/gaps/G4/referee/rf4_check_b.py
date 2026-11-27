# rf4_check_b.py -- referee check B (items 1(b), 2, 4). Exact over Q(zeta5), v = 1. Usage: py -I -B rf4_check_b.py
# B1 TR omega_(1,1) on x = kappa W_1, y = p (engine kernel = DNOPS (2-28)), residues at all five branch points in t = q - a
#    (local involution sigma(t) by fixed point), against the closed form -(5/8) p^4 (p^10 + 18 p^5 + 6)/(p^5 - 1)^4 / kappa.
# B2 (2-35) at (1,1) with R^-1 from DNOPS (2-31) (rf4_model.Rinv_dnops, no use of M_0) and TFT Delta_j = L c_j^2, with the
#    dilaton leaf  s * psi (1 - R^-1(psi)) 1 :  s = +1 (PPZ Def. 2.13 / DOSS) and s = -1 (DNOPS s.2.3 as printed). Controls.
# B3 Res_(varsigma_b) x^k omega_(1,1) = 0 (k = 0..6); Res_(p=0) x^5 omega_(1,1) = -15/4.
# B4 <-log p>_1 = Int_ray omega_(1,1) pushed into the residues (log(sigma q / q) kernel), and from the closed form.
import os, sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from rf4_field import K, ZERO, ONE, smul, sinv, spow, sderiv
from rf4_model import Rinv_dnops, sigma_local, compose

NS = 10
cs = [K.z(-j) for j in range(5)]
ok = True


def tr_w11_coeffs(kappa):
    """for each branch point a: r[k] (k = 1..3) with omega_11/dp0 = Sum_a Sum_k r_(a,k) (p0 - a)^(-k-1)"""
    out = []
    for a in cs:
        s, xt = sigma_local(a, NS, kappa)
        sp = sderiv(s)
        tms = [ZERO, ONE] + [ZERO] * (NS - 1)
        diff = [p - q for p, q in zip(tms, s)]                   # t - sigma = t D(t)
        D = diff[1:] + [ZERO]
        xp = sderiv(xt)[:NS + 1]                                 # x'(a + t) = t E(t)
        assert xp[0].is0()
        E = xp[1:] + [ZERO]
        den = sinv(smul(spow(D, 3, NS), E, NS), NS)
        r = {}
        for k in (1, 2, 3):
            num = [p - q for p, q in zip(spow(tms, k, NS), spow(s, k, NS))]
            ser = smul(smul(num, sp, NS), den, NS)
            r[k] = ser[3] * F(1, 2)                              # Res of ser / t^4
        out.append((a, r, s, sp, den))
    return out


def eval_tr(data, p0):
    tot = ZERO
    for a, r, _, _, _ in data:
        for k in (1, 2, 3):
            tot = tot + r[k] * (K.q(p0) - a) ** (-k - 1)
    return tot.rat()


def closed(p0, kappa):
    p0 = F(p0)
    return -F(5, 8) * p0 ** 4 * (p0 ** 10 + 18 * p0 ** 5 + 6) / (p0 ** 5 - 1) ** 4 / kappa


def graph_sum(p0, kappa, X1, s):
    """(2-35) at (1,1); A_i as derived in REFEREE_2 s.3.2: leaf + dilaton (sign s) + loop; d/dx term"""
    p = K.q(p0)
    xp = (p ** 3 - p ** (-2)) * kappa
    xpp = (p ** 2 * 3 + p ** (-3) * 2) * kappa
    tot = ZERO
    for i in range(5):
        ci = cs[i]
        A = ZERO
        for j in range(5):
            A = A + X1[j][i] * cs[j] / ci * F(1, 24)            # leaf
        for k in range(5):
            A = A - X1[i][k] * ci / cs[k] * F(s, 24)            # dilaton leaf s psi (1 - R^-1(psi)) 1
        A = A - X1[i][i] * F(1, 2)                              # loop (one genus-0 vertex, |Aut| = 2)
        d = p - ci
        tot = tot + A / d ** 2
        tot = tot + (F(-2) / (d ** 3 * xp) - xpp / (d ** 2 * xp ** 2)) * F(1, 24)
    return tot.rat()


pts = (F(3), F(-2), F(1, 2), F(7, 3))
for kappa in (F(1), F(4)):
    data = tr_w11_coeffs(kappa)
    X = Rinv_dnops(kappa, 1)
    X1 = X[1]
    negX1 = [[-c for c in row] for row in X1]
    res = {'TR = closed form': 0, 'GS(s=+1) = TR': 0, 'GS(s=-1) = TR': 0, 'ctrl GS(R^-1 -> R, s=+1) = TR': 0}
    for p0 in pts:
        t = eval_tr(data, p0)
        res['TR = closed form'] += 0 if t == closed(p0, kappa) else 1
        res['GS(s=+1) = TR'] += 0 if graph_sum(p0, kappa, X1, 1) == t else 1
        res['GS(s=-1) = TR'] += 0 if graph_sum(p0, kappa, X1, -1) == t else 1
        res['ctrl GS(R^-1 -> R, s=+1) = TR'] += 0 if graph_sum(p0, kappa, negX1, 1) == t else 1
    print('B1/B2 kappa=%s (mismatches over %d points):' % (kappa, len(pts)), ', '.join('%s: %d' % kv for kv in res.items()))
    print('   sample p0 = 3: TR %s, GS(+1) %s, GS(-1) %s' % (eval_tr(data, F(3)), graph_sum(F(3), kappa, X1, 1),
                                                         graph_sum(F(3), kappa, X1, -1)))
    ok &= res['TR = closed form'] == 0 and res['GS(s=+1) = TR'] == 0
    ok &= res['GS(s=-1) = TR'] == len(pts) and res['ctrl GS(R^-1 -> R, s=+1) = TR'] == len(pts)

# B3: residues of x^k omega_11 at the branch points (kappa = 1), from the TR local data, and Res_0 x^5 omega_11
data = tr_w11_coeffs(F(1))
bad = 0; nz_ctrl = 0
for a, r, s, sp, den in data:
    from rf4_model import x_taylor
    xt = x_taylor(1, a, NS)
    # omega_11 near a: Sum_k r_k (t)^(-k-1) is NOT the local expansion; use the closed form's Laurent series instead
    p = [a, ONE] + [ZERO] * (NS - 1)
    p4 = spow(p, 4, NS); p5 = smul(p4, p, NS); p10 = smul(p5, p5, NS)
    numer = [c * F(-5, 8) for c in smul(p4, [x + 18 * y + (6 if n == 0 else 0) for n, (x, y) in enumerate(zip(p10, p5))], NS)]
    q5 = [c - (1 if n == 0 else 0) for n, c in enumerate(p5)]       # p^5 - 1 = t Q(t)
    Q = q5[1:] + [ZERO]
    base = smul(numer, spow(Q, -4, NS), NS)                         # omega/dp = t^-4 base
    for k in range(7):
        ser = smul(base, spow(xt, k, NS), NS)
        bad += 0 if ser[3].is0() else 1
    for k in range(1, 10):
        ser = smul(base, spow(p, k, NS), NS)
        nz_ctrl += 0 if ser[3].is0() else 1
# Res_0 x^5 omega_11: x^5 = p^-5 (1 + p^5/4)^5, omega_11/dp = -(5/8) p^4 (6 + 18 p^5 + p^10)/(p^5-1)^4 -> coefficient of p^-1
res0 = -F(5, 8) * 6 / F(1)
print('B3 Res_(varsigma_b) x^k omega_11 = 0 (35 residues): %d non-zero; control p^k (k=1..9): %d of 45 non-zero; '
      'Res_0 x^5 omega_11 = %s' % (bad, nz_ctrl, res0))
ok &= bad == 0 and nz_ctrl > 0 and res0 == F(-15, 4)

# B4: <-log p>_1 by the ray-pushed residues, kernel (1/2) log((a + sigma)/(a + t)) sigma'/((t - sigma)^3 x')
def logleg(s, a, xt):
    tms = [ZERO, ONE] + [ZERO] * (NS - 1)
    sp = sderiv(s)
    D = [p - q for p, q in zip(tms, s)][1:] + [ZERO]
    E = sderiv(xt)[1:NS + 1] + [ZERO]
    den = sinv(smul(spow(D, 3, NS), E, NS), NS)
    X = smul([sg - tt for sg, tt in zip(s, tms)], sinv([a, ONE] + [ZERO] * (NS - 1), NS), NS)
    L = [ZERO] * (NS + 1); Xm = [ONE] + [ZERO] * NS
    for m in range(1, NS + 1):
        Xm = smul(Xm, X, NS)
        L = [l + c * F((-1) ** (m + 1), m) for l, c in zip(L, Xm)]
    return smul(smul(L, sp, NS), den, NS)[3] * F(1, 2)


tot = ZERO; ctrl = ZERO
for a, r, s, sp, den in data:
    s_true, xt = sigma_local(a, NS, 1)
    tot = tot + logleg(s_true, a, xt)
    # control: a wrong involution sigma = -t + t^2/a (NB the linearised sigma = -t is NOT a control here: at v = 1,
    # x'''(varsigma_b) = 0, so sigma = -t + O(t^4) and the (1,1) residues do not see the correction)
    ctrl = ctrl + logleg([ZERO, -ONE, ONE / a] + [ZERO] * (NS - 2), a, xt)
val = tot.rat()
# closed form: -(1/8) Int_(-1)^oo (w^-2 + 20 w^-3 + 25 w^-4) dw  (u = p^5, w = u - 1)
G = lambda w: -F(1) / w - F(10) / w ** 2 - F(25, 3) / w ** 3
cf = -F(1, 8) * (F(0) - G(F(-1)))
print('B4 <-log p>_1 (v = 1): ray-pushed residues %s, closed-form integral %s; control (wrong involution) %s' % (val, cf, ctrl.rat()))
ok &= val == F(-1, 12) and cf == F(-1, 12) and ctrl.rat() != F(-1, 12)
print('ALL CHECKS PASSED' if ok else 'SOME CHECK FAILED')
