# -*- coding: utf-8 -*-
"""cr4_run.py -- Stage-4 critic: do the three reductions (evid, poly, gamma) coincide?  Exact over Q, own code only.

Usage:  py -I -B cr4_run.py [Q] [MCAN] [NSP]     (defaults 8 12 40)
  A  set-up checks (periods, model rows, QRR constants)
  B  QRR WKB polynomials: (H), resonances, degrees, ph_n(0) = model, free coefficients
  C  the three readings of [z^5q] c(z): KA (evid referee A, period free), KE (evid Lemma N + Lemma P, row 1 with C0, C1),
     KP (poly s.3, first row with C0): equal on Phi^QRR (= 0) and on a perturbed Phi (= referee-A values), q <= Q, two frames
  D  normalisations: f^can (KA = 0, no QRR input) = f^QRR, connection form ph^can_n(1) = h_n, Ginf closed form
  E  full-Stirling stationary phase at fixed u0 (poly/gamma integrand): Phi^SP = Phi^QRR = Phi^can, flat value c_5(2), controls
Log: logs/cr4_run_<Q>_<MCAN>_<NSP>.log (relative to this folder).
"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import comb
from cr4_core import Field, hcoef, cyclo_powersums, bernoulli, tau_ser, periods, model_rows, Wkb, KA, Gauge, taupows, peval, sinv, smul
from cr4_sp import sp_phi

Q, MCAN, NSP = (int(a) for a in (sys.argv[1:] + ['8', '12', '40'][len(sys.argv) - 1:])[:3])
LOGREL = 'logs/cr4_run_%d_%d_%d.log' % (Q, MCAN, NSP)
LOG = open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8')
T0 = time.time(); FAILS = []

def out(s):
    print(s); LOG.write(s + '\n'); LOG.flush()

def check(name, ok, info=''):
    out(('PASS ' if ok else 'FAIL ') + name + (('   ' + str(info)) if info != '' else ''))
    if not ok: FAILS.append(name)

F = Field()
NMAX = 5 * max(Q, MCAN)
out('cr4_run: Q = %d (z^%d), MCAN = %d (z^%d), NSP = %d; field Q (exact); log %s' % (Q, 5 * Q, MCAN, 5 * MCAN, NSP, LOGREL))

# ---------------------------------------------------------------- A: set-up
per0 = periods(F, 4 * Q + 1)
check('A1 gamma = 7/10, t_c = u + (7/10)u^2 + (41/75)u^3 + (9/20)u^4, f_0 = 1 + (2/625)u^3 + .. (frame (0,0))',
      per0['gamma'] == Fr(7, 10) and per0['tc'][:5] == [0, 1, Fr(7, 10), Fr(41, 75), Fr(9, 20)] and per0['f0'][:4] == [1, 0, 0, Fr(2, 625)])
FR = [(0, 0), (Fr(1, 3), Fr(-2, 7)), (Fr(7, 11), Fr(-3, 5))]
check('A1 C0(0) C1(0) = -1 in the frames (0,0), (1/3,-2/7), (7/11,-3/5)', all(periods(F, 3, a, b)['C0'][0] * periods(F, 3, a, b)['C1'][0] == -1 for a, b in FR))
Fb = model_rows(F, NMAX + 1)
UR_F0 = [1, Fr(-3, 20), Fr(9, 800), Fr(11571, 80000), Fr(-2415609, 6400000), Fr(47735919, 128000000), Fr(52989974037, 25600000000), Fr(-7749624763881, 512000000000)]
check('A2 model F_0 through t^7 = UR s.2 values; F_b(0) = (1,-1,-1,1,1); closure F_0 = -T_0 F_4 (asserted)', Fb[0][:8] == UR_F0 and [f[0] for f in Fb] == [1, -1, -1, 1, 1])
h = hcoef(F, NMAX + 1)
check('A3 h_1 = -3/20, h_5 = -257071/8064000000 (QRR, weight -5)', h[1] == Fr(-3, 20) and h[5] == Fr(-257071, 8064000000))
p = cyclo_powersums(400)
z5 = [m for m in range(1, 400, 2) if p[m] == 0]
check('A3 p_m = Sum_j (1-zeta^j)^-m vanishes for odd m < 400 exactly at the odd multiples of 5 (evid s.8.2)', z5 == list(range(5, 400, 10)), 'zeros %d' % len(z5))

# ---------------------------------------------------------------- B: QRR polynomials
t = time.time(); W = Wkb(F, NMAX, h, 'qrr')
nH = len([n for n in range(1, NMAX + 1) if n % 5])
check('B1 (H): ph_n(1) = h_n for all %d n <= %d with 5 not | n (no free parameter)' % (nH, NMAX), W.H_fail == [], W.H_fail[:5])
check('B2 resonances unobstructed; deg ph_n <= floor(4n/5) for n <= %d' % NMAX, W.res_ok and W.deg_ok)
check('B2 deg ph_n = floor(4n/5) exactly for 3 <= n <= %d (n = 2: deg 0 < 1)' % NMAX, all(W.deg_eq[2:]) and not W.deg_eq[1])
check('B3 ph_n(0) = F_(0,n) (model) for n <= %d: automatic, not a criterion' % NMAX, all(W.ph[n][0] == Fb[0][n] for n in range(NMAX + 1)))
FV = {1: Fr(4537, 3937500), 2: Fr(-177640111151, 31007812500000), 3: Fr(-284328315069009259339, 6226756347656250000000)}
check('B4 f_1, f_2, f_3 (QRR) = the evid / poly / gamma values', all(W.f[m] == FV[m] for m in FV), W.f[1])
out('     B done in %.1f s' % (time.time() - t))

# ---------------------------------------------------------------- C: the three readings of [z^5q] c
t = time.time(); K = 4 * Q + 1; tau = tau_ser(F, K); tpm = taupows(F, tau, 5 * Q, K, -1)
REFA = {1: Fr(1, 3), 2: Fr(6048257071, 24192000000), 3: Fr(780415789166012126329, 3901685760000000000000)}
Wq = Wkb(F, 5 * Q, h, 'qrr'); sh = {m: Fr(1, m + 2) for m in range(1, Q + 1)}
Wa = Wkb(F, 5 * Q, h, 'qrr', fv=sh)                                            # referee-A rule: ph_5q(1) = h_5q + 1/(q+2)
Wb = Wkb(F, 5 * Q, h, 'list', fv={m: W.f[m] + sh[m] for m in sh})            # f_q = f_q^QRR + 1/(q+2) (fixed numbers)
for (a, b) in FR[:2]:
    per = periods(F, K, a, b)
    for lab, Ww in (('QRR', Wq), ('pert. A (ph_5q(1) = h_5q + 1/(q+2))', Wa), ('pert. B (f_q^QRR + 1/(q+2))', Wb)):
        G = Gauge(F, Ww, Q, per, tau, Fb)
        ka = [KA(F, Ww.P[1], q, tpm, Fb[2]) for q in range(1, Q + 1)]
        ke = [G.KE(q) for q in range(1, Q + 1)]; kp = [G.KP(q) for q in range(1, Q + 1)]
        same = all(ka[i][0] == ke[i][0] == kp[i] for i in range(Q))
        reg = all(x[1] for x in ka) and all(x[1] for x in ke) and G.div_ok
        cz = smul(F, [sum(Ww.ph[n]) for n in range(5 * Q + 1)], sinv(F, h[:5 * Q + 1], 5 * Q + 1), 5 * Q + 1)   # Phi(z,1)/h(z)
        pred = all(ka[q - 1][0] == cz[5 * q] for q in range(1, Q + 1))
        if lab == 'QRR': val, tag = all(x[0] == 0 for x in ka), 'all 0 (D = 1 to z^%d)' % (5 * Q)
        elif Ww is Wa: val, tag = all(ka[q - 1][0] == REFA[q] for q in REFA) and all(x[0] != 0 for x in ka), '= evid referee A (rf_pfree) for q <= 3'
        else: val, tag = ka[1][0] == Fr(1, 4) - W.f[1] / 3 and all(x[0] != 0 for x in ka), 'K_2 = 1/4 - f_1/3 (predicted), all non-zero'
        check('C1 frame (%s,%s), %s: KA = KE = KP for q <= %d; REG (KA, KE) and first-row divisibility' % (a, b, lab, Q), same and reg)
        check('C1 frame (%s,%s), %s: values %s; K_q = [z^5q] Phi(z,1)/h(z)' % (a, b, lab, tag), val and pred,
              '' if lab == 'QRR' else [str(x[0]) for x in ka[:2]])
G0 = Gauge(F, Wq, Q, per0, tau, Fb); m5 = G0.g[4][4]
check('C2 m5 = [u^4] g_4 = p\'_(4,1)(0) = -43397/807975000 from the first row (no free coefficient enters)', m5 == Fr(-43397, 807975000))
Ws = Wkb(F, 5 * Q, h, 'list', fv={m: W.f[m] + (m5 if m == 1 else 0) for m in range(1, Q + 1)}); Gs = Gauge(F, Ws, Q, per0, tau, Fb)
check('C2 stored normalisation (f_1 + m5): KA_1 = KE_1 = KP_1 = m5; regularity blind', KA(F, Ws.P[1], 1, tpm, Fb[2]) == (m5, True) and Gs.KE(1) == (m5, True) and Gs.KP(1) == m5 and Gs.div_ok)
t1 = taupows(F, [1] + [0] * (K - 1), 5 * Q, K, -1)
kt = [KA(F, Wq.P[1], q, t1, Fb[2]) for q in (1, 2)]
check('C3 control tau -> 1 (wrong natural variable) breaks KA / REG', any(x[0] != 0 or not x[1] for x in kt), kt[0])
pb = periods(F, K, 0, 0, dg=Fr(1, 7)); Gb = Gauge(F, Wq, Q, pb, tau, Fb)
out('     C4 probe (not a control): vanishing period g + f_2/7 (not a conifold frame): KE = %s, KP = %s (KE sees C1 only via C1(0); KP only C0)'
    % ([str(Gb.KE(q)[0]) for q in (1, 2)], [str(Gb.KP(q)) for q in (1, 2)]))
out('     C done in %.1f s' % (time.time() - t))

# ---------------------------------------------------------------- D: normalisations
t = time.time(); Kc = 4 * MCAN + 1; tauc = tau_ser(F, Kc); tpc = taupows(F, tauc, 5 * MCAN, Kc, -1)
Wc = Wkb(F, 5 * MCAN, h, 'can', ka=(tpc, Fb[2]))
check('D1 f^can (fixed by KA = 0, no QRR input) = f^QRR for m <= %d (z^%d)' % (MCAN, 5 * MCAN), all(Wc.f[m] == W.f[m] for m in range(1, MCAN + 1)))
check('D1 connection form: ph^can_n(1) = h_n for ALL n <= %d (incl. 5 | n): D_0 = 1 mod z^%d, i.e. (F_g) for g <= %d' % (5 * MCAN, 5 * MCAN + 1, (5 * MCAN + 2) // 3),
      all(sum(Wc.ph[n]) == h[n] for n in range(5 * MCAN + 1)))
check('D1 REG of the canonical function (KA reading) for q <= %d' % MCAN, all(KA(F, Wc.P[1], q, tpc, Fb[2]) == (0, True) for q in range(1, MCAN + 1)))
B = bernoulli(5 * MCAN + 2)
bp = lambda n, x: sum(comb(n, k) * B[k] * x ** (n - k) for k in range(n + 1))
lg = [0] + [(-1) ** (l + 1) * bp(5 * l + 1, Fr(1, 5)) / (l * (5 * l + 1)) for l in range(1, MCAN + 1)]
gi = [Fr(1)]
for m in range(1, MCAN + 1): gi.append(sum(l * lg[l] * gi[m - l] for l in range(1, m + 1)) / m)
check('D2 Ginf closed form (gamma s.0 item 6): f^QRR_m = [z^5m] exp(Sum_l (-1)^(l+1) B_(5l+1)(1/5) z^5l/(l(5l+1))), m <= %d' % MCAN,
      all(gi[m] == W.f[m] for m in range(1, MCAN + 1)), 'f_1 = B_6(1/5)/6: %s' % (bp(6, Fr(1, 5)) / 6 == W.f[1]))
out('     D done in %.1f s' % (time.time() - t))

# ---------------------------------------------------------------- E: stationary phase of the full-Stirling MB integrand
t = time.time(); Wz = Wkb(F, NSP, h, 'zero')
check('E1 u0 = 0: Phi^SP(w; 0) = model F_0(w) through w^%d (the model is the conifold limit of the integrand)' % NSP, sp_phi(F, 0, NSP) == Fb[0][:NSP + 1])
for u0 in (Fr(1, 2), Fr(-1, 3), Fr(3, 7)):
    s = sp_phi(F, u0, NSP); qv = [peval(F, W.ph[n], u0) for n in range(NSP + 1)]
    z0 = [peval(F, Wz.ph[n], u0) for n in range(NSP + 1)]
    rat = smul(F, s, sinv(F, z0, NSP + 1), NSP + 1)               # = F_inf^SP(z), z^5m = u0^4m w^5m
    pure = all(rat[n] == 0 for n in range(1, NSP + 1) if n % 5)
    fsp = {m: rat[5 * m] / u0 ** (4 * m) for m in range(1, NSP // 5 + 1)}
    check('E2 u0 = %s: Phi^SP = Phi^QRR (= Phi^can) through w^%d, exact' % (u0, NSP), s == qv)
    check('E2 u0 = %s: Phi^SP / Phi^(0) has only w^5m terms ((H) for SP) and f^SP_m = f^QRR_m = f^can_m, m <= %d' % (u0, NSP // 5),
          pure and all(fsp[m] == W.f[m] == Wc.f.get(m, W.f[m]) for m in fsp))
c5 = 2 ** 20 * sp_phi(F, Fr(1, 32), 5)[5]
check('E3 flat dictionary c_n(L) = L^4n ph_n(L^-5): c_5(2) = 2^20 ph^SP_5(1/32) = 11495493197/31500 (gamma logs)', c5 == Fr(11495493197, 31500))
sb = sp_phi(F, Fr(1, 2), 8, b2shift=1); qb = [peval(F, W.ph[n], Fr(1, 2)) for n in range(9)]
check('E4 control b_2 -> b_2 + 1: SP != QRR first at w^3', [n for n in range(9) if sb[n] != qb[n]][:1] == [3])
sd = sp_phi(F, Fr(1, 2), 8, drop5=True)
check('E4 control drop the Gamma(1+5s) Stirling term: SP != QRR first at w^1', [n for n in range(9) if sd[n] != qb[n]][:1] == [1])
Wx = Wkb(F, 6, hcoef(F, 6, w5=5), 'qrr')
check('E4 control QRR twist weight +5 instead of -5: (H) fails first at n = 1', Wx.H_fail[:1] == [1], Wx.H_fail)
out('     E done in %.1f s' % (time.time() - t))

out('CR4_RUN %d %d %d: %s   (%.0f s)' % (Q, MCAN, NSP, 'ALL CHECKS PASSED' if not FAILS else 'FAILED: %s' % FAILS, time.time() - T0))
LOG.close()
