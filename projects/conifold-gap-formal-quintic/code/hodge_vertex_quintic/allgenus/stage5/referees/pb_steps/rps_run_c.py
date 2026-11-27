# rps_run_c.py -- third part of rps_run.py (T5); executed inside rps_run.py (shares its globals).
# T5: the identity of PB s.2.6 that needs no gauge:  Phi^SP = Sum_b eps_b pi_(c_b) tau^((1-c_b)/4) F_b(tau w),
#     with Phi^SP the hatted SP at X = 1 over Q[[u]] (rps_sp.hat_SP), pi_m from T2, F_b from T4 (Prop. S).
import rps_sp as rsp
log('T5  s.2.6 without the gauge: Phi^SP = Sum_b eps_b pi_(c_b) tau^((1-c_b)/4) F_b(tau w), exact')
t = time.time()
NS = min(NC, M)
nS = min(NS, M)
Bh = rc.hat_B(NS)
PhiSP = rsp.hat_SP(Bh, nS, NS)
check('Phi^SP: ph_0 = 1', PhiSP[0] == [1] + [0] * NS)
degs = []
for n in range(nS + 1):
    nzs = [i for i in range(NS + 1) if PhiSP[n][i]]
    degs.append(max(nzs) if nzs else -1)
okdeg = all(degs[n] <= (4 * n) // 5 for n in range(nS + 1) if (4 * n) // 5 < NS)
check('ph^SP_n has no u-power above floor(4n/5) (within u^%d)' % NS, okdeg, '(degrees %s)' % degs)
check('ph^SP_1 = -3/20 (= h_1) and ph^SP_2 = 9/800', PhiSP[1] == [Fr(-3, 20)] + [0] * NS and PhiSP[2] == [Fr(9, 800)] + [0] * NS)
F5, _, _ = rs.solve_F(nS)
R = rsp.rhs_s26(pi, F5, rs.EPS, rs.C, nS, NS)
first = None
for n in range(nS + 1):
    for i in range(NS + 1):
        if R[n][i] != PhiSP[n][i] and first is None:
            first = (n, i)
check('Phi^SP = Sum_b eps_b pi_(c_b) tau^((1-c_b)/4) F_b(tau w), w^0..w^%d, u^0..u^%d' % (nS, NS), first is None,
      '(first mismatch (w, u) = %s)' % (first,))
for name, kw in (('without the signs eps_b', dict(use_eps=False)), ('without the weights tau^((1-c_b)/4)', dict(use_tau=False)),
                 ('F_b(w) in place of F_b(tau w)', dict(tau_in_F=False))):
    Rc = rsp.rhs_s26(pi, F5, rs.EPS, rs.C, nS, NS, **kw)
    fc = next(((n, i) for n in range(nS + 1) for i in range(NS + 1) if Rc[n][i] != PhiSP[n][i]), None)
    control(name, fc is not None, '(first mismatch (w, u) = %s)' % (fc,))
log('    (%.1f s)' % (time.time() - t))
