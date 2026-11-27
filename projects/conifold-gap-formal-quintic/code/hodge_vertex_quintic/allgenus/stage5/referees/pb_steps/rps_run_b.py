# rps_run_b.py -- second half of rps_run.py (T3, T4, summary); executed inside rps_run.py (shares its globals).
# ------------------------------------------------------------------ T3
log('T3  key lemma: Sum_m pi_m s^(-m) I^C_m(v(u), z) = V on admissible keys E\' + (4b-1)/5 <= %d' % NPI)
t = time.time()
Vp2 = rm.vpowers(Vc_, 4 * (2 * NC + 6 - rc.jmin_for(0, E)) // 5 + 8, NC)
def keys_layer(lo, hi):
    """keys (Ep, b): b = 4 mod 5, lo <= Ep + (4b-1)/5 <= hi, Ep <= E (Ep may be <= 0 for b > 0)."""
    out = []
    b = 4
    while (4 * b - 1) // 5 <= hi + 8:          # margin: b beyond the degree bound must give 0 as well
        b += 5
    b -= 5
    while True:
        sh = (4 * b - 1) // 5
        if lo - sh > E:
            break
        for Ep in range(lo - sh, min(hi - sh, E) + 1):
            out.append((Ep, b))
        b -= 5
    return out
adm = keys_layer(0, NPI)
bad = [k for k in adm if rm.resid_orig(pi, Vp2, k[0], k[1], NPI) != V.get(k, 0)]
nz = sum(1 for k in adm if V.get(k, 0))
npos = sum(1 for k in adm if k[1] > 0)
nneg = sum(1 for k in adm if k[0] <= 0)
check('residue side = V on all admissible keys', not bad,
      '(%d keys: %d non-zero, %d with b > 0, %d with E\' <= 0; first bad %s)' % (len(adm), nz, npos, nneg, bad[:1]))
lay = keys_layer(NPI + 1, NPI + 1)
bad1 = [k for k in lay if rm.resid_orig(pi, Vp2, k[0], k[1], NPI) != V.get(k, 0)]
bad2 = [k for k in lay if rm.resid_orig(pi, Vp2, k[0], k[1], NPI + 1) != V.get(k, 0)]
control('truncation: pi only through u^%d on the layer E\' + (4b-1)/5 = %d' % (NPI, NPI + 1), len(bad1) > 0,
        '(%d of %d keys wrong; e.g. %s)' % (len(bad1), len(lay), bad1[:2]))
check('the same layer with pi through u^%d' % (NPI + 1), not bad2, '(%d keys)' % len(lay))
top = 0
for (Ep, b) in adm:
    sh = (4 * b - 1) // 5
    i = Ep + sh
    one = {m: {k: c for k, c in d.items() if k[0] == i} for m, d in pi.items()}
    if rm.resid_orig(one, Vp2, Ep, b, NPI) != 0:
        top += 1
check('the bound is attained: pi at u^(E\'+(4b-1)/5) contributes on some keys', top > 0, '(%d keys)' % top)
for name, pic in (('sign of pi_-1 flipped', {m: ({k: -c for k, c in d.items()} if m == -1 else d) for m, d in pi.items()}),
                  ('pi_2 and pi_3 exchanged', {**pi, 2: pi[3], 3: pi[2]}),
                  ('pi from the reduction with V -> 1', rm.reduce_forms(beta, [Fr(1)] + [Fr(0)] * NC, NC))):
    badc = [k for k in adm if rm.resid_orig(pic, Vp2, k[0], k[1], NPI) != V.get(k, 0)]
    control(name, len(badc) > 0, '(%d of %d keys wrong)' % (len(badc), len(adm)))
log('    (%.1f s)' % (time.time() - t))
# ------------------------------------------------------------------ T4
log('T4  (S2) and Prop. S (model rows, colour 0) through t^%d' % M)
t = time.time()
F, piv, resid = rs.solve_F(M)
check('Prop. S: the pivot of F_0[n] in the closing relation at t^(n+1) is 4n (n = 1..%d)' % M,
      all(piv[n] == 4 * n for n in piv), '(pivot(1), pivot(%d) = %s, %s)' % (M, piv[1], piv[M]))
check('Prop. S: closing relation satisfied through t^%d' % M, all(r == 0 for r in resid))
check('F_0 = 1 - 3/20 t + 9/800 t^2 + 11571/80000 t^3 - 2415609/6400000 t^4 (UR s.2)',
      F[0][:5] == [1, Fr(-3, 20), Fr(9, 800), Fr(11571, 80000), Fr(-2415609, 6400000)])
check('F_b(0) = eps_b = (+,-,-,+,+)', [F[b][0] for b in range(5)] == list(rs.EPS))
SP = [rs.sp_model(rs.C[b], M) for b in range(5)]
check('(S2): eps_b SP^norm[eta^(c_b) e^(-(1/eta+eta^4/4)/t) d eta] = F_b(t), all b', all(
      [rs.EPS[b] * c for c in SP[b]] == F[b] for b in range(5)))
G0 = [rs.EPS[0] * c for c in SP[0]] + [Fr(0)]
rows, cl = rs.chain(G0, M)
check('the SP rows satisfy the model system (T-relations) and its closing', all(
      rows[b][:M] == [rs.EPS[b] * c for c in SP[b][:M]] for b in range(5)) and cl[:M] == G0[:M])
control('all signs + in (S2)', any(SP[b] != F[b] for b in range(5)))
control('weights c_b -> c_b + 1 in (S2)', any([rs.EPS[b] * c for c in rs.sp_model(rs.C[b] + 1, M)] != F[b] for b in range(5)))
rows_s, cl_s = rs.chain(G0, M, shift=1)
control('system with T_k -> T_(k+1) not satisfied by the SP rows', cl_s[:M] != G0[:M] or any(
      rows_s[b][:M] != [rs.EPS[b] * c for c in SP[b][:M]] for b in range(1, 5)))
log('    (%.1f s)' % (time.time() - t))
exec(open(os.path.join(HERE, 'rps_run_c.py')).read())
# ------------------------------------------------------------------ summary
log('total %.1f s' % (time.time() - t0))
log('ALL CHECKS PASSED' if not fails else 'FAILURES: %s' % fails)
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
with open(os.path.join(HERE, LOGREL), 'w', encoding='utf-8') as fh:
    fh.write('\n'.join(lines) + '\n')
print('log written to ' + LOGREL)
