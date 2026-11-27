# -*- coding: utf-8 -*-
"""ll_eo.py -- differential-level checks of the steps of LOGLEG_PROOF.md on x = W = p^4/4 + u/p, y = p (u = 1), e-basis tensors.
   py -I -B ll_eo.py <out.log> <field: 89 | 61 | exact> <hmax-total>   (all (h,n) with 2h-2+n >= 1 and 2h-2+(n+1) <= hmax-total)
(E0) genus 0 ("Rauch"):  d_v|_x B = -(Lie_xi1 + Lie_xi2) B = Sum_l r_l e_{l,1} (x) e_{l,1},  r_l = 1/(5 c_l^3),  equals Int_0^oo omega_{0,3}(., ., p').
(E1) n = 0:  d_v F_h = (2-2h) F_h(1) (homogeneity) equals Int_0^oo omega_{h,1} = <omega_{h,1}, -log p>.
(E2) EO variational formula at the level of differentials:  d_v|_p omega_{h,n} - Sum_i Lie_xi_i omega_{h,n} = Int_0^oo omega_{h,n+1}(., p').
(E3) boundary terms:  Int_0^oo Lie_xi(p') omega_{h,n+1}(., p') = 0 (Lie only in the integrated slot).
Controls (must FAIL): (E0) with the opposite sign; (E2) without the Lie terms, and with Lie_xi of the opposite sign.
"""
import os, sys, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
from ll_core import QzField, make_fp, P89, P61, setup, Legs, contract_last, dv_p, Lie, sub, fpow, rat
from ll_falsifier import Log


def lie_last(lie, W):
    out = {}
    for key, C in W.items():
        for e, a in lie.table(*key[-1]).items():
            kk = key[:-1] + (e,); v = C * a
            out[kk] = out[kk] + v if kk in out else v
    return {k: v for k, v in out.items() if not v.iszero()}


def main(out, fld, tot):
    log = Log(out); log('ll_eo.py %s %s  (allgenus/stage2/logleg, %s)' % (fld, tot, time.strftime('%Y-%m-%d %H:%M')))
    F = {'89': make_fp(P89), '61': make_fp(P61), 'exact': QzField()}[fld]
    comp = set((h, m) for h in range(5) for m in range(1, 7) if 2 * h - 2 + m <= min(tot, 6)) | (set([(3, 3), (4, 1), (4, 2)]) if tot >= 8 else set())
    cases = [(h, n) for h in range(0, 5) for n in range(0, 6) if 2 * h - 2 + n >= 1 and (h, n + 1) in comp]
    kmx = max(6 * h - 5 + 2 * (n + 1) for h, n in cases); cur, tr = setup(F, 2 * kmx + 8); c = cur.c
    L = Legs(F, c); lie = Lie(F, c); lie_m = Lie(F, c, sign=-1)
    # (E0)
    I3 = contract_last(tr.omega(0, 3), [L.logleg])
    R = {((l, 1), (l, 1)): F.one / (fpow(F, c[l], 3) * 5) for l in range(5)}
    log.check(not sub(F, I3, R), '(E0) Int_0^oo omega_{0,3}(p1,p2,.) = -(Lie_xi1 + Lie_xi2) B = Sum_l e_{l,1}(p1) e_{l,1}(p2)/(5 c_l^3)  [%d terms]' % len(I3))
    log.check(bool(sub(F, I3, {k: -v for k, v in R.items()})), '(E0) control: the opposite sign FAILS (as it must)')
    # (E1)
    for h in range(1, tot // 2 + 2):
        if 2 * h - 1 > tot: continue
        val = contract_last(tr.omega(h, 1), [L.logleg]).get((), F.zero)
        if h == 1:
            log.check((val - F.of(Fr(-1, 12))).iszero(), '(E1) h = 1: Int_0^oo omega_{1,1} = %s = d_v[-(1/12) log v] at v = 1' % rat(F, val))
        else:
            Fh = tr.Fg(h)
            log.check((val - Fh * (2 - 2 * h)).iszero(), '(E1) h = %d: Int_0^oo omega_{%d,1} = %s = (2-2h) F_h(1) = d_v F_h, F_h(1) = %s' % (
                h, h, rat(F, val), rat(F, Fh)))
    # (E2), (E3)
    for h, n in cases:
        if n == 0: continue
        t0 = time.time(); W = tr.omega(h, n); W1 = tr.omega(h, n + 1)
        D = dv_p(F, c, W, h, n); Lw = lie.apply(W); rhs = contract_last(W1, [L.logleg])
        res = sub(F, sub(F, D, Lw), rhs)
        c1 = bool(sub(F, D, rhs)); c2 = bool(sub(F, sub(F, D, lie_m.apply(W)), rhs))
        log.check(not res and c1 and c2, '(E2) (%d,%d): d_v|_p omega - Sum Lie_xi omega = Int_0^oo omega_{%d,%d}(., p\'): %d terms, residual %d; '
                  'controls fail: no-Lie %s, wrong-sign-Lie %s (%.1f s)' % (h, n, h, n + 1, len(rhs), len(res), c1, c2, time.time() - t0))
        bt = contract_last(lie_last(lie, W1), [L.logleg])
        log.check(not bt, '(E3) (%d,%d): Int_0^oo Lie_xi(p\') omega_{%d,%d}(., p\') = 0 (boundary terms vanish)' % (h, n, h, n + 1))
    log.close()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
