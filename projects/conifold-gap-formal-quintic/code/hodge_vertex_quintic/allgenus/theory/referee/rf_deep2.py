# -*- coding: utf-8 -*-
"""rf_deep2.py -- ML_Lambda on the singular lists that probe the DEEP part of T_v (frames that depend on the i >= 1 terms of the
vanishing series: frame(6,0), (7,0), (7,4), (8,0), (8,3), (8,4), ...), which th_mla.py and rf_mla.py skip because one flat term of
their expansion has psi-deficit > 6 (rf_mla_scope.log).  Here a flat term with deficit > 6 is admitted only if it is a pure
tau_0(H)^m correlator, whose value is the gap Phi_h^(m)(1) (Phi_h = c_h v^(2-2h), -(1/12) log v, -(5/2) v^2 log v); every other
flat term comes from MVertex6 with the referee's own M_0 (third prime p3).  Lists with another high-deficit flat term are counted
as not evaluable.  Control C4: the same lists with the frames built from the wrong vanishing series (sign of (-1/(4z))^i flipped).
Input beyond the vertex code: the gap F_h = c_h v^(2-2h) (VERIFIED h <= 5 in the project; c_2..c_5 = 1/1200, 1/25200, 1/180000, 1/660000).
Run: py -I -B rf_deep2.py   Log: rf_deep2.log (appended).
"""
import os, sys, time, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from rf_mla_core import Corr, frames_upto, singular_lists, w
from rf_mla import p3
from rf_flat_core import yc_monomials
CG = {2: Fr(1, 1200), 3: Fr(1, 25200), 4: Fr(1, 180000), 5: Fr(1, 660000)}


def gap(h, m):
    """Phi_h^(m)(1)"""
    if h == 0:
        if m < 3: raise ValueError
        out = Fr(-5)
        for j in range(m - 3): out *= -(j + 1)
        return out
    if h == 1:
        if m < 1: raise ValueError
        out = Fr(-1, 12)
        for j in range(1, m): out *= -j
        return out
    out = CG[h]
    for j in range(m): out *= (2 - 2 * h - j)
    return out


def evaluate(C, h, vecs):
    """(value, number of gap terms used) or None if a non-pure-H flat term has deficit > 6"""
    p = C.p; tot = 0; ngap = 0
    for combo in itertools.product(*[list(v.items()) for v in vecs]):
        cf = Fr(1)
        for (_, (x, e)) in combo: cf *= x
        ins = [kc for (kc, _) in combo]
        if 2 * h - 2 + len(ins) <= 0 or C.deficit(h, ins) < 0: continue
        if C.deficit(h, ins) > C.kmax:
            if all(kc == (0, 1) for kc in ins): val = C.R.of(gap(h, len(ins))); ngap += 1
            else: return None
        else: val = C.flat(h, ins)
        tot = (tot + C.R.of(cf) * val) % p
    return tot, ngap


def main():
    LOG = open(os.path.join(HERE, "rf_deep2.log"), "a", encoding="ascii")
    def log(s): print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    t0 = time.time(); C = Corr(p3())
    log("=" * 100); log("rf_deep2 [%s], p3 = %d" % (time.strftime("%Y-%m-%dT%H:%M:%S"), C.p))
    Ybad = [(b, j, n, cf * (-1) ** (n - j)) for (b, j, n, cf) in yc_monomials(60)]
    for h, NT in ((2, 4), (2, 5), (3, 3), (3, 4), (4, 2), (4, 3)):
        th = time.time(); wmax = 2 * h - 2 + NT + 2
        FR = frames_upto(wmax); FB = frames_upto(wmax, Y=Ybad)
        sens = set(c for c in FR if FR[c] != FB[c])
        st = dict(n=0, ev=0, zero=0, fail=[], noev=0, gapused=0, c4=0, ex=None)
        for combo, m, e in singular_lists(h, NT, FR):
            if e >= 0 or not any(c in sens for c in combo): continue
            st['n'] += 1
            H = [{(0, 1): (Fr(1), Fr(0))}] * m
            r = evaluate(C, h, [FR[c] for c in combo] + H)
            if r is None: st['noev'] += 1; continue
            val, ng = r; st['ev'] += 1; st['gapused'] += (ng > 0)
            if val: st['fail'].append((combo, m, e, C.fmt(val)))
            else:
                st['zero'] += 1
                if st['ex'] is None: st['ex'] = (combo, m, e)
            rb = evaluate(C, h, [FB[c] for c in combo] + H)
            if rb is not None and rb[0]: st['c4'] += 1
            if time.time() - th > 1500: log("   time limit reached at genus %d" % h); break
        log("  genus %d, <= %d insertions: %d singular lists with a deep frame; %d evaluable (%d of them use a gap value), %d vanish, "
            "%d FAIL; %d not evaluable; control C4 non-zero on %d; first vanishing list %s  (%.0f s)"
            % (h, NT, st['n'], st['ev'], st['gapused'], st['zero'], len(st['fail']), st['noev'], st['c4'], st['ex'], time.time() - th))
        if st['fail']: log("     FAILURES (first 5): %s" % st['fail'][:5])
    log("  total time %.0f s" % (time.time() - t0)); LOG.close()


if __name__ == "__main__":
    main()
