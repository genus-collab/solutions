# Polar part test near the conifold. usage: py test_cf.py g field(p|q) N e1 e2 frame [perturb]
import sys, time
from gsum import *

def run(g, fld, N, e1, e2, frame='hodge', perturb=None, percol=True):
    F = FpField() if fld == 'p' else QWField()
    K = needed_K(g)
    t0 = time.time()
    dF = None
    if perturb == 'F+1/9':
        dF = lambda su: L.const(F, F.fr(Fr(1, 9)), N + 2)
    if perturb == 'F+xi/5':
        dF = lambda su: su['xi'] * F.fr(Fr(1, 5))
    D, su = conifold_data(F, e1, e2, N, K, frame, dF)
    per = []
    tot = D.total(g, per)
    tc = su['tc']
    tci = tc.inv()
    pw = L.const(F, F.one, N)
    for _ in range(2 * g - 2):
        pw = pw * tci
    b = beta(2 * g - 2, e1, e2) * Fr(-3) ** (g - 1)
    diff = tot - pw * F.fr(b)
    print('g', g, 'field', F.name, 'eps', e1, e2, 'frame', frame, 'perturb', perturb, 'time %.1f' % (time.time() - t0))
    print(' total val', tot.val, 'prec', tot.prec, ' expected leading', b)
    print(' leading coeff (xi^%d) == expected:' % (2 - 2 * g), F.is0(F.red(tot.coeff(2 - 2 * g) - F.fr(b))))
    pol = [(n, diff.coeff(n)) for n in range(2 - 2 * g, 0)]
    print(' polar part of G_g - (-3)^(g-1) beta/t_c^(2g-2), all zero:', all(F.is0(c) for n, c in pol), ' first bad:', [n for n, c in pol if not F.is0(c)][:3])
    if fld == 'p':
        print(' t_c^0 coeff of difference (rat.recon.):', rr(diff.coeff(0)))
    # per colouring
    if percol and frame == 'hodge':
        bad = []
        worst = 0
        for gs, m, col, t in per:
            nv = len(gs)
            if nv == 1 and sum(c for _, c in m) == 0:
                # vertex graph
                continue
            if t.val < 0:
                bad.append((gs, m, col, t.val))
        print(' per colouring: edge-graph terms with negative valuation:', len(bad), bad[:3], ' #terms', len(per))
        for gs, m, col, t in per:
            if len(gs) == 1 and sum(c for _, c in m) == 0:
                if col == (0,):
                    d0 = t - pw * F.fr(b)
                    print(' vertex col 0: val', t.val, ' minus polar term: polar zero', all(F.is0(d0.coeff(n)) for n in range(2 - 2 * g, 0)))
                else:
                    print(' vertex col', col, 'val', t.val)
    return tot, diff

if __name__ == '__main__':
    g = int(sys.argv[1])
    fld = sys.argv[2]
    N = int(sys.argv[3])
    e1 = Fr(sys.argv[4])
    e2 = Fr(sys.argv[5])
    frame = sys.argv[6] if len(sys.argv) > 6 else 'hodge'
    pert = sys.argv[7] if len(sys.argv) > 7 else None
    run(g, fld, N, e1, e2, frame, pert)
