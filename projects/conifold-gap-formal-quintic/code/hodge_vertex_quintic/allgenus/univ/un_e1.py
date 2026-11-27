# -*- coding: utf-8 -*-
"""un_e1.py -- E1 (regression of the new engine un_engine.Var against numbers on record).

(a) natural rows: the edge/leaf tables rebuilt from the row lists equal those of fq4_split.Split; the run over all inner graphs
    reproduces the stored genus-3 table g3_nat_<p>_nat.json (A3, 42 graphs, u^0..u^3) and the genus-2 values of the notes; all edge
    graphs regular; vertex graph = c_g t_c^(2-2g).
(b) first-order responses (dual numbers, L = 2) to a LEFT dressing of R' by a single constant hatted term
        rows -> (1 + eps zhat^k X_flat) rows,   X_flat = Psi(0) X Psi(0)^-1,  X_ij = Sum Xhat(m, m') zeta^(-(m i + m' j)),
    at [eps^1 u^0] on the one-edge graphs of genus 3 (Mathematica indices 2: loop on g = 2, 5: bridge g1 - g2) and genus 2.
    On record (NATURAL_FORM.md s.1 item 6, nat_kregular_linear.log, ur_euler_linear logs) are the responses to a left dressing of M;
    at leading order in u a left dressing of M by X is a left dressing of R' by -X_flat, so the numbers on record are MINUS ours.
Run: py -I -B un_e1.py [89|61|both]
"""
import os, sys, json, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, PROJ); sys.path.insert(0, HERE)
from un_core import Arith, PRIMES, lift_rows, left_dress
from un_engine import Var, get_base, get_nat, run, collect
from fq3_graphs import outer_graphs
from fq3_rec import Rec
from fq_hat import tc_series
from fq3_msum import tc_expansion3

LOG = open(os.path.join(HERE, "un_e1.log"), "a", encoding="utf-8"); FAILED = []


def log(*x):
    s = " ".join(str(y) for y in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)


def flat_of_fourier(nr, comp):
    """Psi(0) X Psi(0)^-1 as {(a, b): residue} for X_ij = Sum_(m,m') c zeta^(-(m i + m' j))"""
    R = nr['R']; p = nr['p']; Rh = nr['Rh']; z = R.zeta
    X = [[sum(c * pow(z, (-m * i - mp * j) % 5, p) for (m, mp), c in comp.items()) % p for j in range(5)] for i in range(5)]
    Pinv = [[sum(Rh[c][j][0] * e for (c, b2), e in nr['eta_inv'].items() if b2 == b) % p for b in range(5)] for j in range(5)]
    out = {}
    for a in range(5):
        for b in range(5):
            v = sum(Rh[a][i][0] * X[i][j] * Pinv[j][b] for i in range(5) for j in range(5)) % p
            if v: out[(a, b)] = v
    return out


CASES = [   # name, [(k, Fourier components)], expected -(record) on genus-3 graphs (2, 5) or None, expected visibility
    ("zhat H(x)H", [(1, {(3, 3): 1})], (Fr(1, 400), Fr(1, 7200))),
    ("zhat^2 (H(x)1 - 1(x)H)", [(2, {(3, 4): 1, (4, 3): -1})], (Fr(-1, 300), Fr(-7, 14400))),
    ("zhat^3 1(x)1", [(3, {(4, 4): 1})], (Fr(-1, 1200), Fr(-1, 2880))),
    ("Pi-type zhat (H^4(x)H^3 + H^3(x)H^4)", [(1, {(0, 1): 1, (1, 0): 1})], (0, 0)),
    ("Pi-type zhat (H^2(x)1 + 1(x)H^2)", [(1, {(2, 4): 1, (4, 2): 1})], (0, 0)),
    ("zhat^2 (H^4(x)H^2 - H^2(x)H^4)", [(2, {(0, 2): 1, (2, 0): -1})], (0, 0)),
    ("zhat^3 H^4(x)H", [(3, {(0, 3): 1, (3, 0): 1})], (Fr(1, 1200), Fr(1, 14400))),
    ("zhat^3 H^3(x)H^2", [(3, {(1, 2): 1, (2, 1): 1})], (0, 0)),
    ("Euler-dilaton zhat (H - zhat)(x)(H + zhat)", [(1, {(3, 3): 1}), (2, {(3, 4): 1, (4, 3): -1}), (3, {(4, 4): -1})], (0, 0)),
    ("CONTROL zhat (H + zhat)(x)(H - zhat)", [(1, {(3, 3): 1}), (2, {(3, 4): -1, (4, 3): 1}), (3, {(4, 4): -1})], (Fr(1, 150), Fr(7, 7200))),
]


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    for prime in (("89", "61") if which == "both" else (which,)):
        log("=" * 110); log("un_e1: p = 2^%s-1" % prime); ab = ("0", "0")
        for g, K in ((2, 2), (3, 4)):
            t0 = time.time(); base = get_base(prime, g, ab, K); nr = get_nat(prime, g, ab, K); R = nr['R']; p = nr['p']; rec = Rec(R)
            ar = Arith(p, K, 1); v = Var(base, lift_rows(ar, nr['Qh']), 1)
            same_E = all(v.ER[key][i][j] == base.ER[key][i][j] for key in base.ER for i in range(5) for j in range(5))
            same_T = all(v.TR[d][i] == base.TR[d][i] for d in base.TR for i in range(5))
            same_x = all(v.xpow[k] == base.xpow[k] or (v.xpow[k] == [x % p for x in base.xpow[k]]) for k in range(len(base.xpow)))
            check("genus %d: E_R, T_R, xhat rebuilt from the row lists equal the tables of fq4_split.Split" % g, same_E and same_T and same_x, (same_E, same_T, same_x))
            _, tab = run(None, v)
            npol = 2 * g - 2; og = {G['idx']: v.outer_pos[G['key']] for G in outer_graphs(g)}
            bad = {i: c[:npol] for i, c in tab.items() if i != 0 and any(c[:npol])}
            check("genus %d natural run: every outer graph with edges is O(u^%d)" % (g, npol), not bad, bad)
            tc = tc_series(base.C, R, K + 4); V = tc_expansion3(tab[0][:K], tc, R, g)[:npol]
            cg = {2: Fr(1, 1200), 3: Fr(1, 25200)}[g]
            check("genus %d natural run: vertex graph = (%s) t_c^-%d + 0 + .." % (g, cg, npol), V[0] == R.of(cg) and not any(V[1:]), [rec.fmt(x) for x in V])
            if g == 3:
                ref = json.load(open(os.path.join(PROJ, "g3_nat_%s_nat.json" % prime)))['A3']
                mism = [i for i in og if [x % p for x in ref[str(i)][:4]] != tab.get(og[i], [0] * 4)[:4]]
                check("genus 3: all 42 outer graphs equal the stored table g3_nat_%s_nat.json (A3) through u^3" % prime, not mism, mism)
            else:
                check("genus 2: vertex graph [u^0], [u^1] = 1/1200, -7/6000 (notes)", tab[0][:2] == [R.of(Fr(1, 1200)), R.of(Fr(-7, 6000))])
            log("   (%.0f s)" % (time.time() - t0))
        # ---------------- (b) single-term first-order responses
        resp = {}
        for g, K in ((2, 2), (3, 4)):
            base = get_base(prime, g, ab, K); nr = get_nat(prime, g, ab, K); R = nr['R']; p = nr['p']; rec = Rec(R)
            ar = Arith(p, K, 2); Q0 = lift_rows(ar, nr['Qh']); og = {G['idx']: base.outer_pos[G['key']] for G in outer_graphs(g)}
            for name, terms, exp in CASES:
                X = {}
                for k, comp in terms:
                    for (a, b), c in flat_of_fourier(nr, comp).items(): X[(a, b, k)] = ar.mono(c, 0, 1)
                grading = all((a - b - k) % 5 == 0 for (a, b, k) in X)
                v = Var(base, left_dress(ar, Q0, X, base.KZ), 2); _, tab = run(None, v)
                resp[(g, name)] = {i: ar.part(tab.get(og[i], None), 1) for i in og}
                resp[(g, name, 'flat')] = {key: rec.fmt(ar.part(s, 1)[0]) for key, s in X.items()}
                resp[(g, name, 'grading')] = grading
        log("-" * 110); log("first-order responses [eps^1] to rows -> (1 + eps X_flat(zhat)) rows; a response series is listed as [u^0, u^1, ..]")
        for name, terms, exp in CASES:
            r3 = resp[(3, name)]; r2 = resp[(2, name)]
            log("   %s" % name)
            log("      flat entries (a, b, k): %s   Z/5-graded: %s" % (resp[(3, name, 'flat')], resp[(3, name, 'grading')]))
            log("      genus 3: graph 2 (loop on g=2): %s;  graph 5 (bridge g1-g2): %s" % ([rec.fmt(x) for x in r3[2]], [rec.fmt(x) for x in r3[5]]))
            log("      genus 3, all one-edge graphs [u^0]: %s" % {i: rec.fmt(r3[i][0]) for i in sorted(r3) if base_E(3, i) == 1})
            log("      genus 2, graphs 1..7 [u^0]: %s" % {i: rec.fmt(r2[i][0]) for i in sorted(r2)})
            got = (r3[2][0], r3[5][0]); want = tuple(R.of(Fr(x)) for x in exp)
            check("%s: [eps^1 u^0] on genus-3 graphs (2, 5) = %s = MINUS the response on record" % (name, [str(x) for x in exp]), got == want, [rec.fmt(x) for x in got])
        r2 = resp[(2, CASES[8][0])]
        log("   Euler-dilaton at genus 2 (record, M-dressing: vertex -7/128, loop on g1 1/16, bridge g1-g1 -1/128): ours %s" % {i: rec.fmt(r2[i][0]) for i in sorted(r2)})
        # NOTE: the record value -7/128 on the VERTEX graph belongs to the M-dressing (there the total is conserved and the vertex theory
        # changes).  A left dressing of R' with rows a = 2, 3 only does not touch the unit row (leaves, xhat), so the vertex graph
        # cannot respond here; the first version of this check wrongly expected 7/128 there and was corrected (see UNIVERSALITY.md E1).
        vals = {i: rec.fmt(r2[i][0]) for i in r2 if r2[i][0]}
        check("Euler-dilaton at genus 2: [eps^1 u^0] = -1/16 on graph 2 (loop on g1), 1/128 on graph 4 (bridge g1-g1) (minus the record), 0 elsewhere (vertex graph included)",
              vals == {2: "-1/16", 4: "1/128"}, vals)
    log("UN_E1 %s" % ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED))


_E = {}
def base_E(g, idx):
    if g not in _E: _E[g] = {G['idx']: G['E'] for G in outer_graphs(g)}
    return _E[g][idx]


if __name__ == "__main__":
    main()
