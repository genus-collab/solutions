"""Stage-2 critic: end-to-end check of the leg dictionary LD on actual TR correlators (genus 2 and 3 one-point functions).

For every frame X = frame(k,a) = z^k H^a + Sum_(j<k) f_j z^j H of K_0 (a in 0,2,3,4; k <= 3g-2; built from the PiC kernel by
cr2_legcore.pic), take its leg -Lambda(X) (log-free; cr2_legdict L4), pair it with omega_{g,1} of x = W = p^4/4 + 1/p, y = p
computed EXACTLY over Q(zeta5) by the allgenus/tr engine, and compare (-5)^(1-g) x (TR pairing) with the frozen one-point
function <X>_g obtained from the exact one-point table of tr/TR_MODEL.md s.3 (an independent computation: MVertex6 and the
xi-basis extraction).  ML_K predicts 0 whenever the exponent of <X>_g is negative.
Read-only imports: cr2_legcore (this folder); stage2/logleg/ll_core.py, which imports allgenus/tr (tr_field, tr_core, tr_rec, tr_v0).
Usage: py -I -B cr2_e2e.py      (stdlib only; writes nothing; well under a minute)
"""
import os
import sys
import time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', 'logleg')))
import cr2_legcore as C          # noqa: E402
import ll_core as LL            # noqa: E402  (read-only)

# exact frozen one-point functions at v = 1 (TR_MODEL s.3; complete lists for these genera): (b, k) -> value
ONEPT = {2: {(0, 1): Fr(1, 600), (1, 0): Fr(-1, 600), (4, 2): Fr(1, 1200), (3, 3): Fr(7, 5760), (2, 4): Fr(-5, 1152)},
         3: {(0, 1): Fr(1, 6300), (1, 0): Fr(-1, 6300), (4, 2): Fr(1, 12600), (3, 3): Fr(1, 37800), (2, 4): Fr(-1, 151200),
             (1, 5): Fr(53, 1209600), (0, 6): Fr(-1, 46080), (4, 7): Fr(5, 82944)}}
CW = (0, -1, 3, 2, 1)


def frame_coeffs(a, k):
    f = {}
    for i in range(k - 1, -1, -1):
        val = Fr(0)
        for (bb, kk, co) in [(a, k, Fr(1))] + [(1, j, f[j]) for j in range(i + 1, k)]:
            pc = C.pic(bb, kk - i)
            if pc:
                val += pc[0] * co                     # v = 1
        f[i] = -val
    return f


def main():
    t0 = time.time()
    F = LL.QzField()
    cur, tr = LL.setup(F, 2 * (6 * 3 - 5 + 2) + 6)
    legs = LL.Legs(F, cur.c)
    T = C.legs(8)
    nbad, ntest, nnz = 0, 0, 0
    for g in (2, 3):
        W = tr.omega(g, 1)
        for a in (0, 2, 3, 4):
            for k in range(3 * g - 1):
                f = frame_coeffs(a, k)
                leg = {}
                for (bb, kk, co) in [(a, k, Fr(1))] + [(1, j, f[j]) for j in range(k)]:
                    for (j, l, m), x in T[(bb, kk)].items():
                        leg[(j, l)] = leg.get((j, l), 0) - x * co / C.CB[bb]          # v = 1: v^m -> 1
                leg = {key: x for key, x in leg.items() if x}
                C.check(all(l == 0 for (_, l) in leg), "frame(%d,%d): log part" % (k, a))
                tr_val = F.zero
                for (j, l), x in leg.items():
                    r = LL.contract_last(W, [legs.mono(j)]).get((), F.zero)
                    tr_val = tr_val + r * F.of(x)
                C.check(tr_val.is_rational(), "TR pairing not rational")
                tr_val = tr_val.rational() * Fr(-5) ** (1 - g)
                frozen = ONEPT[g].get((a, k), Fr(0)) + sum(f[j] * ONEPT[g].get((1, j), Fr(0)) for j in range(k))
                e = 2 - 2 * g + Fr(4 * k + CW[a] - 4, 5)              # exponent of <frame(k,a)>_g
                ok = C.check(tr_val == frozen, "g=%d frame(%d,%d): TR %s frozen %s" % (g, k, a, tr_val, frozen))
                if e < 0 and e.denominator == 1:
                    ok &= C.check(frozen == 0, "g=%d frame(%d,%d): singular and non-zero" % (g, k, a))
                ntest += 1
                nbad += not ok
                if frozen:
                    nnz += 1
                    print("  g=%d frame(%d,%d): exponent %s, (-5)^(1-g) <leg, omega> = %s = frozen value" % (g, k, a, e, tr_val))
    print(" E2E: %d frame one-point functions (g = 2, 3), %d non-zero, mismatches %d; %.1f s" % (ntest, nnz, nbad, time.time() - t0))
    print("RESULT:", "ALL CHECKS PASSED" if not C.FAILS else "%d FAILURES" % len(C.FAILS))


if __name__ == "__main__":
    main()
