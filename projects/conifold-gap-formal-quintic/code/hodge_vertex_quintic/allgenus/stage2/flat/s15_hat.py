# -*- coding: utf-8 -*-
"""s15_hat.py -- the project's HATTED outer graph sum R'.(M.omega^) at genus 3 with the natural form M = M_0(tau(u) zhat) taken
EXACTLY (not truncated at u^3 as in g3_nat_<p>_nat.json), to higher u-order, with the project's own machinery (read-only imports):
fq3_par.setup (rows, conifold frame, symbol, M_0, HatRowsU), fq3_msum.GraphSummer6 (outer sum), fq3_mvertex.MVertex6.

Homogeneity (PROVED, s15 note s.2 step 3): the R-action is graded by z-degree = cohomological degree, so the degree-D part of
M_0(tau zhat).omega^ is tau^D (M_0.omega^)_D; the vertex tensor of psi-deficit D = 3g-3+n-Sum d is tau(u)^D x the frozen tensor.
This replaces the eps-interpolation (UVertex) of the stored pass; it is checked against the stored tables (A3 [u^0..u^3], valid part).
Modes: 'nat' (M = M_0(tau zhat)), 'M0' (M = M_0 constant: must reproduce g3_state_<p>.json A_M0 on all stored orders).
Run:  py -I -B s15_hat.py <89|61> <K> <prec_abs> <nat|M0> [alpha,beta]      Output: <scratch>/hat_g3_<p>_<mode>_K<K>[_ab..].json, log logs/s15_hat_<p>_<mode>_K<K>[_ab..].log.
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
if PROJ not in sys.path: sys.path.insert(0, PROJ)
from fractions import Fraction as Fr
from fq_fp import P1, P2
from fq_series import Ser
from fq3_par import setup
from fq3_msum import GraphSummer6
from fq3_mvertex import MVertex6, KMAX
from fq3_graphs import outer_graphs3
sys.path.insert(0, HERE)
from s15_flat import SCRATCH


def tau_coeffs(n):
    inv = [Fr(4, 5 * k + 4) for k in range(n)]; t = [Fr(1)]
    for k in range(1, n): t.append(-sum(inv[j] * t[k - j] for j in range(1, k + 1)))
    return t


class MNat(object):
    """duck-typed M for GraphSummer6: only Minv_ser (rows_prime6) and M (EdgeLeaf.T_alt, unused) are read"""
    def __init__(self, R, M0, taupow, K):
        p = R.p; self.is_identity = False; self.M = M0.M
        self.Minv_ser = [[[Ser(R, 0, [M0.Minv[k][j][i] * taupow[k][n] % p for n in range(K)]) for i in range(5)] for j in range(5)]
                         for k in range(KMAX + 1)]


class ScaledVertex(object):
    def __init__(self, R, M0, taupow, K):
        self.R = R; self.mv = MVertex6(R, M0); self.tp = taupow; self.K = K; self.cache = {}

    def tensor(self, g, d):
        key = (g, tuple(d))
        if key in self.cache: return self.cache[key]
        T = self.mv.tensor(g, tuple(d)); D = 3 * g - 3 + len(d) - sum(d); R = self.R; K = self.K
        s = Ser(R, 0, list(self.tp[D]))
        out = [s.scale(t) if t else Ser(R, K, []) for t in T]
        self.cache[key] = out
        return out


def main():
    pr, K, prec, mode = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    p = P2 if pr == "89" else P1
    ab = sys.argv[5].split(",") if len(sys.argv) > 5 else ["0", "0"]
    ftag = "" if ab == ["0", "0"] else "_ab" + "_".join(x.replace("/", "o").replace("-", "m") for x in ab)
    LOG = open(os.path.join(HERE, "logs", "s15_hat_%s_%s_K%d%s.log" % (pr, mode, K, ftag)), "w", encoding="ascii")
    def log(s): print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    log("=" * 100); log("s15_hat: genus 3 hatted outer sum, p = 2^%s-1, K = %d, prec_abs = %d, mode %s, frame %s [%s]" % (pr, K, prec, mode, ab, time.strftime("%Y-%m-%dT%H:%M:%S")))
    t0 = time.time(); st = setup((p, str(Fr(ab[0])), str(Fr(ab[1])), prec, K)); R = st['R']; H = st['H']; M0 = st['M0']
    log("  setup %.0f s; HatRowsU min_prec %d (need >= K = %d), val_ok %s" % (time.time() - t0, H.min_prec, K, H.val_ok))
    assert H.min_prec >= K
    tq = [R.of(x) for x in tau_coeffs(K)] if mode == 'nat' else [1] + [0] * (K - 1)
    tp = [[1] + [0] * (K - 1)]
    for k in range(1, 3 * 3 + 2): tp.append([sum(tp[-1][i] * tq[n - i] for i in range(n + 1)) % p for n in range(K)])
    Mo = MNat(R, M0, tp, K); vt = ScaledVertex(R, M0, tp, K); GS = GraphSummer6(H, Mo, vt)
    out = {}
    for G in outer_graphs3():
        t = time.time(); s = GS.graph(G); co = [s.coef(k) if k < s.prec else None for k in range(K)]; out[G['idx']] = co
        log("  graph %2d (E=%d): %s (%.0f s)" % (G['idx'], G['E'], [R.fmt(x) if x is not None else '?' for x in co], time.time() - t))
    fn = os.path.join(SCRATCH, "hat_g3_%s_%s_K%d%s.json" % (pr, mode, K, ftag))
    json.dump(dict(p=p, K=K, prec_abs=prec, mode=mode, A={str(k): v for k, v in out.items()}), open(fn, "w"))
    log("  written %s ; total %.0f s" % (fn, time.time() - t0)); LOG.close()


if __name__ == "__main__":
    main()
