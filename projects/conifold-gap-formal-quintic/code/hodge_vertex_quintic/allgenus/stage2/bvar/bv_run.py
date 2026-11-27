# -*- coding: utf-8 -*-
"""bv_run.py -- Stage 2 falsifier (c): B-universality.  For the systems of bv_sys.SYSTEMS (B = prod b_a perturbed) recompute tau, v(u),
M = M_0(tau zhat), the Theorem-U gauge, its canonical normalisation and the hatted rows; run the outer graph sum; report regularity
graph by graph and the vertex law c_g (T'')^(2-2g), T'' = -[z^-1 e_0] P Y^C_model = c0 T, T = Int_0^u f, c0 = P_11(0,0) = 1/lambda.
Modes (outputs go to --out, default the scratch folder; logs < 250 lines, one per call):
  prep SYS[,SYS..]                      Theorem U, (G) mod u^6 over Q, symplecticity of the rows, regressions (R0)
  run  SYS g prime [--norm canon|stored] [--ctrl none|dress0|dress2|staleM] [--maxout m] [--nproc n]
Controls: stored = P_00(0,z) = 1 normalisation (visible from genus 4); dress0/dress2 = rows dressed by 1 + (1/7) u^j zhat E_(2,1)
(flat-regular, symplectic, violates annihilation at zhat^1); staleM = correct rows but E_M with the QUINTIC tau (M not recomputed).
Run: py -I -B bv_run.py prep R0,S1,S2,S3,S4,S5 ;  py -I -B bv_run.py run S1 3 89
"""
import os, sys, json, time, argparse, datetime
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path: sys.path.insert(0, HERE)
import bv_sys
from bv_sys import SYSTEMS, Q
from fq_fp import FpRing
from fq_series import Ser
from fq3_rec import Rec
from fq3_msum import tc_expansion3
from un_core import Arith, PRIMES, lift_rows, symplectic_check, left_dress
from un_engine import get_base, get_nat
from un_gauge import check_G
from un_e4_cases import normalise, p11_series
SCR = os.path.join(HERE, "scratch")
CG = {2: Fr(1, 1200), 3: Fr(1, 25200), 4: Fr(1, 180000)}
_G = {}
OUT = [SCR]; LOGF = [None]; FAILED = []


def log(*x):
    s = " ".join(str(y) for y in x); print(s, flush=True)
    if LOGF[0]:
        with open(LOGF[0], "a", encoding="utf-8") as f: f.write(s + "\n")


def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)


def gauge_of(name):
    if name not in _G: _G[name] = bv_sys.gauge(SYSTEMS[name], 12, 24)
    return _G[name]


def normalise_sp(res0, norm):
    """un_e4_cases.normalise, then the constant scalar sigma = 1/p_{0,4}(0) making P symplectic (Prop. S at u = 0:
    p_00 p_03 = p_01 p_02 = p_04^2 = C at u = 0, z^0; sigma = 1 for R0-S4, sigma = 32 for S5)"""
    res, m = normalise(res0, norm); sig = 1 / res['p'][0][4].coef(0)
    if sig != 1:
        res = dict(res); res['p'] = [[(s.scale(sig) if s is not None else None) for s in row] for row in res['p']]
    return res, m, sig


def rows_of(name, norm, prime, g, K):
    res0 = gauge_of(name); res, m, sig = normalise_sp(res0, norm)
    nr = get_nat(prime, g, ("0", "0"), K); R = nr['R']; KZ = 3 * g - 3
    Qh = bv_sys.rows(res, nr['Rh'], R, KZ, K)
    return res, m, Qh, nr


def prep(names):
    for name in names:
        t0 = time.time(); res0 = gauge_of(name); d = SYSTEMS[name]
        log("-" * 100); log("system %s: %s   (gauge over Q, k <= 12, O(u^24): %.1f s)" % (name, d[6], time.time() - t0))
        log("   B = %s" % res0['B'].fmt(2)); log("   tau = %s" % res0['tau'].fmt(3)); log("   v = %s" % res0['v'].fmt(4))
        hol = all(x is None or x >= 0 for inf in res0['info'] for x in inf['vals'])
        check("   Theorem U: every p_{k,b}, k <= 12, is a power series in u (stored normalisation)", hol)
        for norm in ("canon", "stored"):
            res, m, sig = normalise_sp(res0, norm); P11 = p11_series(res['p'])
            p0 = [res['p'][0][b].coef(0) for b in range(5)]
            check("   %s: symplectic scalar sigma = %s; then p_00 p_03 = p_01 p_02 = p_04^2 = 1 at u = 0" % (norm, sig), p0[0] * p0[3] == p0[1] * p0[2] == p0[4] ** 2 == 1, p0)
            ok, c, bad, nchk, rng = check_G(res, 12, 6)
            const = all(P11[k] == 0 for k in P11 if k > 0)
            log("   %s: P_11(0,z) = %s" % (norm, " + ".join("(%s) z^%d" % (P11[k], k) for k in sorted(P11))))
            check("   %s: (G) P Y^C_model(v(u)) = P_11(0,z) Y^C_system on z^%d..z^%d mod u^6 (%d z-powers, exact over Q)" % (norm, rng[0], rng[1], nchk), ok, bad[:2])
            log("   %s: annihilation (P Y^C_model)_+ = 0: %s" % (norm, "YES (P_11(0,z) constant)" if const else "NO (P_11(0,z) has z^5, z^10 terms)"))
        log("   T = Int f = %s;  c0 = P_11(0,0) = %s (canonical, symplectic), so T'' = c0 T" % (res0['T'].fmt(3), p11_series(normalise_sp(res0, "canon")[0]['p'])[0]))
        for prime in ("89", "61"):
            for g in (2, 3, 4):
                K = 2 * g - 2; res, m, Qh, nr = rows_of(name, "canon", prime, g, K); ar = Arith(nr['p'], K, 1); QL = lift_rows(ar, Qh)
                sb = symplectic_check(ar, QL, 3 * g - 3, nr['eta_inv'])
                check("   p=2^%s-1, genus %d: canonical rows symplectic mod (zhat^%d, u^%d)" % (prime, g, 3 * g - 2, K), not sb, sb[:3])
                if name == "R0":
                    ref = nr['Qh']; same = all(Qh[a][i][n] == ref[a][i][n] for a in range(5) for i in range(5) for n in range(3 * g - 2))
                    check("   p=2^%s-1, genus %d REGRESSION: R0 rows = project rows (25 x %d x %d residues)" % (prime, g, 3 * g - 2, K), same)
                    base = get_base(prime, g, ("0", "0"), K); em = bv_sys.em_table(base, bv_sys.tau_list(res, nr['R'], K))
                    check("   p=2^%s-1, genus %d REGRESSION: E_M from the recomputed tau = base E_M" % (prime, g), em == base.EM)


def expand(G, Ts, R, g, K):
    """u^(2-2g) G(u) in powers of T'' (T'' = kap u + O(u^2), Q-series): coefficients of T''^(2-2g) .. T''^-1"""
    kap = Ts.coef(1); Tn = Ts.scale(1 / kap); TnR = Ser(R, Tn.v, [R.of(x) for x in Tn.c]).truncate(K + 4)
    b = tc_expansion3(G[:K], TnR, R, g)[:2 * g - 2]
    # b[j] multiplies Tn^(j-2g+2) = kap^(2g-2-j) T''^(j-2g+2)
    return [b[j] * R.pow(R.of(kap), 2 * g - 2 - j) % R.p for j in range(2 * g - 2)]


def run(name, g, prime, norm, ctrl, maxout, nproc, kext=0):
    from bv_engine import run_single, run_pool
    K = 2 * g - 2 + kext; KZ = 3 * g - 3; npol = 2 * g - 2
    tag = "%s_g%d_%s_%s_%s%s%s" % (name, g, prime, norm, ctrl, "" if maxout is None else "_out%d" % maxout, "" if not kext else "_k%d" % K)
    LOGF[0] = os.path.join(OUT[0], "run_%s.log" % tag)
    log("=" * 100); log("bv_run %s   [%s]" % (tag, datetime.datetime.now().isoformat(timespec='seconds')))
    log("system %s: %s;  normalisation %s;  control %s" % (name, SYSTEMS[name][6], norm, ctrl))
    res, m, Qh, nr = rows_of(name, norm, prime, g, K); R = nr['R']; p = nr['p']; rec = Rec(R); ar = Arith(p, K, 1)
    QL = lift_rows(ar, Qh); tl = bv_sys.tau_list(res, R, K)
    nd = sum(1 for a in range(5) for i in range(5) for n in range(KZ + 1) if Qh[a][i][n] != nr['Qh'][a][i][n])
    log("   rows: %d of %d residue slots (a, i, zhat^n) differ from the project rows (quintic, frame (0,0)); u-window u^0..u^%d" % (nd, 25 * (KZ + 1), K - 1))
    if ctrl in ("dress0", "dress2"):
        j = int(ctrl[-1]); X = {(2, 1, 1): ar.mono(R.of(Fr(1, 7)), j)}
        QL = left_dress(ar, QL, X, KZ); log("   rows dressed by 1 + (1/7) u^%d zhat E_(2,1): flat-regular (mu(1,1) = 0), X^2 = 0, violates annihilation at zhat^1" % j)
    if ctrl == "staleM": tl = None; log("   CONTROL: E_M built with the QUINTIC tau (M not recomputed)")
    sb = symplectic_check(ar, QL, KZ, nr['eta_inv'])
    check("   rows symplectic mod (zhat^%d, u^%d)" % (KZ + 1, K), not sb, sb[:3])
    spec = (prime, g, K, tl, QL, maxout); t0 = time.time()
    if g <= 3: v, tab = run_single(spec); cpu = time.time() - t0
    else: v, tab, cpu = run_pool(spec, nproc, log, sid=tag)
    log("   graph sum: %d inner graphs, wall %.0f s, CPU %.0f s" % (len(v.inner_graphs), time.time() - t0, cpu))
    outer = [dict(E=G['E'], gv=list(G['gv'])) for G in v.inner_graphs]
    shown = [i for i in range(1, len(outer)) if maxout is None or outer[i]['E'] <= maxout]
    rows = []; bad = {}
    for i in shown:
        c = tab.get(i, [0] * K)[:npol]
        if any(c): bad[i] = c
        rows.append((i, outer[i]['E'], outer[i]['gv'], [rec.fmt(x) for x in c]))
    if not bad: log("   RESULT: all %d outer graphs with edges%s are O(u^%d): REGULAR" % (len(shown), "" if maxout is None else " (<= %d edges)" % maxout, npol))
    else:
        log("   RESULT: %d of %d outer graphs with edges are SINGULAR" % (len(bad), len(shown)))
        for i in sorted(bad): log("      outer graph %3d (E=%d, gv=%s): [u^0..u^%d] = %s" % (i, outer[i]['E'], outer[i]['gv'], npol - 1, [rec.fmt(x) for x in bad[i]]))
    v0 = tab.get(0, [0] * K)
    log("   vertex graph [u^0..u^%d] = %s" % (K - 1, [rec.fmt(x) for x in v0]))
    c0 = p11_series(res['p'])[0]; Tpp = res['T'].scale(c0)
    V = expand(v0, Tpp, R, g, K)
    okv = V[0] == R.of(CG[g]) and not any(V[1:])
    log("   vertex graph in T'' = c0 T (c0 = %s): %s  -> %s" % (c0, " + ".join("(%s) T''^%d" % (rec.fmt(V[k]), k - npol) for k in range(npol)),
        "c_%d (T'')^%d + 0" % (g, -npol) if okv else "NOT c_g T''^(2-2g) + 0"))
    tcq = res['C'].tc.truncate(K + 4); Vt = expand(v0, tcq, R, g, K); Vv = expand(v0, res['v'], R, g, K)
    log("   [control] in the quintic t_c: %s" % " + ".join("(%s) t_c^%d" % (rec.fmt(Vt[k]), k - npol) for k in range(npol)))
    log("   [control] in the model v(u): %s" % " + ".join("(%s) v^%d" % (rec.fmt(Vv[k]), k - npol) for k in range(npol)))
    if maxout is None:
        tot = [sum(tab.get(i, [0] * K)[k] for i in range(len(outer))) % p for k in range(K)]
        Vtot = expand(tot, Tpp, R, g, K)
        log("   total, in T'': %s" % " + ".join("(%s) T''^%d" % (rec.fmt(Vtot[k]), k - npol) for k in range(npol)))
    ref = os.path.join(OUT[0], "tab_%s.json" % tag.replace(name + "_", "R0_", 1).replace("_" + norm + "_" + ctrl, "_canon_none"))
    if name != "R0" and os.path.exists(ref):
        rt = {int(i): c for i, c in json.load(open(ref))['table'].items()}
        dif = [i for i in range(len(outer)) if (maxout is None or outer[i]['E'] <= maxout) and [x % p for x in tab.get(i, [0] * K)] != [x % p for x in rt.get(i, [0] * K)]]
        log("   compared with R0 (same g, prime, window): %d of %d outer graphs in range differ (u^0..u^%d)%s" % (len(dif), len([i for i in range(len(outer)) if maxout is None or outer[i]['E'] <= maxout]), K - 1, (": " + str(dif[:12])) if dif else ""))
    json.dump(dict(tag=tag, system=name, g=g, prime=prime, norm=norm, ctrl=ctrl, maxout=maxout, cpu=cpu, rows=rows,
                   vertex=[rec.fmt(x) for x in V], vertex_ok=okv, singular=sorted(bad), table={str(i): c for i, c in sorted(tab.items())}),
              open(os.path.join(OUT[0], "tab_%s.json" % tag), "w"))
    log("BV_RUN %s: %s; edge graphs %s; vertex law %s" % (tag, "checks passed" if not FAILED else "FAILURES %s" % FAILED,
        "all regular" if not bad else "singular %s" % sorted(bad), "holds" if okv else "FAILS"))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("mode"); ap.add_argument("systems"); ap.add_argument("g", nargs="?", type=int, default=2)
    ap.add_argument("prime", nargs="?", default="89"); ap.add_argument("--norm", default="canon"); ap.add_argument("--ctrl", default="none")
    ap.add_argument("--maxout", type=int, default=None); ap.add_argument("--nproc", type=int, default=4); ap.add_argument("--kext", type=int, default=0); ap.add_argument("--out", default=SCR)
    a = ap.parse_args(); OUT[0] = a.out; os.makedirs(a.out, exist_ok=True)
    if a.mode == "prep":
        LOGF[0] = os.path.join(a.out, "prep_%s.log" % a.systems.replace(",", "_"))
        log("=" * 100); log("bv_run prep   [%s]" % datetime.datetime.now().isoformat(timespec='seconds'))
        prep(a.systems.split(",")); log("BV_PREP %s" % ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED))
    else:
        run(a.systems, a.g, a.prime, a.norm, a.ctrl, a.maxout, a.nproc, a.kext)


if __name__ == "__main__":
    main()
