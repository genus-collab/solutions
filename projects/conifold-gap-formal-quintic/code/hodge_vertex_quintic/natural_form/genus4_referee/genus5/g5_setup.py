# -*- coding: utf-8 -*-
"""g5_setup.py -- (a) the patched Split (fast graphs) reproduces stored genus-4 per-inner-graph results; (b) genus-5 setup test.

(a) genus 4, 2^89-1, natural M, K = 6: for the inner graphs whose stored time (parts/g4_89_nat_run.json, read only) is below a
    threshold, Split.graph on the fast graph list equals the stored {outer graph: coefficients} after matching the graphs through
    the original canonical form fq3_graphs.canon.
(b) Split(p, 5, (0, 0), K = 8, 'nat'): set-up time, peak memory, the assertions inside Split.__init__ (symbol valuations, hatted-row
    purity H.val_ok, H.min_prec >= K, M_0 M_0^-1 = 1 through zhat^12), frozen_ok, the u-valuations of E_out, vertex-table times.
Log: g5_setup.log.   Usage: py -I g5_setup.py [a|b|ab] [prime 89|61] [threshold seconds for (a)]
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from g5_lib import make_split, ROOT
from fq_fp import P1, P2

LOG = open(os.path.join(HERE, "g5_setup.log"), "a", encoding="utf-8"); FAILED = []


def log(s):
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


def check(name, ok, detail=""):
    log(("PASS " if ok else "FAIL ") + name + ((": " + str(detail)) if detail else ""))
    if not ok: FAILED.append(name)


def peak_mb():
    try:
        import ctypes
        from ctypes import wintypes
        class PMC(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD), ("PeakWorkingSetSize", ctypes.c_size_t),
                        ("WorkingSetSize", ctypes.c_size_t), ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                        ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]
        pmc = PMC(); pmc.cb = ctypes.sizeof(PMC)
        k32 = ctypes.WinDLL("kernel32"); k32.GetCurrentProcess.restype = wintypes.HANDLE
        psapi = ctypes.WinDLL("psapi"); psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(PMC), wintypes.DWORD]
        psapi.GetProcessMemoryInfo(k32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb)
        return pmc.PeakWorkingSetSize / 2.0 ** 20
    except Exception:
        try:
            import resource
            return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
        except Exception:
            return float("nan")


mode = sys.argv[1] if len(sys.argv) > 1 else "ab"; prime = sys.argv[2] if len(sys.argv) > 2 else "89"
thr = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0; p = P2 if prime == "89" else P1
log("=" * 100); log("g5_setup mode %s, prime 2^%s-1, python %s" % (mode, prime, sys.version.split()[0]))

if "a" in mode:
    from fq3_graphs import legged_graphs, canon
    old = legged_graphs(4, 0, 9); old_pos = {G['key']: i for i, G in enumerate(old)}
    stored = json.load(open(os.path.join(ROOT, "parts", "g4_%s_nat_run.json" % prime)))
    t0 = time.time(); sp = make_split(p, 4, (0, 0), 6, 'nat'); log("   genus-4 setup %.1f s (fast graphs)" % (time.time() - t0))
    to_old = [old_pos[canon(*G['key'])] for G in sp.inner_graphs]
    check("genus 4: the fast list is a relabelling of the original 379 graphs", sorted(to_old) == list(range(379)))
    todo = [t for t in range(379) if stored['times'][str(to_old[t])] < thr]; bad = []; t0 = time.time()
    for t in todo:
        mine = {to_old[pos]: c for pos, c in sp.graph(t).items()}
        ref = {int(k): v for k, v in stored['res'][str(to_old[t])].items()}
        mine = {k: v for k, v in mine.items() if any(v)}; ref = {k: [x % p for x in v] for k, v in ref.items() if any(x % p for x in v)}
        if mine != ref: bad.append(t)
    check("genus 4, 2^%s-1, natural M: Split.graph on the fast graph list = stored per-inner-graph tables for the %d inner graphs with stored time < %.1f s (%.0f s)"
          % (prime, len(todo), thr, time.time() - t0), not bad, bad[:10])

if "b" in mode:
    t0 = time.time(); sp = make_split(p, 5, (0, 0), 8, 'nat'); ts = time.time() - t0
    log("   genus-5 setup %.1f s; K = %d, KZ = %d, %d inner graphs, %d edge keys (a, b), leaves e = 2..%d; peak memory %.0f MB"
        % (ts, sp.K, sp.KZ, len(sp.inner_graphs), len(sp.ER), sp.KZ + 1, peak_mb()))
    check("genus 5: Split.__init__ completed (asserts: symbol valuations, hatted-row purity H.val_ok, H.min_prec >= 8, M_0 M_0^-1 = 1 mod zhat^13)", True)
    check("genus 5: [u^0] E_R(a, b) = Eh_0(a, b) for all (a, b) with a + b <= 11 (frozen_ok)", sp.frozen_ok is True, sp.frozen_ok)
    hist = {}
    for k, v in sp.valO.items(): hist[v] = hist.get(v, 0) + 1
    log("   u-valuation of E_out(a, b) (min over colours): histogram %s" % sorted(hist.items()))
    hist2 = {}
    for k in sp.EO:
        for row in sp.EO[k]:
            for x in row: v = sp.val(x); hist2[v] = hist2.get(v, 0) + 1
    log("   u-valuation of the entries E_out(a, b)[i][j]: histogram %s" % sorted(hist2.items()))
    log("   u-valuation of the leaves T_R(e)[j]: %s" % {e: sorted(set(sp.val(x) for x in sp.TR[e])) for e in sorted(sp.TR)})
    nz = sum(1 for k in sp.EM for row in sp.EM[k] for x in row if x is not None)
    log("   non-zero entries of E_M: %d of %d" % (nz, 25 * len(sp.EM)))
    need = sorted(set((G['gv'][v], len(G['hsAt'][v])) for G in sp.inner_graphs for v in range(G['V'])))
    log("   vertex types (g_v, valence) occurring: %d" % len(need))
    tt = []
    for (g, hs) in need:
        t1 = time.time(); tab = sp.vertex_table(g, hs); tt.append((time.time() - t1, g, hs, len(tab)))
    log("   all vertex tables: %.1f s; slowest (s, g_v, valence, sorted decorations): %s" % (sum(x[0] for x in tt), [(round(a, 1), b, c, d) for a, b, c, d in sorted(tt, reverse=True)[:8]]))
    import math
    def nperm(hp):
        r = math.factorial(len(hp)); cnt = {}
        for x in hp: cnt[x] = cnt.get(x, 0) + 1
        for m in cnt.values(): r //= math.factorial(m)
        return r
    opt = {(g, hs): sum(nperm(hp) for hp in sp.vertex_table(g, hs)) for (g, hs) in need}
    log("   decorations per vertex (ordered half-edge exponents), by (g_v, valence): %s" % sorted(opt.items()))
    log("   peak memory after the vertex tables: %.0f MB" % peak_mb())
log("G5_SETUP %s" % ("ALL CHECKS PASSED" if not FAILED else "FAILURES: %s" % FAILED))
