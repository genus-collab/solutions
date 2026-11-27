# mem_probe.py PKGDIR t minutes [maxout]: peak working set of graph_dp(t) (Windows, ctypes; review tool, not part of the package)
import sys, os, time, threading, ctypes
from ctypes import wintypes
PKG = sys.argv[1]; t = int(sys.argv[2]); mins = float(sys.argv[3]); maxout = int(sys.argv[4]) if len(sys.argv) > 4 else None
sys.path.insert(0, PKG)
class PMC(ctypes.Structure):
    _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD), ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t), ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]
def mem():
    c = PMC(); c.cb = ctypes.sizeof(PMC)
    ctypes.windll.psapi.GetProcessMemoryInfo(ctypes.windll.kernel32.GetCurrentProcess(), ctypes.byref(c), c.cb)
    return c.PeakWorkingSetSize / 2**20, c.WorkingSetSize / 2**20, c.PeakPagefileUsage / 2**20
ST = dict(combos=0, maxstate=0, maxkeys=0)
import g5_engine
orig = g5_engine.combo_dp
def combo_dp(sp, steps, vecs, inner, outer, cap):
    # same as g5_engine.combo_dp, with the state size recorded
    mul = sp.mul; add = sp.add
    cur = {(): {0: [1] + [0] * (sp.K - 1)}}
    for v, conn, keep, loops in steps:
        new = {}; vec = vecs[v]
        for key, sub in cur.items():
            for c in range(5):
                vc = vec[c]
                if vc is None: continue
                facs = [(1 << e, inner[e][c][c], outer[e][c][c]) for e in loops]
                for e, pos, flip in conn:
                    k = key[pos]
                    facs.append((1 << e, inner[e][c][k], outer[e][c][k]) if flip else (1 << e, inner[e][k][c], outer[e][k][c]))
                items = []
                for S, val in sub.items():
                    x = mul(val, vc)
                    if x is not None: items.append((S, x))
                for bit, fi, fo in facs:
                    nxt = []
                    for S, x in items:
                        a = mul(x, fi)
                        if a is not None: nxt.append((S, a))
                        if g5_engine.POPC[S] < cap:
                            b = mul(x, fo)
                            if b is not None: nxt.append((S | bit, b))
                    items = nxt
                    ST['maxitems'] = max(ST.get('maxitems', 0), len(items))
                if not items: continue
                full = key + (c,); nk = tuple(full[k] for k in keep); dst = new.get(nk)
                if dst is None: dst = {}; new[nk] = dst
                for S, x in items:
                    y = dst.get(S); dst[S] = x if y is None else add(y, x)
        cur = new
        n = sum(len(d) for d in cur.values()); ST['maxstate'] = max(ST['maxstate'], n); ST['maxkeys'] = max(ST['maxkeys'], len(cur))
        if not cur: return {}
    ST['combos'] += 1
    return cur.get((), {})
g5_engine.combo_dp = combo_dp
T0 = time.time()
def watch():
    while True:
        time.sleep(20)
        el = time.time() - T0
        print("  %.0f s: peakWS %.0f MB, WS %.0f MB, peak commit %.0f MB, %s" % ((el,) + mem() + (ST,)), flush=True)
        if el > mins * 60: print("TIME CAP, stopping", flush=True); os._exit(0)
threading.Thread(target=watch, daemon=True).start()
from fq_fp import P2
sp = g5_engine.make_split(P2, 5, ("0", "0"), 8, "nat"); print("setup %.0f s, mem %s" % (time.time() - T0, mem()), flush=True)
G = sp.inner_graphs[t]; print("graph", t, "E", G['E'], "V", G['V'], "combos", __import__('math').prod(len(o) for o in g5_engine.decorations(sp, G)), flush=True)
t1 = time.time(); r = g5_engine.graph_dp(sp, t, maxout)
print("DONE graph %d in %.0f s, %d outer graphs; final mem %s; %s" % (t, time.time() - t1, len(r), mem(), ST), flush=True)
