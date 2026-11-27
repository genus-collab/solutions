# -*- coding: utf-8 -*-
"""rv_mut_g3.py -- mutation tests against the genus-2/3 regression of genus4_split.py.

For every mutant (rv_common.MUTANTS) the ORIGINAL function genus4_split.regress is run with get_split replaced by a constructor of the
mutated copy of fq4_split.Split (both primes; genus 3: variants id, M0, nat against the stored tables of the old vertex-tensor code;
genus 2: nat, M0 against the values of the notes).  A mutant is DETECTED if some regression line FAILs or the code crashes.
Usage: py -I rv_mut_g3.py [nproc (default 5)] [name ...]     -> out/mut_g3.json, out/mut_g3.log
"""
import os, sys, json, time, traceback
sys.dont_write_bytecode = True; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rv_common import PROJ, OUT, MUTANTS, write_mutants, mutant_module
from fractions import Fraction as Fr


def run_one(name):
    t0 = time.time(); lines = []
    try:
        import genus4_split as drv
        from fq_fp import P1, P2
        mod = mutant_module(name); cache = {}
        def get_split(spec):
            sp = cache.get(spec)
            if sp is None:
                prime, g, a, b, K, variant = spec
                sp = mod.Split(P2 if prime == "89" else P1, g, (Fr(a), Fr(b)), K, variant); cache[spec] = sp
            return sp
        drv.get_split = get_split
        drv.regress(None, lambda *x: lines.append(" ".join(str(y) for y in x)))
        crash = None
    except Exception as e:
        crash = "%s: %s" % (type(e).__name__, str(e)[:300]); lines.append("CRASH " + crash)
        lines.append(" | ".join(traceback.format_exc().strip().splitlines()[-4:])[:400])
    fails = [l for l in lines if l.startswith("FAIL")]; passes = [l for l in lines if l.startswith("PASS")]
    return dict(name=name, desc=MUTANTS[name][1], crash=crash, nfail=len(fails), npass=len(passes), lines=lines, seconds=time.time() - t0)


def main():
    args = sys.argv[1:]; nproc = int(args[0]) if args and args[0].isdigit() else 5
    names = [a for a in args if not a.isdigit()] or sorted(MUTANTS)
    write_mutants(names)
    import multiprocessing as mp
    t0 = time.time(); res = {}
    path = os.path.join(OUT, "mut_g3.json")
    if os.path.exists(path): res = json.load(open(path))
    LOG = open(os.path.join(OUT, "mut_g3.log"), "a", encoding="utf-8")
    def log(s):
        print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
    log("=" * 100); log("rv_mut_g3: %d mutants, %d processes" % (len(names), nproc))
    with mp.Pool(processes=min(nproc, len(names))) as pool:
        for r in pool.imap_unordered(run_one, names):
            res[r['name']] = r
            verdict = "CRASH (detected)" if r['crash'] else ("DETECTED" if r['nfail'] else "NOT DETECTED")
            if r['name'] == 'baseline': verdict = "control ok" if (not r['crash'] and not r['nfail']) else "CONTROL BROKEN"
            log("%-22s %-18s %2d FAIL / %2d PASS  (%.0f s)  %s" % (r['name'], verdict, r['nfail'], r['npass'], r['seconds'], r['desc']))
            for l in r['lines']:
                if l.startswith("FAIL") or l.startswith("CRASH"): log("      " + l[:230])
            json.dump(res, open(path, "w"), indent=1)
    log("rv_mut_g3 done in %.0f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
