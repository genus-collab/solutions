# -*- coding: utf-8 -*-
"""ig_run.py -- driver: run the independent code (ig_lib) on a list of genus-5 inner graphs, store every finished job, resume.

  py -I -B ig_run.py --prime 89 --variant nat --list @sample_nat89.txt --nproc 13 --deadline 6000
  py -I -B ig_run.py --prime 89 --variant nat --list 149 --shards 10 --nproc 10
  py -I -B ig_run.py --prime 89 --variant nat --list 1,3,5,25,28 --brute 1,3,5,25,28     (colour-space brute force as well)

--list is a comma-separated list of inner graph indices or @file (one index per line, '#' comments allowed).
Results are appended to out/<tag>.jsonl, one JSON object per finished job: {"t", "shard", "res", "aut", "sec", "V", "E"}.
A job already in the file is skipped, so the script can be called again after its --deadline (no state is lost).
The tag is <variant>_<prime> unless --tag is given.
"""
import os, sys, json, time, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
sys.path.insert(0, HERE)
import ig_lib


def read_list(spec):
    if spec.startswith("@"):
        txt = open(os.path.join(HERE, spec[1:])).read()
        items = []
        for line in txt.splitlines():
            line = line.split("#")[0].strip()
            if line:
                items += [int(x) for x in line.replace(",", " ").split()]
        return items
    return [int(x) for x in spec.replace(",", " ").split()]


def done_jobs(path):
    out = set()
    if os.path.exists(path):
        for line in open(path, encoding="utf-8"):
            line = line.strip()
            if not line:
                continue
            d = json.loads(line)
            out.add((d["t"], tuple(d["shard"]) if d["shard"] else None))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prime", default="89")
    ap.add_argument("--variant", default="nat")
    ap.add_argument("--list", required=True)
    ap.add_argument("--nproc", type=int, default=13)
    ap.add_argument("--deadline", type=float, default=6000.0)
    ap.add_argument("--shards", type=int, default=0, help="split each listed graph into this many decoration shards")
    ap.add_argument("--brute", default="", help="also run the colour-space brute force on these indices")
    ap.add_argument("--tag", default="")
    ap.add_argument("--order", default="desc", choices=("desc", "asc"),
                    help="desc = most expensive first (best load balance); asc = cheapest first (use when the deadline may cut in "
                         "and the cheap strata must be complete)")
    a = ap.parse_args()
    assert a.nproc <= 13, "at most 13 worker processes"
    tag = a.tag or ("%s_%s" % (a.variant, a.prime))
    lst = read_list(a.list)
    brute = set(read_list(a.brute)) if a.brute else set()
    path = os.path.join(OUT, "%s.jsonl" % tag)
    log = open(os.path.join(OUT, "%s.log" % tag), "a", encoding="utf-8")

    def say(*x):
        s = " ".join(str(y) for y in x)
        print(s, flush=True)
        log.write(s + "\n")
        log.flush()

    have = done_jobs(path)
    cost = json.load(open(os.path.join(OUT, "cost_ind.json")))
    jobs = []
    for t in lst:
        shards = [(k, a.shards) for k in range(a.shards)] if a.shards else [None]
        for sh in shards:
            if (t, sh) in have:
                continue
            c = cost[str(t)][3] / (a.shards or 1)
            jobs.append((c, (a.prime, a.variant, t, list(sh) if sh else None, t in brute)))
    jobs.sort(key=(lambda x: -x[0]) if a.order == "desc" else (lambda x: x[0]))
    say("=" * 110)
    say("ig_run %s: %d graphs requested, %d jobs to do (%d already stored), nproc %d, deadline %.0f s, started %s"
        % (tag, len(lst), len(jobs), len(have), a.nproc, a.deadline, time.strftime("%Y-%m-%d %H:%M:%S")))
    say("   predicted cost of the queued jobs: %.1f CPU-h (proxy rate 5e-6 s)" % (sum(c for c, _ in jobs) * 5e-6 / 3600))
    if not jobs:
        say("nothing to do")
        return
    T0 = time.time()
    fh = open(path, "a", encoding="utf-8")
    n_done = 0
    stopped = [False]

    def gen():
        for c, j in jobs:
            if time.time() - T0 > a.deadline:
                stopped[0] = True
                say("   deadline reached: %d jobs not started" % (len(jobs) - n_done - a.nproc))
                return
            yield j

    import multiprocessing as mp
    with mp.Pool(processes=min(a.nproc, len(jobs))) as pool:
        for t, sh, o in pool.imap_unordered(ig_lib.job, gen()):
            n_done += 1
            fh.write(json.dumps(dict(t=t, shard=sh, res=o["res"], aut=o["aut"], sec=o["sec"], V=o["V"], E=o["E"],
                                     brute_equal=o.get("brute_equal"), sec_brute=o.get("sec_brute"))) + "\n")
            fh.flush()
            say("   t=%4d%s E=%2d V=%d: %.1f s, %d non-zero outer rows%s   [%d/%d, wall %.0f s]"
                % (t, "" if sh is None else " shard %d/%d" % tuple(sh), o["E"], o["V"], o["sec"], len(o["res"]),
                   "" if o.get("brute_equal") is None else ", brute force equal: %s (%.1f s)" % (o["brute_equal"], o["sec_brute"]),
                   n_done, len(jobs), time.time() - T0))
    fh.close()
    say("ig_run %s: %d jobs finished in %.0f s wall%s" % (tag, n_done, time.time() - T0,
                                                          "; DEADLINE STOP, call again to continue" if stopped[0] else "; queue empty"))


if __name__ == "__main__":
    main()
