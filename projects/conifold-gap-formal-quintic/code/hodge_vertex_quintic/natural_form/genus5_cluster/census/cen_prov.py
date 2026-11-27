# -*- coding: utf-8 -*-
"""cen_prov.py -- step 2 of the census: audit of the task decomposition and the provenance the result files lack.

For every cluster run: recompute the longest-processing-time assignment from data/cost_table_g5.json with my own code
(cen_lib.my_lpt) and compare it slice by slice with the inner graphs actually stored in task<I>of<N>.json; check that every
inner graph sits in exactly one task file; read the task logs for host, Python version, wall time and any "resuming: N done";
read logs/slurm/*.out for the SLURM job and array ids; read logs/validate.log for the machine factor f.

  py -I -B cen_prov.py <out.json>
"""
import os, sys, re, json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cen_lib as L


def read_task_logs(tag):
    out = {}
    d = L.logdir(tag)
    for f in sorted(os.listdir(d)):
        m = re.match(r"^task(\d+)of(\d+)\.log$", f)
        if not m:
            continue
        txt = open(os.path.join(d, f), encoding="utf-8", errors="replace").read()
        rec = dict(file=f, task=int(m.group(1)), ntasks=int(m.group(2)))
        hdr = re.findall(r"task (\d+) of (\d+), python (\S+), host (\S+)", txt)
        rec["headers"] = len(hdr)                        # > 1 means the task was started more than once
        rec["python"] = sorted(set(h[2] for h in hdr))
        rec["hosts"] = sorted(set(h[3] for h in hdr))
        rec["resuming"] = [int(x) for x in re.findall(r"resuming: (\d+) done", txt)]
        rec["setup_s"] = [int(x) for x in re.findall(r"set-up (\d+) s", txt)]
        comp = re.findall(r"COMPLETE: (\d+) inner graphs, (\d+) s", txt)
        rec["complete"] = [[int(a), int(b)] for a, b in comp]
        rec["wall_s"] = max([int(b) for a, b in comp], default=None)
        rec["n_assigned_in_log"] = [int(x) for x in re.findall(r"   (\d+) inner graphs assigned", txt)]
        rec["graph_lines"] = len(re.findall(r"   inner graph +\d+ \(E=", txt))
        out[int(m.group(1))] = rec
    return out


def slurm_index():
    d = os.path.join(L.EV, "logs", "slurm")
    idx = {}
    for f in sorted(os.listdir(d)):
        m = re.match(r"^([A-Z])_(\w+?)\.(\d+)(?:_(\d+))?\.out$", f)
        if not m:
            continue
        txt = open(os.path.join(d, f), encoding="utf-8", errors="replace").read(4000)
        t = re.search(r"genus5_run (g5_\S+):", txt)
        tag = t.group(1) if t else "?"
        key = (tag, m.group(2), m.group(3))
        e = idx.setdefault(key, dict(tag=tag, kind=m.group(2), jobid=m.group(3), array=[], files=0))
        e["files"] += 1
        if m.group(4) is not None:
            e["array"].append(int(m.group(4)))
    out = []
    for e in idx.values():
        a = e.pop("array")
        e["array_ids"] = ("%d..%d" % (min(a), max(a))) if a else None
        e["n_array"] = len(a)
        out.append(e)
    return sorted(out, key=lambda e: (e["tag"], e["kind"]))


def main(outpath):
    rep = dict(runs={}, slurm=slurm_index())
    v = open(os.path.join(L.EV, "logs", "validate.log"), encoding="utf-8", errors="replace").read()
    mf = re.search(r"machine factor f = ([\d.]+)", v)
    hv = re.search(r"g5_validate: python (\S+), (\d+) processes, host (\S+)", v)
    rep["validate"] = dict(machine_factor_f=float(mf.group(1)) if mf else None,
                           python=hv.group(1) if hv else None, nproc=int(hv.group(2)) if hv else None,
                           host=hv.group(3) if hv else None,
                           calibration=re.search(r"calibration \(wall \d+ s\): (\{.*\})", v).group(1),
                           wall_s=int(re.search(r"validation wall time (\d+) s", v).group(1)),
                           n_pass=len(re.findall(r"^PASS ", v, re.M)), n_fail=len(re.findall(r"^FAIL ", v, re.M)),
                           verdict=v.strip().splitlines()[-1])
    for tag in L.TAGS + L.IDTAGS:
        A = L.load_all_pieces(tag)
        per = A["per_file"]
        variant = L.variant_of(tag)
        r = dict(tag=tag, variant=variant, part_files=sorted(per), n_part_files=len(per),
                 cpu_seconds=sum(A["times"].values()), cpu_hours=sum(A["times"].values()) / 3600.0)
        if variant == "id":
            r["mode"] = "local workstation run (parts/<tag>/local.json), no SLURM array"
            r["inner_in_file"] = len(per["local.json"]["inner"])
            r["assignment_check"] = "n/a (single file)"
        else:
            nt = max(int(re.match(r"task(\d+)of(\d+)\.json", f).group(2)) for f in per)
            r["mode"] = "SLURM array, %d tasks" % nt
            r["ntasks"] = nt
            slices, load = L.my_lpt(L.costs_for(variant), nt)
            mism = []
            for i in range(nt):
                f = "task%dof%d.json" % (i, nt)
                got = sorted(per[f]["inner"]) if f in per else None
                want = sorted(slices[i])
                if got != want:
                    mism.append(dict(task=i, n_want=len(want), n_got=(None if got is None else len(got)),
                                     want_not_got=sorted(set(want) - set(got or []))[:10],
                                     got_not_want=sorted(set(got or []) - set(want))[:10]))
            r["assignment_check"] = dict(ntasks=nt, mismatched_tasks=mism, equal=(not mism),
                                         est_hours_heaviest=max(load) / 3600.0, est_hours_mean=sum(load) / nt / 3600.0)
            r["cost_column"] = "est_dp"
            logs = read_task_logs(tag)
            r["logs"] = dict(
                n_logs=len(logs), hosts=sorted(set(h for x in logs.values() for h in x["hosts"])),
                python=sorted(set(y for x in logs.values() for y in x["python"])),
                restarted_tasks=[i for i, x in logs.items() if x["headers"] > 1],
                resumed_tasks={str(i): x["resuming"] for i, x in logs.items() if x["resuming"]},
                wall_max_s=max(x["wall_s"] for x in logs.values()),
                wall_min_s=min(x["wall_s"] for x in logs.values()),
                wall_mean_s=sum(x["wall_s"] for x in logs.values()) / len(logs),
                wall_sum_h=sum(x["wall_s"] for x in logs.values()) / 3600.0,
                setup_s=sorted(set(s for x in logs.values() for s in x["setup_s"])),
                graph_lines_total=sum(x["graph_lines"] for x in logs.values()),
                assigned_total=sum(sum(x["n_assigned_in_log"][:1]) for x in logs.values()),
                all_complete=all(x["complete"] for x in logs.values()))
        # exactly one task file per inner graph
        r["coverage"] = dict(covered=len(A["pieces"]), missing=sorted(set(range(L.NIN)) - set(A["pieces"])),
                             n_duplicate_inner=len(A["dup"]), conflicting_duplicates=A["conflict"],
                             union_is_0_to_4554=(sorted(A["pieces"]) == list(range(L.NIN))))
        rep["runs"][tag] = r
        print("%-24s %-28s assign_equal=%s dup=%d missing=%d CPU=%.1f h wall_max=%s"
              % (tag, r["mode"], (r["assignment_check"] if isinstance(r["assignment_check"], str)
                                  else r["assignment_check"]["equal"]), len(A["dup"]),
                 len(r["coverage"]["missing"]), r["cpu_hours"], r.get("logs", {}).get("wall_max_s")))
        del A
    with open(outpath, "w") as f:
        json.dump(rep, f, indent=1)
    print("written", outpath)


if __name__ == "__main__":
    main(sys.argv[1])
