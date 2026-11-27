# -*- coding: utf-8 -*-
"""ctl_cost.py -- what would an order-by-order genus-5 control cost?  Times graph_dp at --maxout 1, 2 and full
   on the three calibration graphs of the package (data: REF_DP in g5_validate.py: 10.81, 20.92, 3.49 s at full
   on the reference workstation) and scales the quoted 39-46 CPU hours of a full natural run.

   py -I -B ctl_cost.py [prime]
"""
import os, sys, time, json
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "cluster_genus5_quintic")
sys.path.insert(0, PKG)
REF_DP = {808: 10.81, 1109: 20.92, 4511: 3.49}


def main():
    prime = sys.argv[1] if len(sys.argv) > 1 else "89"
    import genus5_run as D
    from g5_engine import graph_dp
    spec = (prime, 5, "0,0", 8, "nat")
    t0 = time.time(); sp = D.get_split(spec)
    print("genus-5 Split set-up: %.0f s (%d inner graphs, K = %d)" % (time.time() - t0, len(sp.inner_graphs), sp.K))
    tot = {}
    for maxout in (1, 2, None):
        s = 0.0
        for t in sorted(REF_DP):
            t1 = time.time(); graph_dp(sp, t, maxout); dt = time.time() - t1; s += dt
            print("   graph %4d (E=%d) maxout=%-4s %8.2f s" % (t, sp.inner_graphs[t]['E'], maxout, dt))
        tot[maxout] = s
    print()
    print("sum over the three calibration graphs: maxout1 %.1f s, maxout2 %.1f s, full %.1f s" % (tot[1], tot[2], tot[None]))
    f = tot[None] / sum(REF_DP.values())
    print("machine factor f = %.2f (this machine / reference workstation)" % f)
    print("ratios: maxout1/full = %.3f, maxout2/full = %.3f, maxout1/maxout2 = %.3f" % (tot[1] / tot[None], tot[2] / tot[None], tot[1] / tot[2]))
    ct = json.load(open(os.path.join(PKG, "data", "cost_table_g5.json")))
    full = sum(v['est_dp'] for v in ct.values()) / 3600.0
    print("estimated CPU time of a full natural genus-5 run on the reference workstation: %.0f h (cost table est_dp)" % full)
    print("   => one --maxout 1 control: about %.1f CPU h; one --maxout 2 control: about %.1f CPU h (reference workstation)"
          % (full * tot[1] / tot[None], full * tot[2] / tot[None]))


if __name__ == "__main__":
    main()
