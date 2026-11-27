import json, os, sys
H = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
for p in ("61", "89"):
    for c in ("circ", "circpi", "circb", "circstar"):
        f = os.path.join(H, "g3_u3_%s_%s.json" % (p, c))
        d = json.load(open(f))
        print(p, c, {k: (len(v) if hasattr(v, '__len__') else v) for k, v in d.items()})
        for k in ("A1", "A2", "A3"):
            if k in d:
                print("   ", k, sorted(int(i) for i in d[k]), set(len(v) for v in d[k].values()))
    g = json.load(open(os.path.join(H, "g3_state_%s.json" % p)))
    print(p, "state keys", {k: (len(v) if hasattr(v, '__len__') else v) for k, v in g.items()})
    q = json.load(open(os.path.join(H, "g3_quad_%s.json" % p)))
    print(p, "quad keys", {k: (len(v) if hasattr(v, '__len__') else v) for k, v in q.items()})
    cj = json.load(open(os.path.join(H, "g3_cubic_%s.json" % p)))
    print(p, "cubic keys", {k: (len(v) if hasattr(v, '__len__') else v) for k, v in cj.items()})
