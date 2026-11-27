import os
import re
def graphs4():
    out=[]
    for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "hodge_collect_genus4.txt")):
        m=re.match(r'^(\d+)\s+\{([\d, ]*)\}\s+\{(.*?)\}\s+(\d+)\s+\d+\s+\d+',line)
        if not m: continue
        gen=tuple(int(x) for x in m.group(2).replace(' ','').split(',') if x)
        edges=[(int(a)-1,int(b)-1) for a,b in re.findall(r'\{(\d+), *(\d+)\}',m.group(3))]
        out.append((gen,edges,int(m.group(4))))
    return out
