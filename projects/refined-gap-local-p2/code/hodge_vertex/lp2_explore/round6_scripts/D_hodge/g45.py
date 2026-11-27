# Fresh per-colouring regularity at genus 4, 5 on graphs with rescaled vertices (rescaled-rescaled edges, loops, genus-0 vertices). Evidence only.
import sys,time
from fractions import Fraction as Fr
from itertools import product
from rows import *
from graphs_check import Vtx,graph_value,recon
GR=[
 ("g4: (3)-(1)",[3,1],[(0,1)]),
 ("g4: (3)+loop",[3],[(0,0)]),
 ("g4: (2)-(2)",[2,2],[(0,1)]),
 ("g4: (2)+loop-(1)",[2,1],[(0,0),(0,1)]),
 ("g4: (2)-(1)+loop",[2,1],[(0,1),(1,1)]),
 ("g4: (2)+2 loops",[2],[(0,0),(0,0)]),
 ("g4: (2)=(1) double edge",[2,1],[(0,1),(0,1)]),
 ("g4: (2)-(0 loop)-? (2)-(0)-(1) w/ loop on 0 impossible; (2)-(0)=(1)",[2,0,1],[(0,1),(1,2),(1,2)]),
 ("g4: (1)-(2)-(1)",[1,2,1],[(0,1),(1,2)]),
 ("g5: (2)=(2) double edge",[2,2],[(0,1),(0,1)]),
 ("g5: (2)+loop-(2)",[2,2],[(0,0),(0,1)]),
 ("g5: (2)-(0)-(2) with loop on (0)",[2,0,2],[(0,1),(1,2),(1,1)]),
 ("g5: (4)-(1)",[4,1],[(0,1)]),
 ("g5: (3)-(2)",[3,2],[(0,1)]),
 ("g5: (3)+2 loops",[3],[(0,0),(0,0)]),
]
e1,e2=Fr(sys.argv[1]),Fr(sys.argv[2])
th=Theory(e1,e2,N=12,K=16)
V=Vtx(th)
for name,gen,ed in GR:
    t1=time.time(); worst=10**9; wc=None; n0=0
    for col in product(range(3),repeat=len(gen)):
        val=graph_value(th,V,gen,ed,col).norm()
        if val.v<worst: worst=val.v; wc=col
    print("%s: min val over %d colourings = %d at %s; %.1fs"%(name,3**len(gen),worst,wc,time.time()-t1)); sys.stdout.flush()
