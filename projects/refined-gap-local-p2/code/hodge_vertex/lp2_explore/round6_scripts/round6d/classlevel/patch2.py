s=open('intersect.py').read()
s=s.replace("""@lru_cache(maxsize=None)
def anc(g,S,X,r):""","""def _deg(r):
    return r if isinstance(r,int) else sum(c*(2*i+1) for i,c in enumerate(r))

@lru_cache(maxsize=None)
def chint(g,a,rs):
    \"\"\"Int_{Mbar_{g,n}} prod psi^a prod_i ch_{2i+1}^{rs[i]}, a sorted tuple, rs a tuple. Mumford, as in Faber's algorithm.\"\"\"
    n=len(a)
    if 2*g-2+n<=0 or g<0: return Fr(0)
    while rs and rs[-1]==0: rs=rs[:-1]
    if sum(a)+_deg(rs)!=3*g-3+n: return Fr(0)
    if not rs: return wk(g,a)
    if g==0: return Fr(0)
    # peel off the first available ch
    i=next(j for j,c in enumerate(rs) if c>0)
    m=2*i+1;l=i+1
    r2=list(rs);r2[i]-=1;r2=tuple(r2)
    from cl_lib import bern
    s=Fr(0)
    s+=chint(g,tuple(sorted(a+(2*l,))),r2)
    for q in range(n):
        b=list(a);b[q]+=2*l-1
        s-=chint(g,tuple(sorted(b)),r2)
    for u in range(0,2*l-1):
        v=2*l-2-u
        sg=(-1)**v
        s+=Fr(sg,2)*chint(g-1,tuple(sorted(a+(u,v))),r2)
        idx=range(n)
        for mm in range(n+1):
            for I in combinations(idx,mm):
                aI=[a[x] for x in I]+[u]
                aJ=[a[x] for x in idx if x not in I]+[v]
                for g1 in range(0,g+1):
                    g2=g-g1
                    if 2*g1-2+len(aI)<=0 or 2*g2-2+len(aJ)<=0: continue
                    # distribute r2
                    for split in _splits(r2):
                        r1,r3,mult=split
                        x=chint(g1,tuple(sorted(aI)),r1)
                        if x:
                            s+=Fr(sg,2)*mult*x*chint(g2,tuple(sorted(aJ)),r3)
    return s*bern(2*l)/math.factorial(2*l)

def _splits(rs):
    from itertools import product as pr
    out=[]
    for t in pr(*[range(c+1) for c in rs]):
        mult=1
        for c,x in zip(rs,t): mult*=math.comb(c,x)
        out.append((tuple(t),tuple(c-x for c,x in zip(rs,t)),mult))
    return out

@lru_cache(maxsize=None)
def anc(g,S,X,r):""")
s=s.replace("""    if sum(p[0]+p[1] for p in S)+sum(X)+r!=3*g-3+n: return Fr(0)""","""    if sum(p[0]+p[1] for p in S)+sum(X)+_deg(r)!=3*g-3+n: return Fr(0)""")
s=s.replace("""        return ch1int(g,tuple(sorted([p[0] for p in S]+list(X))),r)""","""        ee=tuple(sorted([p[0] for p in S]+list(X)))
        return ch1int(g,ee,r) if isinstance(r,int) else chint(g,ee,r)""")
s=s.replace("""    print('intersect ok')""","""    # chint against ch1int and against lambda classes
    assert chint(2,(),(3,))==Fr(1,2880)
    assert chint(2,(2,),(2,))==2*Fr(7,5760)
    # ch_3 on Mbar_2: 3! ch_3=(-1)^{g-1} lambda_g lambda_{g-1}=-lambda_2 lambda_1 (A3.3), Int =-1/5760
    assert chint(2,(),(0,1))*6==-Fr(1,5760)
    # genus 3: Int lambda_3 lambda_2 lambda_1=1/2 Int lambda_2^3 (FP): 1/2 * 1/ (...) known: Int_{Mbar_3} lambda_1^6=1/90720 (Faber)
    # lambda_1 via ch_1: lambda_1^6 on Mbar_3
    assert chint(3,(),(6,))==Fr(1,90720)
    print('intersect ok')""")
open('intersect.py','w').write(s)
