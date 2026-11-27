import re
s=open('gp.py').read()
s=s.replace("def fixed_graphs(g,n,d):","def fixed_graphs(g,n,d,ncol=3):")
s=s.replace("for cols in product(range(3),repeat=V):","for cols in product(range(ncol),repeat=V):")
# value
s=s.replace("""    def value(s,g,ins,d,mode='anc'):
        \"\"\"ins: list of (k,a). Returns the correlator at degree d.\"\"\"
        n=len(ins)
        total=Fr(0)
        for struct,weight in fixed_graphs(g,n,d):
            total+=weight*s.contrib(struct,g,ins,mode)
        return total""","""    ncol=3
    def hval(s,c): return -s.al[c]
    def value(s,g,ins,d,mode='anc',lam=0):
        \"\"\"ins: list of (k,a). Returns the correlator at degree d, with lambda_1^lam of the stabilised curve.\"\"\"
        n=len(ins)
        total=Fr(0)
        for struct,weight in fixed_graphs(g,n,d,s.ncol):
            gs=struct[0]
            pos=[v for v in range(len(gs)) if gs[v]>0]
            if lam==0:
                total+=weight*s.contrib(struct,g,ins,mode,[0]*len(gs))
            else:
                if not pos: continue
                for dist in compositions(lam,len(pos)):
                    rv=[0]*len(gs);mult=math.factorial(lam)
                    for v,r in zip(pos,dist):
                        rv[v]=r;mult//=math.factorial(r)
                    total+=weight*mult*s.contrib(struct,g,ins,mode,rv)
        return total""")
s=s.replace("def contrib(s,struct,g,ins,mode):","def contrib(s,struct,g,ins,mode,rv):")
s=s.replace("val*=(-s.al[cols[legs[m]]])**a","val*=s.hval(cols[legs[m]])**a")
s=s.replace("val*=pref*s.vint(g_v,c,F_S,(),Mk)","val*=pref*s.vint(g_v,c,F_S,(),Mk,rv[v])")
s=s.replace("val*=pref*s.vint(g_v,c,tuple(F_S),tuple(F_X),Mk)","val*=pref*s.vint(g_v,c,tuple(F_S),tuple(F_X),Mk,rv[v])")
s=s.replace("def vint(s,g,col,F_S,F_X,Mk):","def vint(s,g,col,F_S,F_X,Mk,rx=0):")
s=s.replace("key=(g,col,tuple(sorted(F_S)),tuple(sorted(F_X)),tuple(sorted(Mk)))","key=(g,col,tuple(sorted(F_S)),tuple(sorted(F_X)),tuple(sorted(Mk)),rx)")
s=s.replace("rem=dim-base-r\n","rem=dim-base-r-rx\n")
s=s.replace("val=ch1int(g,tuple(sorted(p[0]+p[1] for p in SS)),r)","val=ch1int(g,tuple(sorted(p[0]+p[1] for p in SS)),r+rx)")
s=s.replace("val=anc(g,SS,tuple(sorted(ey)),r)","val=anc(g,SS,tuple(sorted(ey)),r+rx)")
open('gp.py','w').write(s)
s=open('gsum.py').read()
s=s.replace("from intersect import wk\n","from intersect import wk,ch1int\n")
s=s.replace("""    def V(s,i,g,exps):""","""    def V(s,i,g,exps,rv=0):""")
s=s.replace("key=(i,g,exps)\n","key=(i,g,exps,rv)\n")
s=s.replace("ex=dim-sum(exps)\n","ex=dim-sum(exps)-rv\n")
s=s.replace("w=wk(g,tuple(sorted(exps+ds)))","w=ch1int(g,tuple(sorted(exps+ds)),rv)")
s=s.replace("""    def integral(s,g,ins,legsign=-1):""","""    def integral(s,g,ins,legsign=-1,lam=0):""")
s=s.replace("""        for gs,edges,legs,w in graphs(g,n):
            tot=yadd(tot,yscal(s.graph_val(gs,edges,legs,ins,legsign),w))
        return tot""","""        for gs,edges,legs,w in graphs(g,n):
            if lam==0:
                tot=yadd(tot,yscal(s.graph_val(gs,edges,legs,ins,legsign),w))
            else:
                pos=[v for v in range(len(gs)) if gs[v]>0]
                if not pos: continue
                for dist in compositions(lam,len(pos)):
                    rv=[0]*len(gs);mult=math.factorial(lam)
                    for v,r in zip(pos,dist):
                        rv[v]=r;mult//=math.factorial(r)
                    tot=yadd(tot,yscal(s.graph_val(gs,edges,legs,ins,legsign,rv),w*mult))
        return tot""")
s=s.replace("""    def graph_val(s,gs,edges,legs,ins,legsign):""","""    def graph_val(s,gs,edges,legs,ins,legsign,rv=None):""")
s=s.replace("""        dims=[3*gs[i]-3+val[i] for i in range(V)]
""","""        if rv is None: rv=[0]*V
        dims=[3*gs[i]-3+val[i]-rv[i] for i in range(V)]
        if any(x<0 for x in dims): return [Fr(0)]*(s.D+1)
""")
s=s.replace("f=s.V(cols[i],gs[i],tuple(vex[i]))","f=s.V(cols[i],gs[i],tuple(vex[i]),rv[i])")
open('gsum.py','w').write(s)
