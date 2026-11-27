# crosscheck_stored.py: decorated-graph comparison of this code with the stored Mathematica results (read only):
#   (a) conifold_CY_graphterms_polar.m: POINT frame, CY torus, genus 2, 3
#   (b) hvc_summary_g{2,3}_Bm.txt: frame B, fibrewise torus, symbolic x
import sys,os,re,time
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from cylib import *
ROOT=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "")
OUT=os.path.dirname(os.path.abspath(__file__))
fs=open(os.path.join(OUT,"crosscheck_out.txt"),"w")
def pr(*a):
    s=" ".join(str(x) for x in a); print(s,flush=True); fs.write(s+"\n"); fs.flush()
def to_py(s,subs):
    s=s.replace("\\\n","").replace("\n"," ")
    for k,v in subs.items(): s=s.replace(k,v)
    s=re.sub(r'(?<![A-Za-z_\d])(\d+)',r'F_(\1)',s)
    s=s.replace("^","**")
    return s
F_=Fr
def split_top(s):
    """split a {a, b, c} list at top level; returns list of element strings"""
    assert s[0]=='{' and s[-1]=='}'
    out=[];d=0;cur=[]
    for ch in s[1:-1]:
        if ch in '{([': d+=1
        if ch in '})]': d-=1
        if ch==',' and d==0:
            out.append(''.join(cur)); cur=[]
        else: cur.append(ch)
    if ''.join(cur).strip(): out.append(''.join(cur))
    return out
def match_brace(s,i):
    d=0
    for j in range(i,len(s)):
        if s[j]=='{': d+=1
        elif s[j]=='}':
            d-=1
            if d==0: return j
    raise ValueError
def nums(s): return [int(x) for x in re.findall(r'\d+',s)]
def parse_nested(s):
    s=s.replace(" ","").replace("\n","").replace("\\","").replace("{","[").replace("}","]")
    return eval(s)
def canon_all(G,psi,dil):
    """all images of a key under Aut(G)"""
    ed=G['ed'];V=G['V']
    out=set()
    for sig,emap,fl in autos_of(G['gv'],G['A']):
        psi2=[[0,0] for _ in ed]
        for e in range(len(ed)):
            e2=emap[e]
            tgt=(sig[ed[e][0]],sig[ed[e][1]])
            ordr=(0,1) if tgt==ed[e2] else (1,0)
            if ed[e2][0]==ed[e2][1] and fl[e]==1: ordr=(ordr[1],ordr[0])
            psi2[e2][ordr[0]]=psi[e][0];psi2[e2][ordr[1]]=psi[e][1]
        dil2=[None]*V
        for v in range(V): dil2[sig[v]]=tuple(sorted(dil[v]))
        out.add((tuple(tuple(x) for x in psi2),tuple(dil2)))
    return out
def part_a():
    txt=open(ROOT+"conifold_CY_graphterms_polar.m").read()
    txt=txt.replace("\\\n","")
    tx=" ".join(txt.split())
    # split genus sections
    i2=tx.index('2 -> {{{1, 1}, {"genera"'); i3=tx.index('3 -> {{{1, 1}, {"genera"')
    for (samp,tup) in (('S1',(3,4,5,11)),('S3',(17,-5,19,-8))):
        W=cyw(*tup)
        a1,a2,e1,e2=[Fr(x) for x in tup]; a0=e1+e2-a1-a2
        subs={"\\[Alpha][0]":"a0","\\[Alpha][1]":"a1","\\[Alpha][2]":"a2","\\[Epsilon][1]":"e1","\\[Epsilon][2]":"e2"}
        env=dict(a0=a0,a1=a1,a2=a2,e1=e1,e2=e2,F_=Fr)
        for g,lo,hi_ in ((2,i2,i3),(3,i3,len(tx))):
            t0=time.time()
            D=Decorated(g,False)
            rows=Rows(W,3*g-3,N=26+6*(g-2))
            fr=Frame(rows,{1:None,2:None},gmax=g+1)
            gs=GraphSum(D,fr,0).run()
            mine={}
            for r in gs.results: mine[(r['graph'],r['key'])]=r['c']
            n=0;bad=0;nomatch=0
            sec=tx[lo:hi_]
            for m in re.finditer(r'\{\{(\d+), (\d+)\}, \{"genera"',sec):
                gi=int(m.group(1))
                st=m.start()
                # header brace
                j=sec.index('{"genera"',st); j2=match_brace(sec,j)
                head=sec[j:j2+1]
                psi=parse_nested(re.search(r'"psiOnHalfEdges" -> (\{.*?\}), "dilatonLeaves"',head).group(1))
                dil=parse_nested(re.search(r'"dilatonLeaves" -> (\{.*\}), "Aut"',head).group(1))
                k=sec.index('{',j2+1); k2=match_brace(sec,k)
                coefs=split_top(sec[k:k2+1])
                vals=[eval(to_py(c,subs),env) for c in coefs]
                G=D.graphs[gi-1]
                imgs=canon_all(G,psi,dil)
                key=None
                for im in imgs:
                    if (gi,im) in mine: key=im;break
                if key is None: nomatch+=1; continue
                c=mine[(gi,key)]
                n+=1
                for idx,v in enumerate(vals):
                    pw=-(2*g-2)+idx
                    if fr2p(v)!=c.get(pw,0): bad+=1
            pr("(a) point frame CY torus, sample",samp,tup,"genus",g,": decorated graphs compared",n,"(enumerated",len(gs.results),"), unmatched keys",nomatch,", coefficient mismatches",bad,"(polar orders t^-%d..t^-1), %.1fs"%(2*g-2,time.time()-t0))
def part_b():
    for g in (2,3):
        path=ROOT+"hodge_vertex_conifold/hvc_summary_g%d_Bm.txt"%g
        lines=open(path).read().split("\n")
        for x in (Fr(3,5),Fr(7,3)):
            W=weights(0,x,1,x,1)
            D=Decorated(g,True)
            rows=Rows(W,3*g-3,N=26+6*(g-2))
            E=(x+1)
            wB=[fr2p(E),fr2p(-x),fr2p(-1)]
            fr=Frame(rows,{1:wB,2:wB},gmax=g+1)
            gs=GraphSum(D,fr,0).run()
            mine={(r['graph'],r['key']):r['c'] for r in gs.results}
            n=0;bad=0;nomatch=0
            pat=re.compile(r'^\s*-- \{(\d+), (\d+)\} psi (\{.*?\}) dil (\{.*?\})\s+polar (.*)$')
            env=dict(x=x,F_=Fr)
            for ln in lines:
                m=pat.match(ln)
                if not m: continue
                gi=int(m.group(1)); psi=parse_nested(m.group(3)); dil=parse_nested(m.group(4)); pol=m.group(5).strip()
                G=D.graphs[gi-1]
                key=None
                for im in canon_all(G,psi,dil):
                    if (gi,im) in mine: key=im;break
                if key is None: nomatch+=1; continue
                c=mine[(gi,key)]
                n+=1
                if pol=="regular":
                    if any(k<0 for k in c): bad+=1
                else:
                    coefs=split_top(pol)
                    for idx,cf in enumerate(coefs):
                        v=eval(to_py(cf,{}),env)
                        pw=-(2*g-2)+idx
                        if fr2p(v)!=c.get(pw,0): bad+=1
            pr("(b) frame B fibrewise, x =",x,"genus",g,": decorated graphs compared",n,"(enumerated",len(gs.results),"), unmatched",nomatch,", mismatches",bad)
part_b()
part_a()
