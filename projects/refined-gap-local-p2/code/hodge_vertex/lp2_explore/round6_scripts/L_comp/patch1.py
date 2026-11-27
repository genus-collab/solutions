s=open('hodge.py').read()
s=s.replace("""    if g==0 and n==3: return Fr(1)
    if ds and ds[0]==0:""","""    if g==0 and n==3: return Fr(1)
    if g==1 and ds==(1,): return Fr(1,24)
    if ds and ds[0]==0:""")
s=s.replace("""    if chs[-1]>2*g-1: return Fr(0)
    if ds and ds[0]==0:""","""    if chs[-1]>2*g-1: return Fr(0)
    if ds and ds[0]==0 and 2*g-2+n-1>0:""")
s=s.replace("""    if ds and 1 in ds:
        l=list(ds); l.remove(1)
        return (2*g-2+n-1)*chint(g,chs,tuple(l))""","""    if ds and 1 in ds and 2*g-2+n-1>0:
        l=list(ds); l.remove(1)
        return (2*g-2+n-1)*chint(g,chs,tuple(l))""")
open('hodge.py','w').write(s)
