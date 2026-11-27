import re,sys,shutil,os
sys.path.insert(0,'.')
from patch_a import PA
from patch_b import PB
from patch_c import PC
from patch_d import PD
from patch_e import PE
LP=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "")
src=open(LP+"LP2_REGULARITY_PROOF_draft_round4.md",encoding='utf8').read()
lines=src.split('\n')
extra=[]
extra.append((527,527,'(M-a), (M-b). Put N :=',None))
extra.append((786,786,'- Proposition 2.11',r'''- Proposition 2.11 (Givental-Teleman on T_F, mod (C1)-(C3), (C5));'''))
extra.append((788,788,'- Lemma 1.2(c)',r'''- Lemma 1.2(c) (t_c = t_cf, A_cf = C1cf);
- Proposition 2.12 (lem:modrel) and its dictionary, identifying Definition 1.3 at c1 = A with eq:RHb;'''))
P=PA+PB+PC+PD+PE
out=dict()
for s,e,pref,new in extra:
    if new is None:
        l=lines[s-1]; assert l.startswith(pref),l[:60]
        new=l.replace('Put N := ρ P_0j','Put 𝒩 := ρ P_0j').replace('ρ Ehat_0j = N(ψ, w) - N(-w, w)','ρ Ehat_0j = 𝒩(ψ, w) - 𝒩(-w, w)')
        assert '𝒩(ψ' in new
    P=P+[(s,e,pref,new)]
P=sorted(P,key=lambda x:x[0])
# overlap check
for a,b in zip(P,P[1:]):
    assert a[1]<b[0],(a[0],a[1],b[0])
patched=set()
for s,e,pref,new in P:
    assert lines[s-1].startswith(pref),(s,pref,lines[s-1][:60])
    patched.update(range(s,e+1))
# renames on unpatched lines
def spren(l): return re.sub(r'\bF([1-5])\b',r'SP\1',l)
for i in range(223,330):
    if i not in patched:
        lines[i-1]=spren(lines[i-1])
assert 823 not in patched
lines[823-1]=lines[823-1].replace('formal F5','formal SP5')
res=[]
i=1
Pd={p[0]:p for p in P}
while i<=len(lines):
    if i in Pd:
        s,e,pref,new=Pd[i]
        res.append(new); i=e+1
    else:
        res.append(lines[i-1]); i+=1
txt='\n'.join(res)
open('new_draft.md','w',encoding='utf8',newline='').write(txt)
print(len(lines),len(txt.split('\n')),len(txt))
