import os
import difflib,re,sys
root=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "")
a=open(root+"LP2_REGULARITY_PROOF_draft_round6a.md",encoding='utf-8').read().split('\n')
b=open(root+"LP2_REGULARITY_PROOF_draft.md",encoding='utf-8').read().split('\n')
sm=difflib.SequenceMatcher(None,a,b,autojunk=False)
sec=None
secs={}
cur='pre'
for i,l in enumerate(b):
    if l.startswith('#'): cur=l[:50]
    secs[i]=cur
out={}
for tag,i1,i2,j1,j2 in sm.get_opcodes():
    if tag=='equal': continue
    key=secs.get(j1,'end')
    out.setdefault(key,[]).append((tag,i1+1,i2,j1+1,j2))
with open(sys.argv[1],'w',encoding='utf-8') as f:
    for k,v in out.items():
        print(k,v,file=f)
