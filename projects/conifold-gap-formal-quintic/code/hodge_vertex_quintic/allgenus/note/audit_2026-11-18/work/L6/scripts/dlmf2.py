import sys
import re,sys,html
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
s=open(sys.argv[1],encoding='utf-8',errors='replace').read()
# split at equation ids
parts=re.split(r'(id="E[0-9]+")', s)
cur=None
for p in parts:
    m=re.match(r'id="(E[0-9]+)"',p)
    if m: cur=m.group(1); continue
    if cur in sys.argv[2:]:
        alts=re.findall(r'alttext="([^"]*)"', p)
        disp=[html.unescape(a) for a in alts if 'displaystyle' in a]
        print(cur, ' '.join(disp)[:900])
        # also a few non-display alts near (constraints)
        txt=re.sub(r'<[^>]+>',' ',p)
        txt=html.unescape(re.sub(r'\s+',' ',txt))
        print('   text:', txt[:300])
