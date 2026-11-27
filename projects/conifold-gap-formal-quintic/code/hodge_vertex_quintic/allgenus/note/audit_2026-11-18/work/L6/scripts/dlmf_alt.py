import re,sys,html
s=open(sys.argv[1],encoding='utf-8',errors='replace').read()
want=sys.argv[2:] if len(sys.argv)>2 else None
for m in re.finditer(r'<tr[^>]*id="(E[0-9]+)"(.*?)</table>', s, re.S):
    eid=m.group(1); body=m.group(2)
    if want and eid not in want: continue
    alts=re.findall(r'alttext="([^"]*)"', body)
    print(eid, ' || '.join(html.unescape(a) for a in alts)[:700])
