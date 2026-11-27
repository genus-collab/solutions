import os
import json,sys
d=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "round6c_result.json"),encoding='utf-8'))
o=open('r6c.txt','w',encoding='ascii',errors='replace')
for r in d['referees']:
    print('=====',r.get('key'),r.get('verdict'),file=o)
    print('MAJOR:',r.get('major'),file=o)
    for m in r.get('minor',[]): print(' -',m[:300],file=o)
    for m in r.get('checked',[]): print(' c',m[:350],file=o)
    for m in r.get('not_checked',[]): print(' n',m[:350],file=o)
print('fixer',str(d.get('fixer'))[:1000],file=o)
print('missing',d.get('missing'),file=o)
