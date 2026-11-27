import sys, os, difflib, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from closure import parse, norm, refs
NOTE = sys.argv[1]; what = sys.argv[2]; keys = sys.argv[3:]
O = parse(os.path.join(NOTE, 'archive_full_162pp')); N = parse(NOTE)
for k in keys:
    a = O[1][k]; b = N[1].get(k)
    if what == 'stmt':
        ta, tb = a['text'], b['text']
    else:
        ta, tb = '\n=====\n'.join(a['proofs']), '\n=====\n'.join(b['proofs'])
    wa = re.findall(r'\S+', ta); wb = re.findall(r'\S+', tb)
    print('#####', k, what, len(wa), '->', len(wb))
    sm = difflib.SequenceMatcher(None, wa, wb, autojunk=False)
    for op, i1, i2, j1, j2 in sm.get_opcodes():
        if op != 'equal':
            print(' ', op, '| OLD:', ' '.join(wa[i1:i2])[:1500], '| NEW:', ' '.join(wb[j1:j2])[:1500])
