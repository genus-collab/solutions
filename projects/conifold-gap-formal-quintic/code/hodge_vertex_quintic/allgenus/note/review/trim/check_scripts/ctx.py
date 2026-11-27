import sys, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from closure import parse, norm
NOTE = sys.argv[1]
O = parse(os.path.join(NOTE, 'archive_full_162pp'))
pairs = [a.split('=') for a in sys.argv[2:]]
for k, r in pairs:
    t = norm(' '.join(O[1][k]['proofs']))
    for m in re.finditer(re.escape(r), t):
        print('#', k, '<-', r, '::', t[max(0, m.start() - 350):m.end() + 200])
        print()
