# Round 6e: rebuild ROUND6_STATUS keeping the existing 'Passes' section verbatim, except the evidence pointer.
import os
P = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "ROUND6_STATUS_2026-11-02.md")
N = os.path.join(os.path.dirname(os.path.abspath(__file__)), "status_new.md")
old = open(P, encoding='utf-8', newline='').read()
new = open(N, encoding='utf-8', newline='').read()
a = old.index('## Passes on the repaired field-of-definition paragraph')
b = old.index('## Open item')
sec = old[a:b].rstrip('\n') + '\n'
ptr = 'part 3 of the class-level test of report ref28num (SP/lp2_round6/round6d/classlevel/)'
assert sec.count(ptr) == 1
sec = sec.replace(ptr, 'part 3 of the class-level test of report ref28num (round6_scripts/round6d/classlevel/)')
open(P, 'w', encoding='utf-8', newline='').write(new + sec)
print('ok', len(new + sec))
