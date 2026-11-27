import os
P = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "ROUND6_STATUS_2026-11-02.md")
s = open(P, encoding='utf-8', newline='').read()
for old, new in [
    ("| 6c | Conifold Prop 13.2 and s.13 | refCon |", "| 6c | Conifold s.13.0, s.13.3 (5), Prop 13.2(a) | refCon |"),
    ("All of it is mod p or at finite order, and it proves nothing.", "All of it is finite-order evidence and proves nothing."),
]:
    assert s.count(old) == 1, old
    s = s.replace(old, new)
open(P, 'w', encoding='utf-8', newline='').write(s)
print('ok')
