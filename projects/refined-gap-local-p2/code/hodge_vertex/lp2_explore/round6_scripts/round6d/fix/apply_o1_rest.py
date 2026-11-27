import os
import io, sys
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "")
def edit(fn, pairs):
    p = D + fn
    s = io.open(p, encoding="utf-8", newline="").read()
    for old, new in pairs:
        if s.count(old) != 1:
            print("COUNT", s.count(old), old[:90]); sys.exit(1)
        s = s.replace(old, new)
    io.open(p, "w", encoding="utf-8", newline="").write(s)
edit("LP2_REGULARITY_PROOF_draft.md", [
 ("+ Lemma 2.4 (written above) + eq:J=I and the genus-0 QDE, through Proposition 2.6.",
  "+ Lemma 2.4 (written above) + eq:J=I (source at generic α: OPEN item O1, s.0.3) and the genus-0 QDE, through Proposition 2.6."),
 ("hence SP1-SP5, eq:J=I, the genus-0 QDE and the constants of FB:A3. The proof of Theorem 2.7(e)",
  "hence SP1-SP5, eq:J=I (source at generic α: OPEN item O1, s.0.3), the genus-0 QDE and the constants of FB:A3. The proof of Theorem 2.7(e)"),
 ("    Prop 2.11 (Givental-Teleman on T_F)              P mod C     = Theorem FB_B(d); step",
  "    Prop 2.11 (Givental-Teleman on T_F)              P mod C     = Theorem FB_B(d), hence OPEN item O1; step"),
])
edit("FRAMEB_PROOF.md", [
 ("| Thm 2.7(e), steps 1-3 | P mod eq:J=I, genus-0 QDE as a whole theorem;",
  "| Thm 2.7(e), steps 1-3 | P mod eq:J=I (source at generic alpha: OPEN item O1), genus-0 QDE as a whole theorem;"),
])
print("ok")
