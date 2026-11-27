import os
import sys
path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "LC_CY_PROOF.md")
d = open(path, 'rb').read()
orig = d
edits = []
def E(old, new, eol):
    edits.append((old.replace('\n', eol).encode(), new.replace('\n', eol).encode()))
LF, CR = '\n', '\r\n'

# header (line 10)
E("they are listed in s.13.3. No topological recursion is used. Nothing here has been refereed.\n",
  "they are listed in s.13.3. No topological recursion is used. Refereed twice (AI; no major item; the record is\n"
  "s.13.1, whose fixes M1-M6 of 2026-11-03 are applied); no human has read this text.\n", LF)

# s.1.2 preamble (M5)
E("This replaces [PR] s.1.2 (C1)-(C14). The code (rc_lib.m rmatrixAtClosed, Proposition 13.2(b) of [PR]) was not re-read; the generic data are\n"
  "taken from [PR] Prop. 13.2(b)(i) and Lemma 13.3:",
  "This replaces [PR] s.1.2 (C1)-(C14). The generic data are taken from [PR] Prop. 13.2(b)(i) and Lemma 13.3 (the code rc_lib.m, rmatrixAtClosed, was read\n"
  "afterwards in the referee pass and re-read for this revision; it agrees with (C1'), (C2'), (C4'), (C5'), (C7') line by line, s.13.1):", LF)

# (C5') remark on LofU (M5)
E(" (C5') Eigenvalue of H*: L(U) := -a1 Q U/(1 - Q U); L_i := L(U_i).\n",
  " (C5') Eigenvalue of H*: L(U) := -a1 Q U/(1 - Q U); L_i := L(U_i).\n"
  "       The code's LofU (rc_lib.m line 81) is -(Q B U + Q a2)/(1 - Q) = (B U + a2)/tau (as 1 - Q = -Q tau), linear in U, not -a1 Q U/(1 - Q U); the difference is\n"
  "       (B U + a2)/tau + a1 Q U/(1 - Q U) = -LeqU/(tau (1 - Q U)) (numerator: (B U + a2)(1 - Q U) + a1 Q U tau = -B Q U^2 + (B - a1) U + Q (a1 - a2) U + a2 = -LeqU,\n"
  "       using a1 Q U tau = a1 U (Q - 1) and B - a1 = a0), and 1 - Q U is a unit of QA (LeqU(1/Q) = (B - a0)/Q - a1 = -a1 tau != 0 for tau != 0). So the two agree at\n"
  "       both roots and are equal in QA: the L_i, and Delta_i = DeltaL(L_i) (line 128), are the same.\n", LF)

# s.1.3 r_D (M3)
E("r_D := the smallest modulus of a root of delta^2 (r_D = oo if a0 = 0 and E^2 + u = 0).",
  "r_D := the smallest modulus of a root of delta^2 (r_D is finite under (ND): for a0 != 0 delta^2 is\n"
  "  quadratic in tau, and at a0 = 0, where E^2 + u = 4 a1 a2 != 0, delta^2 = E^2 - 4 a1 a2 tau).", LF)

# s.1.5 alg (M4)
E("alg^2 = a1 A delta U_1/(E^3 s* (1 - s*)) (Lemma 3.1(c)).",
  "alg^2 = a1 A delta U_1/(E^3 s* (1 - s*)) = A delta (E - B tau s*)/(E^3 s*) (Lemma 3.1(c), (3.4), (3.4')).", LF)

# (3.4') after (3.4) (M4), CRLF region
E("     a unit of O with alg(0) = 1.\n (d) Tensor rule",
  "     a unit of O with alg(0) = 1. Equivalently, by (b)(ii) in its second form (U_1/(1 - s*) = (E - B tau s*)/a1),\n"
  "        alg^2 = A delta (E - B tau s*)/(E^3 s*),   E - B tau s* = delta + tau L_2,                                                              (3.4')\n"
  "     and E - B tau s* has no zero on D: if B tau s* = E then N(s*) = s* (B tau s* - E - a0 tau) + A = -a0 E/B + A = a1 a2/B != 0. Since delta and s* are units\n"
  "     of O(D) as well (s.1.3), alg^2 is a unit of O(D) and 1/alg is analytic on D.\n (d) Tensor rule", CR)
E("with value 1 (two square roots of a unit of O with the same value at 0 coincide).\n(d) Verbatim [PR]:",
  "with value 1 (two square roots of a unit of O with the same value at 0 coincide). (3.4') is (3.4) with U_1/(1 - s*) = (E - B tau s*)/a1, the second form of (b)(ii);\n"
  "the value of N at s = E/(B tau) is direct.\n(d) Verbatim [PR]:", CR)

# s.3.9 uniqueness reason (M6), CRLF region
E("geometric series E/(E + z). The same scalar h_nw multiplies both idempotents",
  "geometric series E/(E + z) (h_v depends on v only through the odd power sums Sum_a v_a^(1-2k), k >= 1, of the reciprocals, and h_v = 1 iff all of them vanish\n"
  "iff the multiset {1/v_a} equals its negative, by Newton's identities (the power sums of a finite multiset determine it, and the even ones of {1/v_a} and {-1/v_a}\n"
  "agree); so h_{w'} = h_w iff h_{w' - w} = 1 iff w' and w agree up to pairs v, -v). The same scalar h_nw multiplies both idempotents", CR)

# Prop 3.19 pointer to (3.4') (M4), CRLF region
E("i.e. to show R_tr = delta (log alg)_tau; since both are analytic on D with value 1 at 0 and Phi_0 is the unique solution",
  "i.e. to show R_tr = delta (log alg)_tau; since both Phi_0 and 1/alg are analytic on D with value 1 at 0 (1/alg by (3.4')) and Phi_0 is the unique solution", CR)

# Prop 3.19 display (M2), CRLF region
E("    -s* (delta + tau L_2)/U_1 - L_2 + L_2 s*/(1 - s*) - delta_tau + 2 B tau s* L_2/delta = 2 B tau s* L_2/delta,\n"
  "i.e. (with delta_tau = a0 - 2 B s* - 2 B tau s* L_2/delta, whose last term cancels the left side's last term)\n",
  "    -s* (delta + tau L_2)/U_1 - L_2 + L_2 s*/(1 - s*) - delta_tau = 2 B tau s* L_2/delta\n"
  "(the right side is 2 R_tr - 2 delta_tau = 2 B tau s* L_2/delta), i.e. (inserting delta_tau = a0 - 2 B s* - 2 B tau s* L_2/delta, whose last term cancels the right side)\n", CR)

# Prop 6.1 k = 0 (M1), LF region
E("k = 0, by hand. From s.3.11, G_1 = (E/(E + z))(H_1 - (E^2 + u)/(2 E^2)) = [2 E z",
  "k = 0, by hand. From s.3.11 and (4.3) (tau/t = 1 - t/2 + O(t^2)), G_1 = (E/(E + z))(H_1 - 1/2) (= [t^1] r_{0,1} + (u/(2E^2)) G_0, since G = rho r_{0,1} and\n"
  "rho_1 = u/(2E^2)) = [2 E z", LF)

# s.13.1 code record (M5)
E(" - s.1 (C1')-(C14'): definitions. (C1') is taken from [PR] Prop. 13.2(b)(i) and Lemma 13.3, not from rc_lib.m [P mod: that the code's generic superpotential is (C1')].\n",
  " - s.1 (C1')-(C14'): definitions. (C1') is taken from [PR] Prop. 13.2(b)(i) and Lemma 13.3; that it is the code's generic superpotential is now confirmed by reading\n"
  "   rc_lib.m (hodge_vertex_conifold/rc_lib.m, read by the hand referee and re-read for this revision): line 80 spotU = a1 Log[1 - Q U] - a2 Log[U] + (a0 + a2) Log[1 - U]\n"
  "   is (C1') verbatim; line 115 LeqU = Q (a0 + a1) U^2 - (a0 + Q (a1 - a2)) U - a2 is (C2'); line 82 etaSym is (C4'); lines 106-109 ffcClosed =\n"
  "   Exp[-Sum_k (B_2k/(2k (2k - 1))) z^(2k-1) (a1^(1-2k) + (a0 + a2)^(1-2k) - e1^(1-2k) - e2^(1-2k))] is h_nw, nw = (a1, A, -e1, -e2), in the sign convention of (C7');\n"
  "   line 138 rmatrixAtClosed passes the same list ff to both roots; line 81 LofU = -(Q (a0 + a1) U + Q a2)/(1 - Q) equals L(U) of (C5') in QA (remark under (C5')),\n"
  "   and line 128 evaluates Delta at L = LofU [P].\n", LF)
E("Theorem LC_CY: P modulo [PR] Proposition 3.15 and [PR] Lemma 13.3 (both Laplace + Stirling), i.e.\n"
  "   the same analytic inputs as [PR] Theorem LC, plus the generic superpotential (C1').\n",
  "Theorem LC_CY: P modulo [PR] Proposition 3.15 and [PR] Lemma 13.3 (both Laplace + Stirling), i.e.\n"
  "   the same analytic inputs as [PR] Theorem LC (the generic superpotential (C1') is the code's, see the first item).\n", LF)

# referee record at the end of s.13.1
E("h_nw = 1) are replaced by Lemma 3.1(b)(ii) in its two forms, the degree/leading-coefficient argument of Lemma 5.1_CY, a1 a2 - E L_2 = -B tau s* L_2, and the constant line (3.9a).\n"
  "\n### 13.2 What was not re-derived\n",
  "h_nw = 1) are replaced by Lemma 3.1(b)(ii) in its two forms, the degree/leading-coefficient argument of Lemma 5.1_CY, a1 a2 - E L_2 = -B tau s* L_2, and the constant line (3.9a).\n"
  "\n"
  "Referee record (2026-11-03). Two referee passes, both AI, both with no major item. (1) A hand referee re-derived s.1, s.3-s.10, Corollary 12.1 and Remark 0.1\n"
  "line by line against [PR]; verdict: sound with minor fixes. Its items M1-M6 are applied in this revision: M1 the first expression for G_1 in the k = 0 computation\n"
  "of Proposition 6.1_CY (s.6.2; the printed intermediate was [t^1] r_{0,1} = -u z/(E (E + z)(E + 2z)), the numerator and the value were right); M2 the sign of\n"
  "2 B tau s* L_2/delta in the display of the proof of Proposition 3.19_CY (s.3.10; the final line was right); M3 the parenthesis on r_D in s.1.3 (the case a0 = 0,\n"
  "E^2 + u = 0 is excluded by (ND)); M4 the form (3.4') of alg^2, which makes the analyticity of 1/alg on D in Proposition 3.19_CY manifest; M5 the code record above\n"
  "(rc_lib.m read) and the remark on LofU under (C5'); M6 the one-line reason for the uniqueness of the class w in s.3.9. Every changed line was re-derived by hand\n"
  "before the change. (2) An independent numerical referee with its own F_p code (p = 2^127 - 1 and 2^89 - 1; truncated Laurent series in t to t^14, z^23) tested the\n"
  "KEY FORMULA, F1-F3, L1 (k <= 13), L2 (k <= 9), L3 (k <= 9), the flow (7.2), Theorem Q_CY, the regular-colour values (s.8), Corollary 12.1 at a0 = 0, and Remark 0.1\n"
  "at non-CY points, at nine weight points including a1 = a2, u = 0 and e1 = e2; all passed, and four controls failed as required: the class (E + 1, -e1, -e2) in place\n"
  "of w; no h_nw; A -> a2 in the rows; lam_j = 1/(E + j psi). It also observed that at u = 0 the top leaf lam_{k+1} is absent from [t^k](rho q_t) for odd k, and that at\n"
  "E = 9 a2/7 the lam_9 coefficient of [t^8](rho q_t) vanishes, as the sharpness statement predicts (for odd k and E = 2 a2, resp. E = 2 a1, P_k contains the factor\n"
  "(k + 1) a2 - ((k + 1)/2) E, resp. (k + 1) A - ((k + 1)/2) E; P_8 contains 9 a2 - 7 E). Both passes are evidence only, not proof steps. No human has read this text.\n"
  "\n### 13.2 What was not re-derived\n", LF)

# s.13.2 code item (M5)
E(" - [PR] Lemma 8.1 (conjugate encoding in the code): cited; the claim that rc_lib.m's rmatrixAtClosed uses the same scalar h_nw at both roots is [PR] Prop. 13.2(b)\n"
  "   and Lemma 13.3, not re-read in the code.\n",
  " - [PR] Lemma 8.1 (conjugate encoding in the code): cited. That rc_lib.m's rmatrixAtClosed uses the same scalar at both roots is now read in the code (line 138,\n"
  "   s.13.1); that this scalar is the Q = 0 constant h_nw is [PR] Prop. 13.2(b) and Lemma 13.3, not re-proved.\n"
  " - The two referee passes of s.13.1: their derivations and F_p code were not re-run here; only the fixes M1-M6 were re-derived.\n", LF)

for i,(o,n) in enumerate(edits):
    c = d.count(o)
    if c != 1:
        print("EDIT", i, "match count", c, "for", o[:80]); sys.exit(1)
    d = d.replace(o, n)
open(path, 'wb').write(d)
print("applied", len(edits), "edits; bytes", len(orig), "->", len(d))
print("CRLF", d.count(b'\r\n'), "bare LF", d.count(b'\n')-d.count(b'\r\n'), "bare CR", d.count(b'\r')-d.count(b'\r\n'))
