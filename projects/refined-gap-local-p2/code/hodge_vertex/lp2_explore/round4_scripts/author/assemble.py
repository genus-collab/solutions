import io, os, re
D = os.path.dirname(os.path.abspath(__file__))
T = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "LP2_REGULARITY_PROOF_draft.md")
def rd(n): return io.open(os.path.join(D, n), encoding="utf8").read()
s = {n: rd(n + ".md") for n in ["s1", "s2", "s3", "s4", "s5", "s7", "s8"]}
def rep(k, a, b):
    assert a in s[k], (k, a)
    s[k] = s[k].replace(a, b)
rep("s2", "**Theorem 2.1 (R3:Theorem A).**", "**Theorem A (R3:Theorem A).**")
rep("s2", "all of s.2.1 (PF, Theorem 2.1)", "all of s.2.1 (PF, Theorem A)")
rep("s2", "Theorem 2.1 gives", "Theorem A gives")
rep("s2", "- Delta_i^{-1/2}, up to a constant equal to v_i, and the z-only factors h_w(z) and z^{-1}, are accounted for by the identity v S = z^{-1} Sum b_n e_n; h_w and z^{-1} commute with P. \u220e",
    "- Finally, Delta_i^{-1/2} = v_i/sqrt(3 \u03b51 \u03b52 E), so Delta_i^{-1/2} S is a constant times v S = z^{-1} Sum b_n e_n, and the z-only factors h_w(z), z^{-1} commute with P. (R3 said \"Delta_i^{-1/2} commutes with P\"; it does not, as it depends on y. It is absorbed in v S, as here.) \u220e")
rep("s2", "**Lemma 2.1' (R2:Lemma 1; \"Lemma 2.1\" in s.0).**", "**Lemma 2.1 (R2:Lemma 1).**")
rep("s2", "Asym denotes this normalised formal Laplace expansion (BS line ~4245) at a non-degenerate critical point.",
    "Asym denotes this normalised formal Laplace expansion (BS line ~4245) at a non-degenerate critical point. The calculus is written for e^{-\u03c6/z}; BS's e^{\U0001d4b2/z} is \u03c6 = -\U0001d4b2. In dimension one, F5 with z -> -z is how the factors \u0398(-z/a) arise; the difference of Gaussian conventions is the universal constant c of Proposition 2.6.")
rep("s2", "The s_3-integral is F5 with \u03b1 = -E, \u03b2 = -l-1, where \u0398(w, -m) = \u0398(w)/prod_{j<=m}(1 - jw). This gives (s_3^cr)^{-l-1} \u0398(-z/E)/prod_{j<=l+1}(1 + jz/E), v = 1/s_3^cr up to the sign s_3 = -q/V. **The sign (-1)^{l+1} from s_3 = -q/V cancels against (-E)^{-l-1}, so the prefactor carries no l-dependent sign.**",
    "The s_3-integral is computed by the substitution s_3 = -q/V (q the variable of F5), and then F5 with \u03b1 = -E, \u03b2 = -l-1, where \u0398(w, -m) = \u0398(w)/prod_{j<=m}(1 - jw). This gives v^{l+1} \u0398(-z/E)/prod_{j<=l+1}(1 + jz/E), up to a constant. **The sign (-1)^{l+1} from s_3 = -q/V cancels against (-E)^{-l-1}, so the prefactor carries no l-dependent sign.**")
rep("s3", "Let c1 be any invertible section, XX = \u03b8 log c1, and v = (P_0, P_1, P_2) the rows of Definition 1.3 at any colour, in any frame.",
    "Let c1 be any invertible section, XX = \u03b8 log c1, and v = (P_0^B, P_1^B, P_2^B) the unnormalised rows of Definition 1.3 at any colour (the z-only factor h_w does not matter).")
rep("s3", "The P_0-coefficient is -E^3/(27 \u0394_c) \u00d7 \u0394_c. \u220e", "The P_0-coefficient in row 2 is -E^3/27. Self-adjointness of A_0 uses \u03b7^11 = E \u03b7^12 and \u03b7^00 = -(E^3/27) \u03b7^12. \u220e")
rep("s3", "but it produces the pole e1 e2 E^2 (XX(0) + 2/9)/s of Ehat_00 [E]. The proof below does not use closure or the flow. On-shell enters only through Z-infinity, Lemma 2 and H2(d).",
    "but it produces the pole 3 \u03b51 \u03b52 E^2 (XX(0) + 2/9)/t_c of Ehat_00 [E]. The proof below does not use closure or the flow. On-shell enters only through c1 = \u03b8 t_c, i.e. (I1), (I2) of Lemma 4.2, combined with Z-infinity; it is measured by H2(d).")
rep("s5", "(c^(-1) = \u03c3^(-1), S^(-1) = 0, c^(k1) has component structure as in H2)", "(c^(-1) = \u03c3^(-1), S^(-1) = 0, and c^(k1) = 0 for k1 >= 0 by H2)")
rep("s8", "1. **BS genus-0 inputs (C1)-(C3)**, allowed by the user. Through Theorem 2.7 they enter exactly at Proposition 2.10.",
    "1. **BS genus-0 inputs (C1)-(C3)**, allowed by the user. Through Theorem 2.7 they enter exactly at Proposition 2.10. (C3) includes lem:modrel, whose proof in BS is a sketch; see the remark in s.8.3.")
head = rd("s0.md")
out = head.rstrip("\n") + "\n" + "".join(s[k] for k in ["s1", "s2", "s3", "s4", "s5", "s7", "s8"])
io.open(T, "w", encoding="utf8", newline="\n").write(out)
print(len(out.split()), "words")
print(re.findall(r"^#{2,3} .*$", out, flags=re.M))
