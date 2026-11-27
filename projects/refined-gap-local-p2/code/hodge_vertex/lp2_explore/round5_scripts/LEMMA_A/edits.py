import io
def sub(fn, pairs):
    s = io.open(fn, encoding='utf-8').read()
    for a, b in pairs:
        assert s.count(a) == 1, a
        s = s.replace(a, b)
    io.open(fn, 'w', encoding='utf-8', newline='\n').write(s)
sub('chunk2.md', [
("## 3. Dilaton-point calculus [P]\n\nLet u be",
 "## 3. Dilaton-point calculus [P]\n\nThese are Route A Prop. 2.1, Cor. 2.2 and Lemma 3.1 (lemmaA/A_mv/report.md s.2, s.3.1), re-proved here so that this document stands alone, and extended uniformly to parts equal to 1.\n\nLet u be"),
("**(F5) Bilinear identity.** Pair (F3) with <1|Gamma_+(t) tensor <0|Gamma_+(t~) and use (F4):",
 "**(F5) Bilinear identity.** Note S = [z^0] psi(z) tensor psi*(z). Pair (F3) with <1|Gamma_+(t) tensor <0|Gamma_+(t~), move Gamma_+ past psi(z), psi*(z) by the first line of (F4), and apply the two covector identities to Gamma_+(t) Gop|0> and Gamma_+(t~) Gop|1>:"),
("    U = R_B + B,     n_B := 2g - 2 + |B| + k.                                                     (5.0)\n",
 "    U = R_B + B,     n_B := 2g - 2 + |B| + k.                                                     (5.0)\n\nIn detail: the derivatives d/dy_i, i >= 2, in d^B raise the profile from R_B to R' = R_B + {i in B : i >= 2}. With v = #{i in B : i = 1}, one has (d/dy)^v (J_{R'} + theta)^k G_{R'} = (J_{R'} + v + theta)^k G_{R'+1^v}. For R' empty, v = |B| >= 1, and (2.4) applies to the empty-profile part of F_g, with no constant term needed. The powers of x and sigma are x^{|B|} x^k x^{2g-2} and sigma^{|B|} sigma^{g-1}.\n"),
])
sub('chunk3.md', [
("In degree 3g-3 on w_1+w_2+w_3 = 0 the (3,0,0) terms die by lambda_g^2 = 0; Mumford gives lambda_{g-1}^2 = 2 lambda_{g-2} lambda_g, so Int lambda_{g-1}^3 = 2 A_g; and Sum_{a != b} w_a^2 w_b = -3 e_3. Hence Int prod Lambda^vee(w_a) = (-1)^g A_g e_3(w).",
 "On w_1+w_2+w_3 = 0 the degree 3g-3 part of prod Lambda^vee(w_a) integrates to (-1)^{g-1} [ e_3(w) Int lambda_{g-1}^3 + (Sum_{a != b} w_a^2 w_b) A_g ]: the (3,0,0) terms die by lambda_g^2 = 0. Mumford gives lambda_{g-1}^2 = 2 lambda_{g-2} lambda_g, so Int lambda_{g-1}^3 = 2 A_g, and Sum_{a != b} w_a^2 w_b = -3 e_3(w). Hence Int prod Lambda^vee(w_a) = (-1)^g A_g e_3(w); at w = (1, tau, -1-tau), e_3 = -sigma, this gives the stated value."),
("Draft s.8.2 item 2 of the LP2 draft: its first bullet",
 "Item 2 of LP2 draft s.8.2: its first bullet"),
("| g <= 5 at l <= 2, g <= 4 at l = 3, several tau including zeros of Pi | 595 + 150 pairs, 0 mismatches; 9 kinds of control all fail |",
 "| author: g <= 4, 150 pairs (g, mu), 1387 comparisons; referee 2: g <= 5 at l <= 2, g <= 4 at l = 3, 595 comparisons, several tau including zeros of Pi | 0 mismatches; 9 kinds of control all fail |"),
("| s.8.2 at genus 2-5 against the definition; genus 5-7 zeros, Schur side | Xi_4(7) = 193375/214990848 at tau = 2/7 non-zero, as expected out of range | 0 differences |",
 "| s.8.2 against the definition, h = 2, 3, 4; Schur-side zeros at genus 4-7 | Xi_4(7) = 193375/214990848 at tau = 2/7 non-zero, as expected out of range | 0 differences |"),
("| G <= 5, including tau = -1/3, a zero of Pi_3 | 150 relations, 0 failures; 4 control kinds break it |",
 "| G <= 5, including tau = -1/3, a zero of Pi_3 | 186 relations, 0 failures; 5 kinds of control break it |"),
])
