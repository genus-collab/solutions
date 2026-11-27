# Referee report on `gaps/G5`, part 2 of 2: items 3 and 4, evidence, corrections, sources

Continuation of `REFEREE.md` (object, labels, s.0 verdict table, items 1-2). Same rules.

## 3. Item 3: `cor:mai-natfact` (the natural factorisation and `Psi(u)`)

### 3.1 The definition
`Psi(u) := Qhat(u, zhat = 0)` and `R~ := Psi^-1 Qhat = Psi^T eta^-1 Qhat`. I re-read the docstring of `natural_form/referee/ref_rows.py`. Its
Teleman form `Rt_ji = (1/5) Sum_(a+b = 3 mod 5) Qh_(b,0)^(j) Qh_a^(i)` is exactly `(Psi^T eta^-1 Qhat)_(ji)`. Its forward relations at `zhat = 0`
(`Qh_0 = zeta^i Phi_i/C0`, `kh = (1/C1, C0^2 C1^2, 1/C1, 1/C0)`) give `Qh_a(0) = zeta^(i(a+1)) w_a`. This agrees with Lemma 3.A, with DA_7
item 7, with UR s.4 (Corollary F, where `Psi(u)` appears undefined) and with OU3 s.3, s.6d. The definition is right, and it is the code's.

### 3.2 Lemma 3.A, Lemma 3.B, parts (a)-(c): checked line by line
* *Lemma 3.A.* From `C_1 kappa^LR = (C_1, L^5/(C_0^2 C_1^2), C_1, C_0, C_0)` I recomputed
  `X_a(0) = (L/C_0, zeta^j L^2/(C_0 C_1), zeta^(2j) C_0 C_1/L^2, zeta^(3j) C_0/L, zeta^(4j))`. The rows are `Q^(j) = zeta^j X`. The hatting gives `w`, and `w(0) = eps`.
  The orthonormality sum is correct.
  *Minor.* `Q^(j)_a = zeta^j X_a` is the `Delta^(1/2)` convention (UR s.1, "the normalised rows `Q_a^(i) = zeta^i X_a`"; ref_rows: `Q = P/Delta^(1/2)`;
  NOTATION convention 7). Equivalently, row 0 is `Delta_j^(-1/2)`, the unit. Cite it at that step.
* *Lemma 3.B.* The `s`-exponent, divisibility by 5 (`nu_a = -a mod 5`), the bound `mu >= max(0, (4k-4)/5)` and the zero set
  `{k = 0} U {(1,1)}`: all correct. This is the table of UR s.5. The bound `floor((5n+4)/4)` is correct as an a priori bound (UR's sharper
  `floor((5n+2)/4)` uses (R1), (R2); G5 correctly leaves it to `rmk:mai-degree`).
* *(a).* Correct: `Q^m_b = varsigma^(nu_b - 1) Rhat^mod_b(z/varsigma^4)` and `z/varsigma^4 = tau zhat` (UR s.4; outline `thm:gau-R0` (i)).
* *(b).* Correct. `p_(1,1)(0) = 0` follows from `lem:gau-R1` ("column 1 of `P(0,z)` is `P_11(0,z) e_1`"), and `Ghat(0,0) = 1` from Lemma 3.A. That
  `Qhat` can be evaluated at `u = 0`, i.e. lies in `Q(zeta)[[u]][[zhat]]`, is a consequence of (a); say so in a clause.
* *(c).* The algebra is correct.
  *Minor.* In the symplecticity line, G5 states only that hatting by `Lambda` preserves `eta^-1`. The model rows also need
  `Lambda'^-1 eta^-1 Lambda'^-1 = eta^-1`, where the `varsigma`- and `tau`-powers cancel on `eta`-pairs (`nu_a + nu_(3-a) = 2`). Add that line.
* *Remark (i), the G8 line.* Correct (`A P_0 = P_0 A~`, the ratios `b~_a/b_a` at `u = 0`, and `lim u b_1 = lim u b~_1 = -1`). It agrees with the
  resolution of the owning cluster G2 (`gaps/G2/GAPS_G2_4.md` s.7), so it is not needed twice.
* *Erratum (s.3.4).* Correct. OU3 uses the same index convention (`R' = R^ M_0^-1 = Psi + ..`: `Psi` is flat x colour). By hand,
  `(Psi_0 eta^-1 Psi_0^T)_ab = zeta^(3(a+1)) delta_ab`, which is not `eta^-1`; VERIFIED (control in `ref5_item3`). The correct relation is
  `Psi^T eta^-1 Psi = 1`, equivalently `Psi Psi^T = eta`, and OU3's conclusion ("drops out of every edge") holds with it. Also replace OU3
  s.6d's "orthogonal" by "an isometry from the colour basis (identity) to the lowered flat components (`eta^-1`)".

### 3.3 Evidence of my own (`ref5_item3.py`, `ref5_item3_main.py`, `logs/ref5_item3_U9_K12.log`; mod `2^61 - 1`, `zeta` in `F_p`)
* *Construction.* From scratch, in two frames, `(0, 0)` and `(1/3, -2/7)`:
  * the hatted rows, built from the HPF polynomials of my recursion (free coefficients `[z^5m] Ginf`, i.e. `D = 1` by `thm:gam-D` and
    `cor:gam-free`), the PF4 Frobenius basis (my recursion; the `7/10` re-derived as in part 1 s.2.4), the forward relations and colour covariance;
  * then `Psi`, `R~`, `M_0(tau zhat)` from my own T-system, and `G = R~ M^-1`.
* *Results through `u^8`, `zhat^11`, in both frames:*
  * the rows close (`4 -> 0`), i.e. (CYC) in the frame;
  * Lemma 3.A, `Psi^T eta^-1 Psi = 1`, and (b);
  * `G(u, 0) = G(0, zhat) = 1`;
  * no coefficient of `[u^n] G` beyond `zhat^(floor((5n+4)/4))`;
  * `R~(-zhat)^T R~(zhat) = 1` and `G(-zhat)^T G(zhat) = 1`.
* *Controls, failing as they must:*
  * `M_0` on the left: violations at `(u^1, zhat^3), (u^2, zhat^4), (u^3, zhat^5)`;
  * `M_0(zhat)` without `tau`: the same;
  * `gamma = 7/10 + 1/3`: the rows do not close;
  * OU3's relation: false.
* *Observation for `rmk:con-why` (VERIFIED).* Shifting `f_1` by `1/7`, i.e. `D = 1 + delta z^5 + ..`, leaves (a), (b), the polynomiality of `G` and
  `G(0, zhat) = 1` intact (`G -> G c(z)` with `c` scalar in `z^5 = u^4 zhat^5`). The first failure appears in `G(-zhat)^T G(zhat) = 1`, at
  `(u^8, zhat^10)`, where `c(z) c(-z) != 1` at `z^10`. So of the conclusions of `cor:mai-natfact` only the symplecticity of `G` sees `D = 1`.
  G5's proof uses `D = 1` exactly where it should (`thm:gau-R0` (iv)); the note must not present (a), (b) as tests of `D = 1`.
* *Notation (minor, rule 6).* `w_a` clashes with the leg weights `w_b` (NOTATION_2 table 7) and with bold `w_i`; `m_a` clashes with the order
  index `m`. `Lambda` (leg map) is inherited from NOTATION's own `Ghat` entry and is for the editor.
* **Verdict on item 3: accepted-with-corrections** (all minor: the attribution of `Q = zeta^j X`, the `Lambda'` line, the `u = 0` clause, notation).

## 4. Item 4: `prop:con-deg` step 5 (the triviality of the resonance)

### 4.1 Checked line by line
* *Set-up.* `D = -(1/5)(1 - x^5) theta_x`, `M = x - (z/5)(1 - x^5) theta_x`, `q = (1 - x^-5)/5^5` and `P = M^5 - 1 + (x^-5 - 1) Pi`: correct.
* *(L0).* Correct.
* *(L1).* Correct. `N = (1 - x^5) o (-theta_x/5)` and multiplications commute, so every word containing an `N` has the left factor; the
  commutation identity G5 also writes is true but not needed.
* *(L2)* (ER's lemma) and *(L3)*: correct; I recomputed the `z^1` parts `-(m+2)(1-x^5) x^(m+4)` and `+3 x^(m+4)`.
* *Steps 4-5.* Correct; `x Q[[x]] ∩ Q[x, x^-1] = x Q[x]`.
* *Dictionary.* `rho_n = x^(4n) r_n(x^-5)` with POLY's `r_n` is correct. It follows from `P(s^-1 Phi) = s^-6 HPF_full(Phi)`, `HPF = u^-1 HPF_full` and
  `[w^(n+1)] HPF = (1-u)[(4n - 5 theta_u) ph_n + r_n]`; VERIFIED for `n <= 30` with an independently built hatted operator.
* *Grading, Step 6, the colour-`i` reduction (`x' = zeta^i x`).* Correct.
* *Sources.* ER s.4's sentence ("W stops at `s^-6`") and EVID s.2.1 (b) (the left factor `(1-u)`) are as described. That the two ingredients
  are written nowhere together is fair.

### 4.2 F4.1 (required): the missing degree bound, with proof
The claim "`A_n` ... is the unique solution in `x^(4n+1) Q[[x^-5]]`", and with it the translation in Step 6 into `deg ph_n <= floor(4n/5)`, needs
`deg A_n <= 4n + 1`, i.e. `deg RHS_n <= 4n + 5`. (L0)-(L3) do not give this: the crude count gives `4n + 10`. G5 states it nowhere (its script
tests it).

**(L4) (PROVED).** Every term `z^j x^a theta^b` of `P` has `a <= 4j`.
*Proof.*
1. Give `x` weight 1, `z` weight `-4` and `theta` weight 0 (so `zhat` has weight 0). Then `M = M_1 + M_-4` with `M_1 = x + (z/5) x^5 theta` and
   `M_-4 = -(z/5) theta`, and `5M + kz = 5 M_1 + (weight -4)`.
2. A word in these factors that contains a factor of weight `-4` has weight `<= 4 - 4 = 0`. Hence `M^5` and `Pi` have the same weight-5 part
   `M_1^5`, and all their other parts have weight `<= 0`.
3. `P = M^5 - 1 + x^-5 Pi - Pi`. The weight-5 parts cancel; `x^-5 Pi` and `1` have weight `<= 0`.
4. All weights occurring are `= 0 mod 5`, so `P` has weight `<= 0`. QED
*Consequence.* If every `A_k z^k` (`k < n`) has weight `<= 1`, then `RHS_n z^(n+1)` has weight `<= 1`. So `deg RHS_n <= 4n + 5`,
`deg rho_n <= 4n` and `deg A_n <= 4n + 1`, and the induction closes inside `x^(4n+1) Q[[x^-5]]`. (Equivalently, cite the hatted form of `lem:gau-HPF`.)
VERIFIED as an operator identity: the top exponent is `a = 4j` exactly, for `j = 1..5` (`ref5_item4`).
*Text for s08:* add (L4) to step 3, next to (L1), and the sentence "`deg rho_n <= 4n` by (L4)" to Step 5.

### 4.3 Other points
* **F4.2 (attribution).** G5 says the argument proves EVID's "VERIFIED-only" (b'). ER s.7 item 3 had already relabelled the upper bound and
  "`deg rho_n < 4n/5` for `5 | n`" PROVED, with the terse step 5. Write: "completes ER's proof (ER s.4, s.7 item 3) of EVID s.2.1 (b')".
* **F4.3 (notation, R8 (e), rule 6).**
  * R8 (e) asks for `L`, not `x`.
  * The local `M` clashes with `M = M_0(tau zhat)`, `P` with the gauge, `Pi` with the tower `Pi_i`, `E` with the edges `E_R` and the matrix units
    `E_(ab)`, and `rho_n` with `rho_3` of NOTATION table 3.
  * Rename them and declare them local (the outline's DD1 text has the same clashes).
* *Evidence of my own* (`ref5_item4.py`, `logs/ref5_item4_N40.log`; exact over `Q`).
  * *Independent representation.* `P` is built as a normal-ordered differential operator (`theta^b x^c = x^c (theta + c)^b`), and (L0)-(L4) and the
    grading are checked as operator identities; (L2) is checked on `x^m`, `m <= 80`.
  * *Recursion, `n <= 40`*, for free coefficients 0 and `(-1)^m (3m+1)/11`:
    * the identities of the recursion hold;
    * `rho_n` lies in `x Q[x]` with `deg rho_n <= 4n`, and `[x^0] rho_n = 0`;
    * `A_n` lies in `x Q[x]` with `deg A_n <= 4n + 1`;
    * the hatted `HPF_full`, built from `A_e`, vanishes through `zhat^41`;
    * `ph_5` agrees with POLY s.1.
  * *Controls.*
    * `P + z^2 x^3` keeps (L0) and (L2)-(L4) and the grading, but breaks (L1): `RHS_1` is not divisible by `1 - x^5`, and the hatted resonance
      is obstructed at `n = 5` (coefficient `6/25`).
    * `P + z^2 (1 - x^5) x^-2` keeps (L1) but breaks (L2): the resonance is obstructed at `n = 5` (`1814/3125`).
    * So neither the valuation lemma alone (ER) nor the left factor alone makes the resonance trivial; G5's combination is needed.
* **Verdict on item 4: accepted-with-corrections** (F4.1 required; F4.2, F4.3 minor).

## 5. Evidence record

| script | log | content | result |
|---|---|---|---|
| G5's three scripts (copies in scratch) | compared with `../logs/*.log` | re-run | identical |
| `ref5_item1.py 12 12` | `logs/ref5_item1_N12.log` | item 1 (a)-(c), Morse-coordinate SP mod `2^61-1`, jets | all PASS; 4 controls fail |
| `ref5_item3_main.py 9 12` | `logs/ref5_item3_U9_K12.log` | item 3 on actual data, two frames, `u^8`, `zhat^11` | all PASS; 4 controls fail; `f_1` observation |
| `ref5_item4.py 40` | `logs/ref5_item4_N40.log` | item 4, operator identities, recursion `n <= 40`, hatted HPF | all PASS; 2 controls fail |

Run times are 2-5 s each, single processes. The logs print no paths. Everything here is evidence; no verdict above rests on it, except
the counterexample to (1.1) as stated, which is also a two-line computation.

## 6. Consolidated corrections for the editor
1. *`thm:con-Pa` (i), (1.1):* restrict the amplitude class to `Q[sigma, z]` (RL-C2; F1.1); keep Step 3 as written.
2. *`thm:con-Pa`:* resolve the forward reference to `def:gam-c` and `prop:gam-rational` (F1.2). Add "over `A_L(zeta)`, then Galois invariant" in
   Step 4 (F1.3). Rename the local letters (F1.4).
3. *`cor:mai-gap`, definitions:*
   * replace the sentence on the primitive of `1` by the automorphism statement (F2.1);
   * drop the extra `|Aut Gamma|^-1` in `F_g^NN` and in Step 2, since `Cont_Gamma` contains it (F2.2);
   * print the covariance formula for the Galois descent and the `Lambda'`-isometry line in Step 1.
4. *Scaling remark:* use the count of `lem:flh-comm`'s proof, not its statement (F2.3).
5. *`rmk:mai-constant`:* replace the pointer "Ex. 4.7" by "Ex. 1.3; s.4.2, Prop.-Def. 4.5 (1)", to be checked in the PDF (F2.4). Notation:
   `-Xi`, register `Cf`, rename `lambda`, `mu`, `sigma` (F2.5).
6. *`cor:mai-natfact`:* cite the `Delta^(1/2)` convention for `Q = zeta^j X`; add the `Lambda'^-1 eta^-1 Lambda'^-1 = eta^-1` line and the `u = 0`
   clause; rename `w_a`, `m_a`. Do not present (a), (b) as tests of `D = 1` (s.3.3). The G8 line belongs to G2 (it agrees).
7. *`prop:con-deg`:* add (L4) with the weight proof (F4.1); attribute the completion of ER's step 5 (F4.2); use `L`, rename `M`, `P`, `Pi`,
   `E`, `rho_n` (F4.3).
8. *`design/SOURCE_ERRATA.md`:* G5's three entries (OU3 s.3; PB s.2.2; ER s.4's undefined "W") are correct. Add OU3 s.6d ("orthogonal").

## 7. Overall verdict, sources, compliance
* **Overall: ACCEPTED WITH CORRECTIONS.** The four resolutions are PROVED after the corrections of s.6. Items 1 and 4 need the two
  mathematical repairs, F1.1 and F4.1; they are written out above, and the conclusions stand.
* *Conditionality, exactly:*
  * Corollary B is conditional on Theorem A, (T), (QRR), `thm:gam-D`, R1 and the symplecticity of `M_0`;
  * `cor:mai-natfact` is conditional on `thm:gam-D`, R1, the symplecticity of `M_0`, `thm:gau-R0`, `prop:gau-symp` and `lem:gau-R1`;
  * `thm:con-Pa` and `prop:con-deg` are unconditional, modulo appendix A, `lem:con-form` and the s03 lemmas named.
* *Project sources read* (passages named in the text):
  * `allgenus/stage4/poly/POLY.md` s.0-s.3;
  * `allgenus/stage4/evid/EVID.md` s.2.1 and `.../evid/referee/REFEREE.md` s.4, s.7;
  * `allgenus/stage5/pb/PB_AUDIT.md` s.2.0-s.2.2;
  * `allgenus/stage5/lemmas/FORMAL_SP_LEMMAS_2.md` s.3.2-s.5.4 and `.../lemmas/referee/REFEREE.md` s.3.4, `_2.md` s.7.3;
  * `allgenus/stage5/chart/CHART_PROOF_3.md` s.3.10;
  * `allgenus/stage5/referees/pb_steps/REFEREE.md` s.4;
  * `allgenus/stage2/flat/S15_PROOF.md` s.1-s.2, s.5;
  * `natural_form/theory/UNIFORM_REDUCTION.md` s.1, s.2, s.4, s.5;
  * `natural_form/referee/REFEREE_NAT.md` s.3 (c)-(e) and the `ref_rows.py` docstring;
  * `natural_form/genus4_referee/review/REVIEW.md` s.2;
  * `OU_EXTENSION_genus3.md` s.3, s.6d; `CONIFOLD_SYMBOL_genus2.md` (Psi passage);
  * `gaps/G1/GAPS_G1.md` s.0; `gaps/G2/GAPS_G2_4.md` s.7.
* *Web* (public pages only; no personal data sent): `arxiv.org/abs/1303.1043`, `ar5iv.labs.arxiv.org/html/1303.1043`,
  `arxiv.org/abs/2609.13921`, `arxiv.org/html/2609.13921v2`. The reading went through a text-extraction tool, so the PDF check that G5
  recommends for CGYZ26 remains advisable.
* *Compliance.*
  * New files only, all in `NOTE/gaps/G5/referee/` (this report, four scripts, three logs); G5's scripts were run only as copies in the
    scratch folder.
  * Python 3.14 stdlib, `py -I -B`, at most one process at a time, no cluster, no Mathematica.
  * No topological recursion was used anywhere.
