# Referee report on gap cluster G1 (set-up, s02): part 1 of 2

Author-facing; not part of the note. Task key `gaps/G1/referee`. Object: `NOTE/gaps/G1/GAPS_G1.md` .. `GAPS_G1_5.md`,
`g1_cyclic_check.py` and its log (`NOTE` = `allgenus/note/`), checked against OUTLINE_3 s02, OUTLINE s.0.2, OUTLINE_REVIEW R1, R2, R7,
R8, DECISIONS D2, NOTATION (+`_2`), LABELS (+`_2`), DIGEST_D2_4 (QRR-1..QRR-5), the project sources the prover names, and the literature
it cites, opened in full text by the referee (reading record: part 2, s.10). Part 2 (`REFEREE_2.md`): items 3-5, the design
consequences, the evidence on record. Labels: PROVED, VERIFIED, CONJECTURAL, SKETCH, OPEN. AI-written; no human has read it. No existing
file was modified; the referee's scripts and logs are in this folder.

## 0. Overall verdict: **accepted-with-corrections**

| item | prover | verdict | reason (details in the section named) |
|---|---|---|---|
| 1. (M) | CITED, Givental Cor. 11.7 | **accepted-with-corrections** | The pinned corollary is printed with a misprint in exactly the factor that matters ("Sum_alpha l_alpha", not "Sum_alpha lambda_alpha"). GAPS_G1 s.1.2 corrects it silently. Pin Theorem 11.6, which is printed correctly, and add the two-line derivation. Translation and hypotheses check out; confirmed by the referee's own degree-1 localisation (s.1). |
| 2. `prop:set-cyclicLR` | PROVED | **accepted** | Checked line by line; complete and correct for all `q`. Independent re-implementation to `q^30` and a degree-1 localisation of the equivariant invariant behind `kappa_3`; both have controls. Only notation corrections (s.2). |
| 3(a). passage to `A^QRR_n(q=0) = h_n` | PROVED | **accepted** | Correct, given the convention of `def:set-rows` (part 2, s.3). |
| 3(b). uniqueness at `q = 0` | CITED (+ Step 1 PROVED) | **accepted-with-corrections** | Quotes verified verbatim. Four citation gaps remain (part 2, s.4): (i) the QDE property of the GJR R-matrix, needed in Step 1; (ii) Mumford's relation, which nobody read; (iii) the Teleman/PPZ convention, asserted with "i.e."; (iv) the non-vanishing that gives injectivity, for which the paraphrase of Teleman 8.10 is inaccurate. |
| 4. `lem:set-dilaton` | PROVED | **accepted-with-corrections** | Correct. State the domain hypothesis on `y = q + z1`; use `\varkappa` (R8 (i)); pin `Int psi = 1/24`. The hypothesis `d(0) = 1` must come from a proved statement (part 2, s.5). |
| 5a. symplecticity of `M_0` | PROVED | **accepted** | Correct. VERIFIED by the referee to `t^24`; the controls show that self-adjointness carries the proof (part 2, s.6). |
| 5b. homogeneity of the model CohFT | PROVED | **accepted** | Correct degree count and fifth-root invariance; one naturality sentence to add (part 2, s.7). |

No item is rejected and none needs to be re-opened as OPEN. Two corrections are **required** before the LaTeX is written: the pin of (M)
(s.1.1) and the citation chain of 3(b) (part 2, s.4.2-s.4.5). The others are editorial.

One observation affects the design (part 2, s.4.6). Steps 1-5 of 3(b) never use (QRR). They derive `R(q=0) = diag(h_(w_i))` to all
orders in `z` from pinned literature, so the main session may decide to state (QRR) as a proposition rather than a quoted input.

On the prover's "residual reliance" paragraph:
- For Givental 1996, the issue is not only a terse proof history: the printed corollary is wrong as printed (s.1.1).
- For Teleman, the problem is larger than "injectivity read from 8.7 (ii)". The text also attributes to 8.10 a statement 8.10 does not make (part 2, s.4.5).

## 1. Item 1, the genus-0 mirror theorem (M): accepted-with-corrections

### 1.1 The pinned statement: misprint in Corollary 11.7 (required correction)

The referee read Givental alg-geom/9603021v2 p.23-27 and p.35-42 as page images. For p.40-41, the content streams were also extracted
at font level (`ref_g1_pdfglyph.py`, `.log`). Findings:
- The set-up of s.9 matches GAPS_G1 s.1.2 exactly (p.23-27):
  - `Y_(2,d)`, `E_d`, `S_d(hbar) = E_d/(hbar + c_1^(0))`, `S(t, hbar)`;
  - the `d = 0` convention, `-p` the equivariant Chern class of the Hopf bundle, `phi_i`;
  - `T x T'`, `S'`, `E'_d`, `Z_i`, and `<phi_i, S'> = e^(lambda_i t/hbar)(l lambda_i - lambda') Z_i`.
- Theorem 11.6 (p.40) and `Z^*_i` (p.39) are also as transcribed. Its factor (2) prints `lambda`, cmmi code `0x15`.
- **Corollary 11.7 (p.41) is not as transcribed.**
  - Operation (2) is printed `f(e^t) exp{[g_1(e^t)(Sum_alpha l_alpha) - Sum_a lambda'_a g_(l_a)(e^t)]/(hbar f(e^t))}`. Here `l_alpha` is the Latin letter `l` (cmmi code `l`) with subscript `alpha` (code `0x0b`). It is not `lambda_alpha` (code `0x15`): see `ref_g1_pdfglyph.log`.
  - Operation (1) is printed `T = t + Sum_alpha l_a [..]`. The index mismatch there is harmless.
  - GAPS_G1 s.1.2 transcribes (2) as `(Sum_alpha lambda_alpha)` and says nothing.
- **Why it matters.** Read literally, `Sum_alpha l_alpha` is either undefined (`l` carries the index `a = 1..r`) or equals `Sum_a l_a = n + 1 = 5`. The second reading gives, at the formal quintic point, the factor `C_0 exp(5 g_1/(z C_0)) != C_0`, which makes (M) false. At order `q` the J-function would acquire `+600`; control A1 of `ref_g1_loc.log` refutes this. The note must not cite Cor. 11.7 as printed.
- **The correct form follows from Theorem 11.6 in two lines (PROVED; the referee's derivation).**
  - Put `s^*_i(t) = e^(lambda_i t/hbar) Z^*_i(e^t, hbar)`. This is the definition on p.41; compare `Z^*_i` on p.39. Put `s_i(T) = e^(lambda_i T/hbar) Z_i(e^T, hbar)`, with `Q(e^t) = e^T`.
  - Theorem 11.6 gives `s^*_i(t) = f E_i e^(-lambda_i (T - t)/hbar) s_i(T)`, with `E_i` its factor (2).
  - Using `T - t = Sum_a l_a (g_(l_a) - g_1)/f`, the exponent is `(1/(f hbar)) [g_1 Sum_alpha lambda_alpha - Sum_a lambda'_a g_(l_a) + lambda_i g_1 (Sum_a l_a - (n+1))]`. The last term vanishes in the CY case.
  - So the factor is `f exp{[g_1 Sum_alpha lambda_alpha - Sum_a lambda'_a g_(l_a)]/(hbar f)}`, which is Cor. 11.7 with `lambda_alpha` in place of `l_alpha`.
  - Givental's own proof of Cor. 11.8 (p.42) fits this reading and not the literal one: at `lambda' = 0`, `Sum lambda_alpha = 0` the transformation is "the change ... followed by the multiplication by `f(e^t)`" alone.
- **Correction (required).** Cite "[MR1408320, Theorem 11.6] (Calabi-Yau case `n = 4`, `r = 1`, `l_1 = 5`, `lambda' = 0`)". Print the two-line rewriting for `s_i` above, with the remark that Corollary 11.7 states it with the misprint `Sum_alpha l_alpha` for `Sum_alpha lambda_alpha`. The editor should add the misprint to `design/SOURCE_ERRATA.md`.

### 1.2 Translation (GAPS_G1 s.1.3): checked, correct

The five items were checked line by line.
- Item 1: `f = C_0`; `5(g_5 - g_1) = I~_1`. This is also the mirror map printed in LP18 p.3.
- Item 2: correct once the factor is corrected.
- Item 3: correct.
- Item 4: correct.
  - `Y_(2,d)` is two-pointed (p.23). Forgetting `x_1` gives `pi_*(1/(hbar - psi_0)) = hbar^-1 (hbar - psi_0)^-1`.
  - The eta-dual push-forward is `ev_*(.)/e(O(5))`.
  - `E = ev^* O(5) + E'`, with surjective evaluation because `f^* O(5)` is globally generated, so `E'` is a bundle of rank `5d`.
  - `<phi_i, S'> = S'|_(p_i)`.
- Item 5: correct.

**Independent evidence (VERIFIED, exact; `ref_g1_loc.py`, `.log`, check A).**
- Method: localisation on `Mbar_(0,1)(P^4, 1) = P(S) -> G(2,5)`, with isolated fixed points (line `L_ij`, marking at `p_i`).
- It gives `J_1|_(p_i) = Sum_(j != i) prod_(m=1..5)((5-m) lambda_i + m lambda_j) / (prod_(k not in {i,j})(lambda_j - lambda_k) (z + lambda_i - lambda_j))`.
- This equals `z A_1 - 120 z - 770 lambda_i - 120 e_1(lambda)`, the `q^1` term of Theorem 11.6 rewritten for `s_i`.
  - Here `A_1 = prod_m(5 lambda_i + m z)/prod_j(lambda_i - lambda_j + z)`.
  - Tested at three random rational weight vectors with `e_1 != 0`, and at `lambda_i = zeta^i` in `Q(zeta)`, as rational functions of `z`.
  - By hand: the right side is `z A_1` minus its polynomial part `120 z + 770 lambda_i + 120 e_1`, i.e. its partial-fraction expansion.
- Controls fail as required: (A1) Cor. 11.7 read literally (`+600`); (A2) the `e_1` term dropped, at generic weights; (A3) `E_d` in place of `E'_d`.
- This checks every convention of s.1.3 at order `q`: `lambda' = 0`, fibre weight `5 lambda_i`, tangent weights, `E'_d`, the eta-dual basis, and the signs.

### 1.3 Hypotheses (GAPS_G1 s.1.4): agree

The following are correct:
- the Calabi-Yau case;
- `lambda' = 0` explicitly allowed (p.41-42);
- genus-0 convexity;
- polynomiality in `lambda`, which makes the specialisation `lambda_i = zeta^i` a ring homomorphism;
- the side remark that the higher-genus edge weights do not vanish at the zeta point.

The proof of Thm 11.6 rests on the Main Lemma (p.35), whose proof (end of s.11, p.43-48) neither the prover nor the referee read. That is
acceptable for a CITED statement, but it is one more reason to add an exposition as a second source (s.1.4).

### 1.4 The other candidates (GAPS_G1 s.1.5)

- **Coates-Givental math/0110142v2 Theorem 2** (p.9; read by the referee): the `S^1` acts on the fibres only, and `X` is compact without a torus. Agreed: not a literal pin.
- **"No Corollary 5.1."** The referee searched every page's text: CG numbers its statements globally (Theorems 1-2, Corollaries 1-8), so the claim is literally true. But CG has two statements the prover did not read (reading record p.1-9):
  - Corollary 5 (p.11-12): the Birkhoff-factorisation description of `J_(e,E)` and the mirror map;
  - Corollary 6 (p.12, s.9): the convex-case "mirror theorem" at `lambda = 0`.

  These are probably what PLAN_2026-11-03 meant. Both are non-equivariant on `X`, so the conclusion stands. The sentence in GAPS_G1 s.1.5 should name them.
- **LP18** p.3-4 (definition (5), the specialisation `lambda_i = zeta^i`, footnote 4) and **GJR18** p.18-19 and p.26: read; as the prover says.
- **Expositions** (abstract pages only, NOT read) that a reader may check for a T-equivariant statement before adding them as second sources:
  - Pandharipande, "Rational curves on hypersurfaces [after A. Givental]", math/9806133 (Bourbaki talk);
  - Bini, De Concini, Polito, Procesi, math/9805097;
  - Lian, Liu, Yau, "Mirror principle I", alg-geom/9712011 (its abstract speaks of equivariant Euler classes).
- **arXiv versions:** v1 is a withdrawn 1 KB placeholder; v2 (48 pp.) is the paper (abstract page). The IMRN numbering was not compared, so the prover's `\draftnote` stands.

### 1.5 Verdict

**Accepted-with-corrections.** The printed statement of (M) (GAPS_G1 s.1.1) is unchanged, but its citation changes as in s.1.1.

## 2. Item 2, `prop:set-cyclicLR` (R1): accepted

### 2.1 Line-by-line check

Every step was recomputed; all are correct.
- **(G0-a):** the degree count `dim Mbar_(0,n)(P^4,d) = 5d + n + 1`, `rk E = 5d + 1`, so the degree is `Sum a_j - n`; symmetry through coordinate permutations.
- **Lemma A:** at `zeta`, `e_1 = .. = e_4 = 0`, `e_5 = 1`; the classical term `5 delta`.
- **Lemma B:** Frobenius and commutativity give `kappa_a = kappa_(2-a)`.
- **Lemma C:** correct. One strength deserves a sentence in the note: only derivatives along the line `t = T H` are used. `s_(a+1)` comes from `s_a` through the QDE with `b = 1`, so the J-function is never needed off the line.
- **Lemma D:**
  - `PF5 I_i = 0` recomputed;
  - `lem:set-pf5` from `D q = q (D + 1)`;
  - the constants of `Q[[q]][ell]` are `Q`;
  - the `q^0` part commutes with `D`.
- **Lemma E:** `a_4 = 1 - 3125 q`, `a_3 = -6250 q`, from `prod_(k=1..4)(5D + k) = 625 D^4 + 1250 D^3 + 875 D^2 + 250 D + 24`.
- **Lemma F:** `m_(nj) = delta_(nj) mod q`.
- **Proof, steps 1-5:**
  - the Wronskian identities `W(C_0 f) = C_0^4 W(f)` and `W_D = C_1^6 W_(d_T)`;
  - `det = kappa_1^2 kappa_2`, and the sign of `kappa_1`;
  - the leading coefficient of `Mop`, and `Mop c_0..c_4`;
  - `g = C_0/C_1`.
- **Corollary 2.10:**
  - (a) the eigenvalues, the idempotents, `eta(e_i, e_i) = kappa_1 mu_i^-3`, `Delta_i = C_0^2 zeta^(3i)/L^2` and `D u_i = zeta^i L`;
  - (b) the PF5 closing for every QDE solution, in any differential `Q[[q]]`-algebra, since `Mop = g PF4` is an identity of operators;
  - (c) the ring membership.

### 2.2 Evidence (VERIFIED; never a proof step)

- **The prover's script.** `g1_cyclic_check.py`, re-run from a scratch copy, gives a log identical to the one on record.
- **The referee's own implementation** (`ref_g1_birkhoff.py`, `.log`). It keeps the formal log `ell` explicit and uses a recursion different from the prover's `w`-shift trick. Results, to `q^30` with 15 periods:
  - every `kappa_a` is ell-free and has the stated closed form;
  - the closing `phi^(5)_(k+5) = phi^(0)_k` holds through `k = 9`;
  - the operator identity `d_T^2 kappa_1^-1 d_T^2 C_0^-1 = (C_0/C_1) PF4` holds coefficient by coefficient;
  - `PF4 c_k` is as in Lemma D.
  - Controls fail as required: (c1) non-equivariant closing; (c2) `O(4)`; (c3) the operator identity with `kappa_1 = 1/(C_0^2 C_1^3)`; (c4) `kappa_3 = C_1/C_0`.
- **Degree-1 localisation** (`ref_g1_loc.log`, check B). On `Mbar_(0,2)(P^4,1) = P(S) x_G P(S)` at `lambda_i = zeta^i`:
  - `<H,H>_(0,2,1) = 2875 = 5 [q] kappa_1`;
  - the genuinely equivariant `<H^3,H^4>_(0,2,1) = -3250 = 5 [q] kappa_3`, with `kappa_3 = C_0/C_1 = 1 - 650 q + ..`;
  - the guess `C_1/C_0` fails;
  - the Z/5 grading forces `<H^2,H^3>_(0,2,1) = <H^4,H^4>_(0,2,1) = 0` at zeta, while `<H^2,H^3>_(0,2,1)` is non-zero at generic weights.

  This is independent of (M) and of the Birkhoff recursion, and confirms the equivariant content of the proposition at order `q`.

### 2.3 Corrections (editorial; R8)

The proof's local symbols clash with symbols of the note; each needs a local letter declared at the top of the proof:
- `s_a` (and Givental's `s_i`): `s = u^(1/5)`;
- `g`: the genus;
- `W`: NOTATION already has `Wr` for the Wronskian;
- `w`: bold `w_i`, and the sources' `w` for `zhat` (NOTATION_2 s.9);
- `c_k`: `c_g`, `c_n(L)` and sans `c`;
- `phi_k`, `phi^(a)_k`: CYC's `phi_i` and LMI's `phi` (NOTATION_2 s.9);
- `Mop`: the large-radius analogue of `M_cyc` (`def:gau-conc`); name the relation in a remark;
- `X_a`, `Fhat`, `ctilde_m`;
- `d_T` and `ell`, as the prover already says (GAPS_G1_5 s.6 item 8).

On GAPS_G1_3 s.2.11: the referee read GJR18 (13) on p.26. It is printed `tau-dot = H + q d/dq(I_1/I_0)`. Whether the misplaced `H` is a factor
depends on whether GJR's `I_1` contains `log q`; in either case the intended operator is `C_1 H *`. Calling it a misprint is fair. Keep it
out of the note (corroboration only).

### 2.4 Verdict

**Accepted.** PROVED for all orders in `q`, from (M) and the genus-0 axioms (G0), with no finite check in the proof.
