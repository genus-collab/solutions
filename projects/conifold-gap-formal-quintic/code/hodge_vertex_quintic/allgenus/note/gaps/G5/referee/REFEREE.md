# Referee report on `gaps/G5` (connection problem and main theorem, s08 and s10), part 1 of 2

Author-facing; not part of the note. Key `gaps/G5/referee`. Object, read in full: `NOTE/gaps/G5/GAPS_G5.md`, `_2`, `_3`, `_4`, the scripts
`g5_degree_bound.py`, `g5_flat_hatted.py`, `g5_natfact_checks.py` and their three logs. The three scripts were re-run from copies in the
scratch folder: all three logs are reproduced exactly (up to line endings). Checked against `OUTLINE_7.md` (s08), `OUTLINE_8.md` (s09-s10),
the s02-s03 entries of `OUTLINE_3.md`, `OUTLINE_4.md` that the cluster uses, `NOTATION.md` (+`_2`), `LABELS.md` (+`_2`),
`design/OUTLINE_REVIEW.md` (R1, R6-R8), `design/DECISIONS.md`, the digests DD1 (`_4`, `_5`), DA (`_2`, `_3`, `_4`, `_7`, passages), DL (`_7`),
the project sources listed in part 2 s.7, and the literature (PPZ15, CGYZ26), opened by me. Part 2 (`REFEREE_2.md`): items 3 and 4, the
evidence record, the consolidated corrections, sources and compliance. Labels: PROVED, VERIFIED (exact; range, field or prime stated),
CONJECTURAL, SKETCH, OPEN; finite checks are evidence, never proof. Own code (stdlib Python, `py -I -B`): `ref5_item1.py`, `ref5_item3.py`,
`ref5_item3_main.py`, `ref5_item4.py` in this folder; logs in `logs/`. Nothing outside this folder was modified. AI-written; no human has read it.

## 0. Verdict: ACCEPTED WITH CORRECTIONS

| # | item | prover's resolution | verdict | decisive findings |
|---|---|---|---|---|
| 1 | `thm:con-Pa` step 5 (domain of P(a)) | (b) PROVED | **accepted-with-corrections** | The comparison in `Q((s))` (ii), HPF (iii) and polynomiality (iv) are correct and complete. In (i), identity (1.1) is FALSE as stated (amplitudes with coefficients in `A_L[z]`); it is true for amplitudes with `D`-constant coefficients, which is all that Step 3 uses (F1.1, required). Forward reference s08 -> s09 (F1.2, required, editorial); minor F1.3-F1.4 |
| 2 | `cor:mai-gap` (Corollary B) and the constant of [CGYZ26, Conj. 1.5] | (b) PROVED, with (c) and (a) | **accepted-with-corrections** | Derivation correct. PPZ15 numbering and hypotheses confirmed in the text; REVIEW s.2.3 is a complete argument for `n = 0`; CGYZ26 dictionary and constant confirmed. Corrections: a false sentence on the choice of the primitive of `1` (F2.1), the `|Aut Gamma|` factor counted twice (F2.2), the use of `lem:flh-comm` in `NN` (F2.3), one misdirected pointer into CGYZ26 (F2.4), notation (F2.5) |
| 3 | `cor:mai-natfact` (natural factorisation, `Psi(u)`) | (b) PROVED | **accepted-with-corrections** (minor) | Definition of `Psi(u)` correct and identical to the code's Teleman form; Lemmas 3.A, 3.B, parts (a)-(c), the G8 line and the erratum are correct; VERIFIED independently on actual data with five controls. Minor wording and notation (part 2 s.3) |
| 4 | `prop:con-deg` step 5 (triviality of the resonance) | (b) PROVED | **accepted-with-corrections** | (L0)-(L3) and the induction are correct. The existence and uniqueness clause in `x^(4n+1) Q[[x^-5]]` needs the degree bound `deg RHS_n <= 4n + 5`, which is nowhere stated; I supply it, (L4) below, with a proof (F4.1, required). Notation and one attribution (part 2 s.4) |

No item is rejected and none needs to be reclassified as OPEN. Two corrections are mathematical (F1.1, F4.1). Both are repaired in
this report by a written argument, and in both cases the conclusion of the prover stands. Conditionality (as the prover says, and
correctly): items 2 and 3 use R1 (`prop:set-cyclicLR`, `cor:gau-rowsNN`) and the symplecticity of `M_0`. Both are owned by cluster G1,
which resolves them as PROVED (`gaps/G1/GAPS_G1.md` s.0, items 2 and 5a); that resolution is not refereed here. Item 2 also uses Theorem A,
`thm:gam-D`, (T) and (QRR) as stated in the outline.

## 1. Item 1: `thm:con-Pa` step 5

### 1.1 What the item had to close, and whether the description is right
The outline (`OUTLINE_7.md`, `thm:con-Pa`) asks to name the domain of the flat-hatted comparison. Its source is PB s.2.2, which reads
"on `0 < u < 1`" (`stage5/pb/PB_AUDIT.md` l.104; confirmed). The outline suggests "as rational functions" (`lem:fsp-param` (i)). CP3 s.3.10
(`stage5/chart/CHART_PROOF_3.md`) argues exactly so and invokes Theorem P(a) for the polynomiality of `ph^SP_n` ("By (L1c)(i) both sides are
rational in the parameter; `ph^SP_n` is a polynomial (Theorem P(a))"). Used inside P(a), that would be circular, as G5 s.1.1 says. The
device of a comparison inside `Q((s))` is the one of PS s.4.2 (d) and s.4.5 (`stage5/referees/pb_steps/REFEREE.md` l.209, l.224;
confirmed). I searched HVQ (`allgenus/**`, `natural_form/**`) for another statement of the comparison ring for P(a) and found none.

### 1.2 Part (i), the flat identity in `A_L`: checked line by line
* *Step 1.* `F(L+t; eps) = S(L) + eps L + eps t + Sum_(k>=2) S^(k)(L) t^k/k!`: the gradient at `t = 0` is `eps`, in the ideal `(eps)`, and the
  Hessian `S''(L) = -Delta` is a unit of `A_L`. So the formal IFT (FSL2 s.3.2, `lem:fsp-ift`) applies over `A_L[[eps]]`. The ODE
  `c' = c(c^5 - 1)/5` is correct. `exp(eps D)` is a ring homomorphism, and the formal ODE with initial value `L` has one solution. Hence
  `psi_eps = exp(eps D)` on `A_L`. Correct.
* *Step 2.* FL4 (FSL2 s.4.5) is applied with base ring `A_L`, parameter `eps`, perturbation `delta = eps (L + t)` (in `eps A_L[[t]][[eps]]`, as
  Corollary T' requires) and an amplitude free of `eps`. The bookkeeping `(Delta(eps)/Delta)^(-1/2) (psi_eps(r)/r)^(1/2) = psi_eps(L)/L` uses
  `r/Delta = L^2` and the unique square root with constant term 1; it is correct. So is `d/d eps F(c(eps); eps) = c(eps)`.
* **F1.1 (required).** The identification of the bracket with `e^(..) (psi_eps(L)/L) psi_eps(c^[g])` rests on the sentence "the Taylor data
  at `c(eps)` of ... `g A_St/A_St(c(eps), 0)` are the images under `psi_eps` of the Taylor data at `L`". That sentence is true when the
  coefficients of `g` are rational. It is false for coefficients in `A_L`, which s.1.3 admits: `psi_eps` moves them, while the deformed
  integrand does not. For `g = a(L) g_0` with `g_0` in `Q[sigma, z]` the two sides of (1.1) differ by `z (D a) L c^[g_0]`. *Counterexample:* for
  `g = q`, `(zD + L)(L c^[q]) - L c^[sigma q] = z (Dq) L c = z q L c`, which is nonzero (its `z^1` coefficient is `qL`).
  VERIFIED (`ref5_item1.py` (b), control "g = q(L)": (1.1) fails at every tested point).
  *Correction.* State (1.1) for amplitudes `g` in `Q[sigma, z]`, i.e. with `D`-constant coefficients. This is the hypothesis `D_q a = 0` of RL-C2,
  which Step 3 cites. Step 3 needs nothing more: the amplitudes `sigma^m` and `prod_k (5 sigma + kz)` are in `Q[sigma, z]`, and `q` is applied
  only at the end, by `A_L`-linearity of `SPn`. G5 does this correctly; its "q acting by multiplication on the left" is exactly the
  right order. With the restriction, (1.1) and Step 3 are PROVED.
* *Step 3.* `PF5^L(L c) = L(c^[sigma^5 - 1] - c^[q prod_k (5 sigma + kz)])`: correct, with the factors applied in sequence to amplitudes in `Q[sigma, z]`.
* *Step 4.* I checked directly that the second Mellin-Barnes form has the Taylor data of `SS_0` at orders `>= 2`: its phase differs by a
  term linear in `sigma`, absorbed in the critical-point normalisation. Its normalised amplitude equals that of `def:gam-c`: the half-densities
  differ by the constant `(-1)^(1/2)`, and `SS_2k = bk_k[(5 sigma)^(1-2k) - Sum_j (sigma - zeta^j)^(1-2k)]` in both forms. So both forms have
  the same `c^[g]`. The identity `exp(S'(L+t)) = q (5 sigma)^5/(sigma^5 - 1)` is correct (equal at `t = 0` by `q(5L)^5 = L^5 - 1`, with equal
  logarithmic derivatives). The translate and `lem:fsp-transl` (b) with `p = z` are admissible.
  **F1.3 (minor).** "`z/(sigma - zeta^j)` lie in `z A_L[[t]]`" is true only over `A_L(zeta) := A_L (x) Q(zeta)`, since `(L - zeta^j)^-1` is not in
  `A_L`. The functional equation holds factor by factor there, and the product over `j` is Galois invariant, so `b_z` lies in
  `A_L[[t]][[z]]`. Say so in one clause.
* The holomorphic alternative is correct.

### 1.3 Parts (ii)-(iv): checked line by line, correct and complete
* *(ii) Step 1.* `iota` is injective and `D`-equivariant; I recomputed `D(s^-1) = (1 - s^5)/(5 s^6) = iota((L^6 - L)/5)`.
* *(ii) Step 3.* `d/dX S'(L/X) = 5uX^4/(1 - uX^5)` and `-s^4 WWhat'(X) = (L/X^2) log((1 - uX^5)/(1 - u))`, the latter termwise in `u`: recomputed.
  The two-homomorphism argument is sound. Its domain is `Q[s^(+-1), X^(+-1), (1 - s^5 X^5)^-1]`. Homomorphism 1 goes through
  `Q(s)[[X-1]]` into `Q((s))[[X-1]]`. Homomorphism 2 goes through `Q[X^(+-1)]((s))`, re-expanded at `X = 1` coefficientwise in `s`; this is a
  ring homomorphism because each `s`-coefficient is a finite sum. The two agree on the generators, and on the inverse by uniqueness of
  inverses.
* *(ii) Steps 4-5.* The normalised half-density `(1 - u)^(1/2) (1 - uX^5)^(-1/2)`, the Stirling parts (`sum over j` gives `5 Sum_P binom(2k-2+5P, 5P)
  u^P X^(5P)`), Lemma H with `mu = s^4` (the outline's `mu = L^-4`), and (U1) for `Q[[u]] -> Q((s))`: all correct.
* *(iii).* The chain rule `(zD + L)(s^e F) = s^(e-1) A_e F` is correct: I re-derived `A_e = 1 + zhat(1-u)((4/5) theta_zhat - e/5 - theta_u)`. So is
  `PF5^L(s^-1 Phi) = s^-1 HPF(Phi)` with POLY's normalisation of HPF: `s^-1 HPF = s^-6 HPF_full`, `HPF_full := A_-5..A_-1 - u -
  ((1-u)/3125) prod_j (5 A_-j + j u zhat)`, and `HPF_full` vanishes mod `u` as an operator. Finally, `Q[[u]][[zhat]]` embeds in `Q((s))[[zhat]]`.
* *(iv).* `lem:gau-poly` needs only one power-series solution with `ph_0 = 1` (POLY s.2 Lemma 2; outline `lem:gau-poly`). Its "resonances
  for every choice" clause is not used, so there is no circularity. The order `thm:con-Pa` -> `prop:gam-hatflat` is correct.

### 1.4 Further points
* **F1.2 (required, editorial).** (i) is stated for `c(z; L)` of `def:gam-c`, and its proof cites `prop:gam-rational`. Both are in s09, while
  `thm:con-Pa` is in s08 and s09 cites `thm:con-Pa`. This is not circular: `prop:gam-rational` uses only `lem:fsp-univ` and `lem:fsp-const`.
  But it is a forward reference that G5's list of inputs (s.1.6) does not flag. The editor should either move `def:gam-c` and
  `prop:gam-rational` before `thm:con-Pa` (e.g. into `sec:con-form`), or state (i) for the s08 flat form (the second Mellin-Barnes summand of
  `lem:con-form` (a)-(b), which Step 4 uses anyway) and identify it with `c(z; L)` by `lem:con-form` (d) in s09. The rationality of the
  coefficients of that form is `prop:gam-rational`'s argument verbatim.
* **F1.4 (minor, rule 6).** Several local letters clash with note symbols:
  * `c(eps)` (critical point) clashes with `c(z; L)` and `c_n(L)` in the same proof;
  * `Delta` clashes with `Delta_i`;
  * `g` (amplitude) clashes with the genus;
  * `psi_eps` clashes with `psi`, the psi-classes and the WKB `psi`;
  * `e(eps)` clashes with sans `e`;
  * `iota` clashes with CGYZ26's `iota` in s10 (item 2).
  Rename them and declare them local.
* *Evidence* (`ref5_item1.py`, `logs/ref5_item1_N12.log`; exact over `F_p`, `p = 2^61 - 1`). The series `c^[g]` is computed by an algorithm
  other than G5's: the Morse coordinate `y` with `S(L+t) - S(L) = (S''(L)/2) y^2`, then Lagrange inversion. G5 uses Wick contraction of
  `e^(S_3/z)`; the two agree by `lem:fsp-cov`. Results:
  * (a) `c_n(L_0) = L_0^(4n) ph_n(L_0^-5)` for `n <= 12` at 12 random points, with `ph_n` from my own recursion and the free coefficients
    `[z^5] Ginf = 4537/3937500` and `[z^10] Ginf = -177640111151/31007812500000`.
  * (b) (1.1) for `g = sigma^m` (`m = 0..4`, `n <= 8`, 6 points), using jets `L = L_0 + eps` mod `eps^2`.
  * (c) `PF5^L(L c) = 0` through `z^7` at 3 points, using jets mod `eps^6`.
  * Controls, all failing as they must: `f_m` shifted; `L_0^(4n-1)`; `g = q(L)` in (b); no half-density in (c).
* **Verdict on item 1: accepted-with-corrections** (F1.1, F1.2 required; F1.3, F1.4 minor).

## 2. Item 2: Corollary B (`cor:mai-gap`) and the constant of [CGYZ26, Conj. 1.5]

### 2.1 (a) FOUND: `natural_form/genus4_referee/review/REVIEW.md` s.2.3
I read s.2.1-s.2.4. Steps 0-5 form a complete combinatorial proof, for `n = 0`, of the regrouping `F_g(R' M)` = the outer sum with
vertices `M.omega^`, under its hypotheses H1-H4:
* the edge and leaf splitting;
* orbit-stabiliser on marked graphs;
* the exact sequence of automorphism groups;
* the vanishing of the unstable genus-0 pieces by the degree count with leaves `O(z^2)`.

It is labelled ARGUED there. Its H1 requires `R'(0) = 1`. G5 applies it, through PPZ, to `X = Cf(R~) M^-1` (`X(0) = 1`), and only afterwards
inserts the constant isometry by `lem:flh-frame`. That order is the correct one. The FOUND claim is accurate, and citing PPZ instead of
transcribing REVIEW is the better choice.

### 2.2 (c) CITED: [PPZ15], arXiv:1303.1043 (opened; abstract page: v3 of 2015-03-18 is current; text read in the ar5iv rendering)
* *Numbering, confirmed in the text:*
  * Rem. 2.1: the difference quotient is a power series under the symplectic condition.
  * Def. 2.2: the R-matrix action. `Omega` sits at the vertices, `R^-1(psi)` at the legs, the edge is `(eta^-1 - R^-1(psi') eta^-1 R^-1(psi'')^t)/(psi'
    + psi'')`, and each graph carries `1/|Aut Gamma|`.
  * Prop. 2.3: `R Omega` is a CohFT.
  * Prop. 2.4: the R-matrix action is a left group action.
  * Def. 2.5: translations.
  * Prop. 2.12: the CohFT with unit.
  * Def. 2.13: `R.Omega = R T Omega` with `T(z) = z 1 - z R^-1(z) 1`.
  * Prop. 2.14: the unit-preserving action is a left group action.
  G5's numbers are right.
* *Hypotheses.* `R` must lie in the group of PPZ's (5): `R(0) = 1` and symplectic. This holds for `X` and for `M` (Step 2), and `Cf omega` and
  `M.(Cf omega)` are CohFTs with unit.
* *Coefficient ring.* PPZ work over `Q` (s.0.5: a finite-dimensional `Q`-vector space). The transfer to `Q(zeta)((s))` rests on the remark that
  the proofs are identities of finite polynomial expressions. That remark is G5's own; the note should present it as a remark after
  the citation, not as part of it.

### 2.3 (b) PROVED: checked line by line
* *Definitions (s.2.2).* `Cf(W) = varpi_0 D varpi_1 - varpi_1 D varpi_0 = varpi_0^2 D t_c`: recomputed. So is `D Cf(i_j) = Cf(i_(j+1))` (uses `D Lam' = 1`);
  on `i_2` the relation is the concomitant, with `lem:gau-iso`.
* **F2.1 (correction).** The sentence "any primitive of 1 gives the same value on `Lam`-free elements" is false as worded. `i_0` is `Lam`-free,
  but `Cf(i_0) = varpi_1 - Lam' varpi_0` changes by `-kappa varpi_0` under `Lam' -> Lam' + kappa`. Correct statement: changing the primitive by `kappa`
  composes `Cf` with the automorphism `i_j -> i_j - kappa c_j`, `Lam -> Lam + kappa` of `NN`. That automorphism fixes `C_0`, `I_1`, hence `W`, `C_1`,
  `Delta_i^(1/2)`, the rows and `F_g^NN`; so the values there do not depend on the primitive. Since G5 fixes `Lam' = log(1 - u)`, nothing else changes.
* *Galois descent.* From `X^(i)_a(z) = zeta^(ia) X^(0)_a(zeta^-i z)` (forward relations), `Q^(i)_(a,n) = zeta^(i(a+1-n)) Q^(0)_(a,n)`. With
  `Delta_i^(1/2) = zeta^(4i) C_0/L`, `zeta -> zeta^k` sends colour `i` to colour `ki`, so the colour sum is invariant. G5's sentence is right; print
  the covariance formula.
* `Psi^T eta^-1 Psi = 1` holds in `NN` exactly: `m_0 m_3 = m_1 m_2 = m_4^2 = 1`, with `C_1 = W/c_0^2` a unit. So `R~` lies in `NN[[z]]` and
  `R~(0) = 1`. The division-free edge is a sound device: it makes `F_g^NN` definable without knowing that `R~` is symplectic in `NN`, and it
  equals PPZ's edge for symplectic `X` (Rem. 2.1).
* *Step 1.* `Cf(R~)` is symplectic: from `Cf(Q) = P_can Q^m` (`thm:gau-R0` (iv), `D = 1`), `P_can(-z)^T eta P_can(z) = eta`, and
  `Q^m(-z)^T eta^-1 Q^m(z) = 1`. The last needs `Lambda'^-1 eta^-1 Lambda'^-1 = eta^-1` (`nu_a + nu_(3-a) = 2`) and `M_0` symplectic; write this
  line. Then use `Psi Psi^T = eta`.
* *Step 2.* Correct. The S15 Lemma 1 transcription is checked: with `R' = Psi X`, `R'(-psi)^T eta^-1 R'(-psi') = X(-psi)^T X(-psi')`, and the leaves agree
  because row 0 of `Psi` is `Delta_j^(-1/2)`.
* **F2.2 (correction).** In S15 s.1.1, the source of `def:flh-outer`, `Cont_Gamma` already contains `|Aut Gamma|^-1`. Then `F_g^NN := Sum |Aut Gamma|^-1
  Cont_Gamma` (s.2.2) and `F_g^cf = Sum |Aut Gamma|^-1 T_Gamma` (s.2.4) count it twice. Resolve G5's parenthesis: drop the factor in both places.
* *Steps 3-4 (Theorem A, finite sum).* Correct.
* *Frame remarks (s.2.5).*
  * Complex frames: correct. The `s`-valuation of the images is bounded below independently of `(alpha, beta)`; `c_0^-1` and `W^-1` go to series
    with constant terms `1`, `-1`; the coefficients are polynomial.
  * Scaling: I checked that `sigma` commutes with `D` on every generator. On `c_j`: PF4 is linear with coefficients in `Q[L^(+-5)]`. On `i_j`,
    `j < 2`: Leibniz. On `i_2`: `sigma(i_3) = lambda i_3 + (lambda - mu)(Lam c_3 + 3 c_2)` follows from the bilinearity of the concomitant and the unit
    coefficient of `i_3`. Also `sigma(W) = mu lambda W`, so `sigma` preserves the localisation. Further, `Cf' = Cf o sigma`; the rows scale by
    `diag(mu^-1, lambda^-1, lambda, mu, 1)` (from `C_1 kappa^LR`), an `eta^-1`-isometry; and `Delta^(1/2) -> mu Delta^(1/2)`, `1_omega -> mu^-1 1_omega`.
* **F2.3 (minor).** `lem:flh-comm` is stated for the R-action. In `NN` the symplecticity of `R~` is not known, which is the reason for the
  division-free edge. So use the count of its proof rather than the statement: `sigma` fixes `E_(R~)`, multiplies each leaf by `mu^-1` and
  each vertex by `mu^(2g_v - 2 + n_v)`, and `Sum_v (2g_v - 2 + n_v) = 2g - 2 + #leaves`. The result `sigma(F_g^NN) = mu^(2g-2) F_g^NN` stands.

### 2.4 CGYZ26 (arXiv:2609.13921; abstract page: v1 2026-09-12, v2 2026-10-31 current; v2 read in the arXiv HTML rendering through a text-extraction tool)
* *Confirmed:*
  * Def. 1.1: `K = Q(q)`, `D = q d_q`, and `R` is generated by `I_0^(+-1)`, `(D tau)^(+-1)`.
  * Def. 1.2: `omega_0(0) != 0`; `omega_1` vanishes simply; the log period has the form `omega_1 log Delta + O(1)`; admissibility means that the
    relations are preserved; `iota : R -> R^coni` inside `C((Delta))`.
  * Ex. 1.3: `Delta = 1 - 5^5 q`, PF operator `D^4 + (Delta - 1) prod_j (D + j/5)`, `omega_0 = 1 + O(Delta^3)`, `omega_1 = Delta + (7/10) Delta^2 + O(Delta^3)`.
  * Rem. 1.4: `omega_0 + c_1 omega_1 + c_2 omega_2` with `c_i` in `Q` and `omega_1` fixed.
  * Conj. 1.5: the displayed constant `(-1)^(g-1) B_2g/(5^(g-1) 2g(2g-2) tau_c^(2g-2))` and the space `C[[tau_c]]`.
  * Thm 2: `F_g(Q(q))` lies in `R`, and the gap holds for the quintic.
  * Def. 1.8 and Thm 5 print `(-1)^g` with a constant `kappa`, as G5 item 6 says. Rem. 1.9 names the formal-quintic gap without stating it.
* **F2.4 (minor).** G5 writes "Thm 2 is stated for 'normalized `omega_1`' (Ex. 4.7: `Delta + O(Delta^2)`)". In my reading, Thm 2 says "normalized as in
  Section 1.1"; s.1.1 does not fix the scale; the explicit `omega_1 = Delta + O(Delta^2)` is Prop.-Def. 4.5 (1) in s.4.2; and Ex. 4.7 concerns
  the hypergeometric examples (`omega_p`). Replace the pointer by "Ex. 1.3; s.4.2, Prop.-Def. 4.5 (1)", to be confirmed in the PDF. G5's
  claims rest on Ex. 1.3 and Conj. 1.5 and are unaffected.
* *Dictionary: confirmed.*
  * The power-series solutions are `a f_0 + b f_1 + c f_2`, so `omega_0 = f_0` and `omega_2 = f_2`.
  * `omega_1 = varpi_1`. Its `7/10` I re-derived: the logarithmic solution `varpi_1 log u + O(1)` exists for `7/10` and not for `7/10 + 1/3`
    (`ref5_item3`, log).
  * Ex. 1.3 is the frame `(0, 0)`; Rem. 1.4 is `(alpha, beta) = (c_1, c_2)` in `Q^2`.
  * On the common subring `S` (generated by `c_0^(+-1)`, `(W/c_0^2)^(+-1)` over `Q[L^(+-5)]`), `iota o ev` and `Cf` are differential homomorphisms
    that agree on the generators, hence on `S`.
  * The constant is identical.
* *Differences 1-6:* accurate. I checked `-5/144`, the quintic's constant-map term `chi |B_4 B_2|/(4g(2g-2)(2g-2)!)` at `chi = -200`, `g = 2`. I
  did not check the formal quintic's `-19/9600` (an appendix-C number).
* **F2.5 (minor, R8 (h), rule 6).** Write `Lam'` as `-Xi` (NOTATION table 8). Register `Cf` (and `NN`, with its clash with DP's sans `N`, R8 (h)).
  The scalars `lambda`, `mu` and the endomorphism `sigma` clash with `lambda_i`, the `lambda` of s09, `mu(b, k)` and the Mellin variable `sigma`: rename them.
* **Verdict on item 2: accepted-with-corrections** (F2.1, F2.2 corrections; F2.3-F2.5 minor).
