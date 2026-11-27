# Referee report on gap cluster G1 (set-up, s02): part 2 of 2

Continuation of `REFEREE.md` (overall verdict, items 1-2). This part covers:
- s.3 item 3(a); s.4 item 3(b); s.5 item 4; s.6 item 5a; s.7 item 5b;
- s.8 the prover's design consequences;
- s.9 the evidence on record; s.10 the referee's reading record.

## 3. Item 3(a), passage to `A^QRR_n(q = 0) = h_n`: accepted

The argument is correct.
- `eta(1, e~_i) = Delta_i^(-1/2) = zeta^i L/C_0` under NOTATION_2 s.10 item 7.
- At `q = 0` we have `u = 1`, `L = C_0 = 1`, `zhat = z`.
- Hence `Phi_i(z, 1) = h_(w_i)(z)`, and `A_n = L^(4n+1) ph_n(L^-5) = [z^n](C_0 Q^(0)_0)`.
- `ph_0 = 1` follows from `R = 1 + O(z)`.
- The diagonal entries of `R(0)` do not depend on the sign of `e~_i`. The value of `Phi_i` does, and it uses the sign fixed in NOTATION_2.

Two remarks for the writer:
- **The convention is an input.** The proof uses (QRR) exactly as the design states it: "Givental's convention" means BOTH legs `R^-1(psi)` and that the columns `R e~_i e^(u_i/z)` solve `z d_T v = H * v` (`def:set-rows`; DIGEST_D2_4 QRR-1). This is a convention built into (QRR), not something G1 proves. If (QRR) is derived as in s.4.6, the convention identity must be cited (s.4.6).
- **Adjacent step (not a G1 item).** `def:set-rows` also states the covariance `Phi_i(zhat, u) = Phi_0(zeta^-i zhat, u)`, which with the new definition of `Phi_i` needs a proof. Route:
  - PF5 is invariant under `z -> zeta z`, since `z^5` is;
  - `u_i = zeta^i u_0`, from `D u_i = zeta^i L` with the constants of `lem:set-canon`;
  - formal WKB solutions with a given exponential are unique up to a factor `c(z)` (`lem:gau-wkb`);
  - `c` is fixed at `q = 0` by (QRR) and `lem:set-hrot`.

  The s02/s03 writer should place it.

## 4. Item 3(b), uniqueness at `q = 0`: accepted-with-corrections

### 4.1 Quotes: verified verbatim (referee's reading)

- **Teleman 0712.0160v3:**
  - Theorem (2) (p.11; the paragraph number "(1.15)" was not seen by the referee);
  - 8.2 and 8.4 Proposition (p.39);
  - 8.5 Corollary, (8.6) with `Zbar[h] := (nth co-power of) alpha^g exp{Sum h_(2j-1) ((2j)!/B_(2j)) mu_(2j-1)}`, and 8.7 Proposition (i)-(ii) (p.40);
  - the Riemann-Roch identity `mu_(2j-1) = (B_(2j)/(2j)!) kappa_(2j-1)` over `M_g` (p.40);
  - 8.10 Proposition (p.41).
- **PPZ15 1303.1043v3:**
  - Def. 0.3 with axiom (iii) (p.9);
  - (5) and the symplectic condition (p.19);
  - Def. 2.2 and the convention sentence (p.20);
  - Prop. 2.12, Def. 2.13, Prop. 2.14 (p.27-30).
- **GJR18:** s.4.1 and Thm 4.1 (p.18-19). **LP18:** eq. (26) (p.23).
- **Numbering:** v3 has "8.2 Proposition" but no section heading "8.2". The ar5iv claim was not checked.

### 4.2 Step 1 (PROVED): correct

The `q^0` argument is right, and robust: if the R-matrix solved the QDE only after `z -> -z`, the conclusion (diagonal, `exp(odd)`)
would not change.

**Correction C3b-1 (required).** Step 1 is applied to the R-matrix of GJR18 Thm 4.1, but GJR18 s.4.1 states only that this R-matrix is
symplectic. The note must cite that it solves the QDE (the flatness system):
- Teleman 8.7 (ii): the theories are `T_z^-1 E(z) exp h(z) T_z (I_A)`, with `E` a solution of the system of 8.2 (ii) (cf. 8.4);
- `exp h` is constant and diagonal, so `E exp h` solves the same system.

### 4.3 Step 2: correct

Recomputed:
- `H_g = W_i^(g-1) prod_(a in w_i) Lv(a)/a^g`;
- `Lv(a)/a^g = c_(-1/a)(E) = exp(-Sum_k (2k-2)! ch_(2k-1)(E) a^(1-2k))`;
- `W_i = Delta_i(0) = zeta^(3i)`, since `prod_(j != i)(zeta^i - zeta^j) = 5 zeta^(4i)`.

Two additions:
- **Correction C3b-2 (required).** Mumford's relation `c(E) c(E^vee) = 1` (with `ch_(2k)(E) = 0`) is used without a source read in full; the prover's reading record lists Faber-Pandharipande as not read. Pin it (Mumford 1983 s.5, or [MR1728879] s.2.1) after reading.
- **One sentence to add.** "The `q^0` part is the degree-0 part" also needs the CohFT-class version of an argument the prover gives only for invariants (Lemma A). In degree 0, insertions of `T H` give classes pulled back along forgetful maps, so their push-forwards vanish; and the divisor equation holds for the classes.

### 4.4 Step 3: correct, but the convention needs a pin

**Correction C3b-3 (required).** The identification of Teleman's `T_z^-1 exp h(z) T_z (I_A)` with `exp(h).omega` in the convention of
PPZ15 Def. 2.13 is asserted with "i.e.". Justify it by both of the following:
- **PPZ15 p.20, after Def. 2.2.** Their choice (inverse R-matrix on the legs, a left action) is, quoting PPZ, "the same convention is used by Givental and Teleman in their papers".
- **The referee's genus-1 check, which needs nothing from Teleman.** In PPZ's convention, take the rank-one theory `exp(Sum_k a_k ch_(2k-1)(E))` and `R = exp(Sum r_(2k-1) z^(2k-1))`.
  - Then `Int_(Mbar_(1,1)) (R.omega)_(1,1)(1) = r_1 (Int kappa_1 - Int psi_1) + r_1/2 = r_1/2`. The second term comes from the self-loop graph: `|Aut| = 2`, edge term `r_1`.
  - On the Hodge side `Int_(Mbar_(1,1)) a_1 lambda_1 = a_1/24`.
  - So `r_1 = a_1/12 = a_1 B_2/2!`, as `h_(2k-1) = a_k B_(2k)/(2k)!` requires; a mismatch of conventions would flip the sign of every odd coefficient, so this order suffices.
  - For colour 0, `a_1 = -Sum_(a in w_0) a^-1 = -9/5`, so `r_1 = -3/20 = h_1`. This is the value that HPF forces on the rows (`rmk:set-qrrsign`, `lem:gau-HPF`).

The prover's own sign cross-check rests on LIT s.1.3 (Pandharipande 1712.02528, not read). Replace it by the two items above.

### 4.5 Steps 4 and 5

**Step 4** is correct:
- for fixed `(g, n)` and degree, the graph sum is polynomial in coefficients that are power series in `q`;
- a diagonal R-matrix acting on a direct sum of rank-one TFTs keeps the colours (the translation `T` decomposes as well).

**Step 5: correction C3b-4 (required).**
- **The paraphrase of 8.10 is inaccurate.** The text says "8.10 Proposition (rank one, p.41): the classification datum is an invariant of the theory". 8.10 says that formula (8.9) describes all rank-one CohFTs, and that flat-vacuum theories are those with `a_j = 0` for `j > 0`. It does not say that the parameters `b_(2j-1)` of (8.9) are determined by the theory.
- **What injectivity needs.** The injectivity actually needed is: for every `k >= 1`, `ch_(2k-1)(E) != 0` in `H^*(Mbar_(g,n))` for some `(g, n)`. Then `exp(Sum_k d_k ((2k)!/B_(2k)) ch_(2k-1)(E)) = 1` forces `d_k = 0`: take logarithms (the class is unipotent), separate degrees, and use `B_(2k) != 0`.
- **The prover's route is sound but mis-attributed.** It uses the Riemann-Roch identity over `M_g`, which Teleman states on p.40, and `kappa_(2k-1) != 0` in `H^*(M_g)` for large `g`. The non-vanishing of the odd `kappa` classes in the stable range is due to Miller (1986) and Morita (1987); Madsen-Weiss prove more. None of these was read.
- **A more elementary route.** `ch_(2g-1)(E) = (-1)^(g-1) lambda_(g-1) lambda_g/(2g-1)!` on `Mbar_g`, from Newton's identities and Mumford's relation; the referee checked it by hand for `g = 1, 2, 3`. Add a pinned non-zero evaluation, e.g. Faber's `Int_(Mbar_g) lambda_(g-2) lambda_(g-1) lambda_g`, or Getzler-Pandharipande's `Int_(Mbar_(g,1)) psi_1^(g-1) lambda_(g-1) lambda_g` (not read).
- **Draftnote.** Keep ONE `\draftnote` naming exactly this non-vanishing input until a statement is read in full.

### 4.6 Observation for the main session: (QRR) follows

Steps 1-5 do not use (QRR). They derive `R(q = 0) = diag(h_(w_i))` to all orders in `z` (the `z^(5k)` parts included) from:
- GJR18 Thm 4.1, with Teleman 8.7 (ii) for the QDE property;
- LP18 (26);
- Mumford's relation;
- Teleman 8.7 (i), with the PPZ convention (s.4.4);
- the non-vanishing of s.4.5.

The quoted sources of (QRR), GJR18 (17) and GJR17 (19)-(20), only refer to Givental 2001 and Coates-Givental (DIGEST_D2_4 QRR-1). So
once the open citations of s.4.3-s.4.5 are pinned, the main session may state (QRR) as a proposition with pinned citations, and remove it
from the conditionality of Theorem A.

One caveat: (QRR) also packages the convention of `def:set-rows` (s.3). In that case it must be cited separately: the CohFT R-matrix
with legs `R^-1(psi)` is the one whose columns `R e~_i e^(+u_i/z)` solve the QDE, for instance through GJR18 (14) and the Birkhoff
factorisation of Lemma 5.2. The sign of `z` is pinned by `h_1 = ph_1 = -3/20` (s.4.4). This is a design decision, not a G1 gap.

### 4.7 Verdict

**Accepted-with-corrections** (C3b-1 to C3b-4).

## 5. Item 4, `lem:set-dilaton`: accepted-with-corrections

**(i) Correct.**
- PPZ15 Def. 0.3 (iii) (p.9); Arbarello-Cornalba (1.7) (p.6) and "kappa_0 = 2g - 2 + n" (p.5), read by the referee.
- AC's `kappa_i = pi_*(K^(i+1))`, `K = c_1(omega(Sum D_i))`, is PPZ's convention (p.4, footnote 2).
- (1.7) already has the `psi_j` of `Mbar_(g,n+1)` on the left, which handles `psi_j psi_(n+1) = (p^* psi_j) psi_(n+1)`.

**(ii) Correct.** For the model TFT, `omega_(1,1)(1) = Sum_i Deltahat_i^(1/2) Deltahat_i^(-1/2) = 5`. Pin `Int_(Mbar_(1,1)) psi_1 = 1/24`
(standard, e.g. Witten 1991; not read), as the prover itself flags.

**(iii) Correct, recomputed.**
- the rising factorial `(N)_b` and the binomial series;
- the `h = 0` unstable terms vanish (`psi = 0` on `Mbar_(0,3)`);
- `Fbar_1(kappa q) = Fbar_1(q) - (5/24) log kappa`, as in OUTLINE_3.

**Corrections.**
- **C4-1 (required in the statement).** Make the domain a hypothesis. `Fbar_h` is a formal function of `y = q + z 1`, and the expansion at `q/kappa` needs `y` formal (or in `m V[[z]]`) together with `kappa in 1 + m`. The prover says this in the proof; it belongs in the statement.
- **C4-2.** Notation:
  - `kappa` becomes `\varkappa` (R8 (i)); it clashes with `kappa_a`;
  - `Lam` needs a local letter, since `Lambda` is the leg map of s06.
- **C4-3.** The hypothesis `d(0) = 1`.
  - Current support: the prover's source is the bracketed "consistent with" remark in S15 (stage2/flat/S15_PROOF.md, line 90); the identity `d = tau^(1/4) varpi_0` is there only VERIFIED (theory/referee/rf_flat.log).
  - What the writer should do: write the one-line derivation. S15 (I4), which is proved there, gives `omega^Q = d.omega^m`, so `d` is the ratio of the unit coefficients of the two TFTs. Feed in the two normalisations: `xhat = varpi_0` for the quintic (`def:set-hat`) and `varsigma = s tau^(-1/4)` for the model (`def:set-tau`). This should give `d = tau^(1/4) varpi_0` and hence `d(0) = 1`; the referee did not write out this line.
  - Status: until that line is written, `d(0) = 1` is an input of the chain, not of the lemma.
- **C4-4.** The journal reference of Arbarello-Cornalba: J. Algebraic Geom. 5 (1996) 705-749, from the referee's memory, to be confirmed.

## 6. Item 5a, symplecticity of `M_0`: accepted

**Each step recomputed; all correct.**
- (1) `d/dv = (5 varsigma^4)^-1 d/d varsigma`, `zeta^(5i) = 1`, `nu_b = -b mod 5`.
- (2) `eta^-1 B` is symmetric iff `bm_(a-1) = bm_(3-a)`; all five `a` checked.
- (3)-(4) `dU_k/dv = -zeta^k varsigma^-1`, `dt/dv = -(4/5) t/v`; the ODE `(zeta^j - zeta^i) S + (4/5) t^2 S' = 0` and the lowest-coefficient argument.
- (5) `eps_a eps_(3-a) = 1` for every `a`, hence `Psi_0^T eta^-1 Psi_0 = 1`, `Psi_0 Psi_0^T = eta`.

**Evidence (VERIFIED, exact in `Q(zeta)`; `ref_g1_model.py`, `.log`).**
- The referee re-derived the T-system: `T_(nu_b - 1) F_b = bm'_b F_(b+1)`, `bm' = (1,-1,1,-1,-1)`.
- `F_0 = 1 - (3/20) t + (9/800) t^2 + (11571/80000) t^3 + ..`, and `F_b(0) = eps_b`.
- `S_ij(t) = delta_ij` for all `i, j` through `t^24`.
- Controls fail as required:
  - the diagonally rescaled, gauge-equivalent system `bm' = (2,-1,1,-1,-1/2)`, which is not `eta`-self-adjoint, with the same `eta`;
  - the pairing without `t -> -t`;
  - the pairing `a + b = 2 mod 5`.
- The first control shows that the self-adjointness hypothesis of step (2) carries the proof, not just the shape of the system.

**Notation.** `P_ij` (clashes with the gauge `P`), `B` (the "never" list of NOTATION table 4) and `S_ij` need local letters.

## 7. Item 5b, homogeneity of the model CohFT: accepted

**(a) Recomputed graph by graph.** Per ingredient:
- a vertex gives `varsigma^(2h_v-2+n_v+m_v)`;
- a leg `psi^k` gives `varsigma^(-4k)`;
- an edge monomial gives `varsigma^(-4(a+b+1))`;
- a kappa-leg gives `varsigma^(-1-4(m-1))`.

With `Sum_v (2h_v - 2 + n_v) = 2h - 2 + n`, every term of degree `D` carries `varsigma^(2h-2+n-4D)`.

**(b)** Correct: the exponent is `5e`.

**(c)** Correct. Under `varsigma -> zeta^-1 varsigma`:
- `Deltahat_i^(1/2) -> Deltahat_(i+1)^(1/2)`;
- `Rhat(zeta^-1 t) = Lam^-1 Rhat(t) Pi`, with `Lam = diag(zeta^(b+1))` and `Pi` the cyclic shift of colours;
- hence `M_0(zeta^-1 t) = Pi^-1 M_0(t) Pi` and `Psi^m_v -> Psi^m_v Pi`.

Add one sentence on naturality: `A_*(R.omega) = (A R A^-1).(A_* omega)` for an isometry `A`, which is immediate from PPZ15 Def. 2.2 and 2.13.

**(U)** The unit is correct: the lowered components of `Psi(Sum_i Deltahat_i^(-1/2) e~_i)` are `5 delta_(b3)`, so it is `phi_0`.

**Notation.** The `Lam` and `Pi` used here are the referee's local letters.

## 8. On GAPS_G1_5 s.6 (consequences for the design)

The referee agrees with items 1-8, subject to the following.
- **Item 1.** Pin Theorem 11.6 with the derivation of s.1.1, not "Cor. 11.7".
- **Item 3.**
  - Add the citations C3b-1 to C3b-3.
  - The single `\draftnote` should name the non-vanishing input `ch_(2k-1)(E) != 0` (s.4.5), not "injectivity read from 8.7 (ii)".
- **Item 6.** The removal of the five steps from the R7 list is agreed, with these marks:
  - `thm:set-qrr` (b) keeps one `\draftnote` until C3b-2 and C3b-4 are pinned;
  - `lem:set-dilaton` is conditional on `d(0) = 1` (C4-3);
  - "(QRR) unchanged" is a decision for the main session (s.4.6).
- **Item 7.** Add Mumford 1983 (DIGEST_D2_6 lists MR717614) once read; record the Givental misprint in `design/SOURCE_ERRATA.md`.
- **Item 8.** Extend the notation list with REFEREE s.2.3 and s.5, s.6 above.

## 9. Evidence on record (this folder; exact, finite, evidence only)

| file | what | result |
|---|---|---|
| `ref_g1_loc.py` / `.log` | degree-1 localisation: (M) at order `q` (generic weights and zeta), `<H^a,H^b>_(0,2,1)` at zeta | all PASS; 6 controls fail as required |
| `ref_g1_birkhoff.py` / `.log` | `prop:set-cyclicLR` with explicit `ell`, `q^30`, 15 periods; operator identity | all PASS; 4 controls fail as required |
| `ref_g1_model.py` / `.log` | `M_0(-t)^T M_0(t) = 1` through `t^24` in `Q(zeta)` | all PASS; 3 controls fail as required |
| `ref_g1_pdfglyph.py` / `.log` | glyphs of Givental Thm 11.6 (2) and Cor. 11.7 (1)-(2), p.40-41 | Cor. 11.7 prints `l_alpha`, Thm 11.6 prints `lambda` |

The scripts use Python 3.14 stdlib only and run as `py -I -B <script>` from the referee's scratch folder; the logs print relative names
only. The PDF reader takes the arXiv PDF as an argument and reads it as data. The prover's `g1_cyclic_check.py`, re-run from a scratch
copy, reproduces its log exactly.

## 10. Referee's reading record (full text, as PDF page images unless stated)

Read in full text:
- Givental alg-geom/9603021v2: p.23-27, p.35-42; font-level text of p.40-41.
- Teleman 0712.0160v3: p.11, p.39-41.
- PPZ15 1303.1043v3: p.3-4, 6-7, 9-10, 19-20, 27-30.
- Arbarello-Cornalba alg-geom/9406008v1: p.4-6.
- LP18 1803.01409v2: p.3-4, 23.
- GJR18 1812.11908v1: p.18-19, 26.
- Coates-Givental math/0110142v2: p.9-12, plus a text search of all pages for numbered statements.

Abstract pages only: alg-geom/9603021 (versions), math/9806133, math/9805097, alg-geom/9712011.

NOT read: Givental p.43-48 (the proof of the Main Lemma); Guo-Ross; Zagier-Zinger; Mumford 1983; Faber-Pandharipande; Madsen-Weiss;
Miller; Morita; Faber; Getzler-Pandharipande; Witten; GJR18 Lemma 5.2. Project sources consulted:
- DECISIONS; OUTLINE, OUTLINE_3, OUTLINE_REVIEW;
- NOTATION and NOTATION_2; LABELS and LABELS_2 (grep);
- DIGEST_D2_4, DIGEST_D2_6 (grep);
- PLAN_2026-11-03 "Key facts";
- S15_PROOF.md and theory/referee (grep for `d`).
