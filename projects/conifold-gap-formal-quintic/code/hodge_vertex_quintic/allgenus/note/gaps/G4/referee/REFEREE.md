# Referee report on gap cluster G4 (s06, the model), part 1: verdict, item 1, item 2 (key gaps/G4/referee, 2026-11-15)

Author-facing working note; not part of the note; nothing here enters a `.tex` file as is. AI referee, AI-checked only; no human has
read it. Labels as in the project: PROVED, VERIFIED (exact computation; range and field stated), CONJECTURAL, SKETCH, OPEN.
Object: `NOTE/gaps/G4/GAPS_G4.md`, `_2` .. `_5`, `g4_field.py`, `g4_rtr.py` (+ `g4_rtr_10.log`), `g4_w11.py` (+ log), `g4_period.py`
(+ `g4_period_5.log`), `pdftext.py`; all read in full. Part 2 (`REFEREE_2.md`): items 3 and 4, the findings F1-F10, the referee's
own exact checks, the list of corrections, files. No file outside this folder was created, modified, moved or deleted.

**Read for this report.** Design: `OUTLINE_6.md` (sec:mod-tr .. sec:mod-V), `OUTLINE_3.md` item 10 (`def:set-model`,
`prop:set-modrows`, `def:set-M0`, `def:set-modcohft`), `design/OUTLINE_REVIEW.md` (R3, R7), `design/DECISIONS.md` (D5, D6),
`design/DIGEST_C2_7.md` s.7 (G1-G3, W5), the NOTATION/LABELS rows used. Project sources: TRM (`allgenus/tr/TR_MODEL.md`, all),
TRRF (`allgenus/tr/referee/REFEREE.md`, all), XY (`allgenus/xyswap/XYSWAP_MODEL.md`, all), XYR (kernel paragraph, its l. 37-38), LL2
(`stage2/logleg/LOGLEG_PROOF_2.md` s.3), LG (`natural_form/lg/LG_MODEL.md`, symplecticity lines), UR l. 375. Literature, read in
the TeX source of the arXiv version named (downloaded to a scratch folder; equation and statement numbers recounted from the
source with the reading aids `eqcount.py`, `eo_assemble.py`, `eo_appA_count.py` of this folder, whose counts were corrected by hand
where a `\substack{..\\..}` produced spurious line breaks):
* DNOPS = arXiv:1509.06954v2 (25 Dec 2016; J. Inst. Math. Jussieu 18 (2019)): s.2.1 Remark 2.1; s.2.3; s.2.5 (2-26)-(2-29);
  s.2.6 (2-30)-(2-35); s.3 Thm 3.1; s.4 Thm 4.1, Cor. 4.2, Prop. 4.3; s.5 Thm 5.1 and its proof, Prop. 5.2; s.8 (8-1), Thm 8.1.
* DOSS = arXiv:1211.4021: s.2.3 (leaf, dilaton-leaf, edge, vertex weights; the graph formula), s.4 Theorem 4.1.
* PPZ = arXiv:1303.1043v3 (Introduction = section 0): Def. 2.2, Def. 2.5, Def. 2.13, Prop. 2.14.
* EO = arXiv:math-ph/0702045v4: s.4.1 (4-1)-(4-9); s.4.4 Thms 4.1-4.7, Cor. 4.1 (4-26)-(4-29); App. A, (A-1), (A-5), (A-25).
* ABDKS = arXiv:2212.00320v3 (30 Jul 2024, the current version; JEMS 2025): s.1.1 (1.1); s.7, Thm 7.1, (7.5)-(7.12), Rems 7.2-7.4.
* LPSZ = arXiv:1504.07439: its Theorem `thm:CohFT` and `alpha^Top` (only for the TFT normalisation; context).

## 0. Overall verdict: **accepted-with-corrections**

All four resolutions are correct mathematics: I checked every step by hand and re-ran the decisive identities with independent code
along different computational routes (part 2, s.3). One correction is **required**: G4 copies a sign misprint of DNOPS s.2.3 into
the statement of the dilaton leaves (item 1 (b)) and calls the discrepancy with the outline a "sign convention"; it is not, and the
outline's sign is the right one (finding R-1, decided by an exact computation at `(1,1)`). The other corrections are wording and
pinning. After R-1 the four `\draftnote`s D6/G2, G3, W5 and D5/G1 of OUTLINE_6 can be deleted. Nothing in G4 needs the label OPEN.

| item | G4's resolution | verdict | in one line |
|---|---|---|---|
| 1 `prop:mod-Rtr` (D6, G2; R3) | PROVED mod DNOPS/DOSS | **accepted-with-corrections** | (a), (2-32), `lem:set-M0sympl` and the transfer are correct line by line; the dilaton leaf in (b) has the wrong sign (R-1, required); the degree sentence of the transfer is wrong as worded (R-2); the operative quoted statement should be pinned as DOSS Thm 4.1 via DNOPS s.2.6 (R-3) |
| 2 `lem:mod-xfun` (G3, M6) | PROVED (self-contained; EO pin) | **accepted** | Lemma T and the per-branch-point deduction are complete; the EO pins (Thm 4.4 = (4-23), (4-24); Cor. 4.1 = (4-26)-(4-29); (A-25)) are right; F2 is a genuine precision |
| 3 `thm:mod-Lreg`, `lem:mod-abdks`, `prop:mod-period`, `lem:mod-vanper` (W5) | PROVED mod ABDKS Thm 7.1, EO Thm 4.2 | **accepted-with-corrections** (minor) | every step checked; Step 4's finite identity replaces the Borel argument correctly; the edge weight (7.7) is mis-transcribed (R-4, harmless) |
| 4 `thm:mod-V` (ii), base `h = 1` (D5, G1) | PROVED (finite derivation by hand) | **accepted** | the local formula and the numbers re-derived by hand; one implicit step (zero residues of `omega_(1,1)`) to be said (R-6) |

## 1. Item 1: `prop:mod-Rtr` (GAPS_G4_2) -- accepted-with-corrections

### 1.1 The quoted statements, checked in the source
* DNOPS (2-31) reads `R^-1(zeta^-1)^j_i = -(2 pi zeta)^(-1/2) Int B(w_i, w_j)/dw_i|_(w_i=0) e^((x(w_j) - x(c_j)) zeta)` with the
  coordinates `x = x(c_i) - w_i^2/2` of s.2.6; (2-32) reads `Sum_k (R^-1(zeta^-1))^i_k Delta_k^(-1/2) = sqrt(zeta/2pi) Int dy(w_i)
  e^((x(w_i) - x(c_i)) zeta)`; (2-30) `Delta_i^(-1/2) = dy/dw_i(0)`; (2-34) `xi^i = Int^z B(w_i, .)/dw_i|_0`; (2-35) as G4 states.
  G4's transcription (with `zeta -> lam`) is exact, and so is its index reading (upper index = row = output; leaves carry
  `R^-1(psi) e_alpha = Sum_beta (R^-1)^beta_alpha e_beta`). The kernel of (2-28) is `Int_z^(sigma z) B / (2 (omega_(0,1)(sigma z) -
  omega_(0,1)(z)))`, equal to the engine kernel (both orientation and denominator flip), as `conv:mod-kernel` says.
* The sentence of s.8 before Thm 8.1 is quoted correctly. DNOPS justify "reconstruct a CohFT" in the proof of Thm 5.1 by the
  compatibility of Thm 4.1, which for a rational curve with `x(oo) = oo` is Prop. 4.3 (1) (`eta = 0`); for `f = 4 W_v` this holds
  (`dy/dx = p^2/(4(p^5 - v))` has poles only at the branch points and no residue at `oo`). Thm 8.1 itself (the Frobenius-manifold
  identification) is not used by G4's route, as F7 says. Remark 2.1 (Givental's convention is the transpose) is quoted correctly.
* R3: at `a = (0,0,0,0,4v)`, `f = 4 W_v`, `f' = 4(p^5 - v)/p^2`, `f''(p_c) = 20 p_c^2 != 0`, critical values `5 p_c^4` pairwise
  distinct. Correct.

### 1.2 Line-by-line check
* `lem:mod-laplace`. (L1) re-checked on `w^m`, `m != -2` (odd `m >= 1`, odd `m <= -1` including `m = -1`, even `m`); (L3) re-derived
  (`w_j(v)(k p') = k^2 w_j(v')(p')`, the sign fixed by the `v`-independent root `sqrt(Lambda)`), giving `L_v[p^m dp](lam) = k^(m-1)
  L_(v')[p'^m dp'](k^4 lam)`; (L4) re-checked on `w^(2n) dw`. The finite-part value `G_(-2) = -lam` agrees with the contour
  integral over `R + i eps` (one integration by parts), and it is the reading under which `R^-1 = I + O(lam^-1)` holds; it agrees
  with DOSS's prescription (`delta_ij` from `dw dw'/(w - w')^2` plus the transform of the regular part). The "Meaning" sentence is
  heuristic only (`w` is a local coordinate, the real line leaves the disc); the definition of record is the formal one (R-7).
* Step 1: correct; the additive constant in `xi^i` is irrelevant (`L_j[C dx_kappa] = C G_1 = 0`). Step 2: correct
  (`-kappa/Lambda = 1/5`).
* Step 3, re-derived: (L1) with `F = p^m` gives (Mod2) `m z J_(m-1) - J_(m+3) + v J_(m-2) = 0`, `z = -1/(kappa lam)`; (L3) at `k = 1`
  gives `5 v d_v J_m = (m-1) J_m + 4 lam d_lam J_m`; (L4) with `x_kappa - x_kappa(varsigma_j) = kappa(p^4/4 + v/p) - (5/4) kappa
  varsigma_j^4` and (Mod2) at `m+1` give `z d_v J_m = -J_(m-1) + J_m/varsigma_j`; with `d_v[W_v(varsigma_j(v))] = 1/varsigma_j` this is
  `z d_v Yt_m = -Yt_(m-1)`. The five relations `z d_v Y_b = bmod_b Y_(b+1)` follow (checked one by one, incl. `Y_2 = -v Yt_(-2)`).
  With `d_v U_j = -zeta^j/varsigma`, `d_v t = -(4/5) t/v`: `z d_v[e^(U_j/z) varsigma^(nu-1) F(t)] = e^(U_j/z) varsigma^(nu-2) T^(j)_(nu-1) F`,
  `T^(j)_k = -zeta^j + (t/5)(k - 4t d/dt)`, so the relations are `F_1 = T_(-1) F_0`, `F_2 = -T_(-2) F_1`, `F_3 = T_2 F_2`, `F_4 = -T_1
  F_3`, `F_0 = -T_0 F_4` (the powers of `varsigma` cancel, `bmod_1 = -1/v` supplying `varsigma^-5`). Initial values and the colour
  rotation: correct. Uniqueness: on `t^n`, `(1 + T_0 T_1 T_2 T_(-2) T_(-1)) t^n = -4 n c^4 t^(n+1) + O(t^(n+2))`, `c = zeta^j`, so the
  closing relation fixes `f_n` from `f_(n-1) .. f_(n-4)`; existence and uniqueness given `F_0(0)` follow from this alone.
* Step 4 (`m_b = -b` mod 5, `varsigma_i^(1 - m_b) varsigma^(m_b - 1) = zeta^(-i(b+1))`) and Step 5: correct. Result `X(lam) =
  M_0(t)^T`, `t = -1/(kappa varsigma^4 lam)`, hence `R_TR(w) = M_0(w/(kappa varsigma^4))`: **PROVED** (mod the quoted DNOPS/DOSS data).
* `lem:set-M0sympl`: correct. `(eta^-1 A)_(ac) = (1/5) bmod_(3-a) delta_(c, 4-a)` is symmetric iff `bmod_(3-a) = bmod_(a-1)`, true
  for `a = 0..4`; `nu_a + nu_b = 2` on the `eta`-pairs `(0,3), (1,2), (4,4)` (indices mod 5); `Delta' = -(zeta^j - zeta^i)/varsigma`;
  the induction and `Psi_0^T eta^-1 Psi_0 = 1` (`eps_b eps_(3-b) = 1` for all `b` mod 5) are right. The formal manipulations take
  place in `e^(Delta/z) C(v)((z))` with `d_v` acting by Leibniz: legitimate. It proves what the sources only VERIFIED (LG s.3,
  through `zhat^40`). GAPS_G4.md s.1 should write `eta(H^a, H^b) = 5 delta_(a + b = 3 mod 5)` (R-5).
* (b): (2-32) holds identically (both sides `Lambda^(-1/2) varsigma^-1 Rhat^mod_0^(i)(t)`): correct. **The dilaton leaf is not**: R-1.
* s.6, transfer `kappa = 4 -> 1`. `lem:mod-scale` is correct (same branch points and involutions, kernel `/ (kappa mu)`, induction on
  `2h - 2 + n`; `F_h` through `Phi_ydx -> kappa mu Phi_ydx`). (i): the TFT normalisation `omega(e_j^(x n)) = Delta_j^((2h-2+n)/2)` is
  the right one (it is DOSS's vertex weight `Delta^((2g-2+val)/2)` in their graph formula, and LPSZ's `alpha^Top` in their
  coordinates; it is also forced by `y -> a y` on the Airy curve), but G4 gives no reference for it (R-3). (ii), (iv): correct.
  (iii): the final formula `alpha^4 = 2^(2h-2+n) Sum_D 4^-D [alpha^1]_D` is right, the explanation is not (R-2).

### 1.3 Findings for item 1
* **R-1 (required). The dilaton leaf has the sign of PPZ, not the printed sign of DNOPS s.2.3.** DNOPS s.2.3 decorates each dilaton
  leaf by `-psi (Id - R^-1(psi)) e_1`. PPZ's unit-preserving action (Def. 2.13: `R.Omega = R T Omega`, `T(z) = z 1 - z R^-1(z) 1`;
  the translation legs of Def. 2.5 carry `T(psi)` itself, not `R^-1(psi) T(psi)`) puts `+psi (1 - R^-1(psi)) 1` there, and so does
  DOSS: their dilaton-leaf weight `[z^(k-1)](-R(-z)^i_1)`, `k >= 2`, sums to `psi (1 - R(-psi)) 1`, and DOSS's `R(-z)` is DNOPS's
  `R^-1(z)` (DNOPS Remark 2.1; DOSS's ordinary leaves carry `R(-z)`). TR decides between the two: `rf4_check_b.py` computes
  `omega_(1,1)` on `x = kappa W_1` by the recursion itself (residues at all five branch points, local involution by a fixed point,
  no symmetry) and the graph sum (2-35) at `(1,1)` with DNOPS's own `R^-1_1` from (2-31) (not from `M_0`). The `(1,1)` graph sum
  (derivation in part 2, s.3.2) matches TR at all four test points for `kappa = 1` and `kappa = 4` with `+T`, and fails at all four
  with DNOPS's printed `-T`; the control `R^-1 -> R` fails too. (TRM s.2.3 has the right sign: "omega_(0,1) = y dx produces exactly
  the frozen dilaton leaf T(z) = z (1^ - M_0^-1(z) 1^)"; the TRRF and MVertex6 comparisons used it.) Hence:
  - GAPS_G4_2 s.1 (b) and s.4: replace "`-psi (I - M_0^-1(psi/(kappa varsigma^4))) 1hat`" by
    "`psi (1hat - M_0^-1(psi/(kappa varsigma^4)) 1hat)`, the translation `T` of [PPZ15] Def. 2.13";
  - GAPS_G4_2 s.7, third bullet: delete "the same translation with the opposite sign convention for the leaf". With legs
    `R^-1(psi)` and edges as in PPZ Def. 2.2 the sign of `T` is not a convention; the outline's `z(1hat - M_0^-1(z) 1hat)` is right;
  - GAPS_G4.md s.1 and GAPS_G4_2 s.7, second bullet: "the R-action of DNOPS s.2.3 = [PPZ15] Def. 2.2" -> "the leaves and edges of
    DNOPS s.2.3 are those of [PPZ15] Def. 2.2; with dilaton leaves the action is [PPZ15] Def. 2.13 (with Def. 2.5)";
  - record in `design/SOURCE_ERRATA.md`: external misprint, DNOPS arXiv v2 s.2.3, dilaton-leaf sign (journal version not checked).
  Nothing else changes: (a), (2-32), the transfer and the identification `alpha = Omega^mod_v` hold verbatim with PPZ's `T`, which
  is the form the note uses everywhere else (`def:set-graphsum`, NOTATION_2 s.9 item 1).
* **R-2 (recommended). Transfer, s.6 (iii).** "(push-forwards preserve degree)" is false as written: forgetting a dilaton point lowers
  the degree by one, gluing an edge raises it by one. The correct bookkeeping: under `R -> R(./4)` a leg term `psi^k` gets `4^-k`, a
  dilaton term `psi . r_k psi^k` (degree `k + 1` upstairs, `k` on `Mbar_(h,n)`) gets `4^-k`, an edge term `psi'^a psi''^b` (degree
  `a + b` upstairs, `a + b + 1` downstairs) gets `4^-(a+b+1)`; so every contribution of degree `D` on `Mbar_(h,n)` gets `4^-D`. State
  it once as `[R(./c).omega]_(h,n) = Sum_D c^-D [R.omega]_(h,n)^(D)`.
* **R-3 (recommended). Pin the operative quoted statement.** What item 1 uses is: for this compact curve with simple branch points,
  `R` from (2-31) in Givental's group and (2-32) holding, TR equals (2-35) with `alpha = R.omega_TFT` (unit-preserving action, TFT
  `Delta^((2h-2+n)/2)`). In DNOPS this is not a numbered statement but the s.2.6 recollection of DOSS plus the sentence before
  Thm 8.1 (justified via Thm 4.1 / Prop. 4.3 (1) in the proof of Thm 5.1). Cite it as "[DOSS] Thm 4.1, as recalled in
  [Dunin-Barkowski:2015caa] s.2.6, with Prop. 4.3 (1) for the compatibility", plus DOSS's vertex weight for the TFT. DOSS is in the
  root file `refGW/refs.bib` as `DuninBarkowski:2012bw` (CMP 328 (2014)); BIBPLAN does not list it yet. The substance of (Q) is unchanged (it is the
  theorem DNOPS rest on), but its wording is not: the (Q) entry "[Dunin-Barkowski:2015caa] Thm 8.1" names a statement that the
  chain does not use (part 2, s.5 item 3).

### 1.4 Evidence (VERIFIED; independent of G4's code)
`rf4_check_a.py 9` (`rf4_check_a_9.log`), exact over `Q(zeta5)`: (2-31) computed by residue extraction `[w^m] f = Res f w^(-m-1)`
in the coordinate `t = p - varsigma_j` (no series reversion; G4 reverts), `kappa = 1, 4, -1, 1/2`, through `lam^-9`: `X = M_0(t)^T`
(0 mismatches of 250 each), (2-32) with its right side computed independently (0), the IBP form (1) (0), `M_0(-t)^T M_0(t) = 1` and
`Rhat(-t)^T eta^-1 Rhat(t) = 1` (0); controls fail (no transpose 160, `t -> -t` 125, `eta^-1 -> 1/5` fails). `rf4_check_b.py`: R-1.
These are evidence; the proof is G4_2 s.3-s.6 as checked above.

## 2. Item 2: `lem:mod-xfun` (GAPS_G4_3) -- accepted

* Lemma T: `prod_i (xi - p^i(X)) = xi^5 - 4 X xi + 4 v = 4 xi (x(xi) - X)` (the fibre equation of `x = W_v`), so `Sum_i B(p^i(X), xi)
  = dX dx(xi)/(X - x(xi))^2` (this is EO (A-1)); its integral from `sigma_b q` to `q` is `dX[1/(X - x(q)) - 1/(X - x(sigma_b q))] = 0`
  (this is EO (A-5)). The contour representation of the recursion is valid for `p_1` outside `U` provided the discs `U_b` are
  chosen disjoint from `J` (implicit in "p_2 .. p_n off the branch points"; say it, R-8); the integrand has in `U_b` only the pole
  at `varsigma_b`; the finite sum over sheets commutes with the contour integrals; `Tr omega` is meromorphic on `P^1` and vanishes
  on an open set. Correct, and it covers every `(h, n)` with `n >= 1`, `2h - 2 + n > 0` (also `(0,3)`).
* Per-branch-point deduction: the fibre over `X_b` is `varsigma_b` (double) and three unramified points (critical values pairwise
  distinct; `p = 0, oo` map to `oo`); `omega(., J)` is holomorphic there (EO Thm 4.2, or (ii) with shrinking discs, slot 1, then
  symmetry EO Thm 4.6); the trace formula for residues gives `Res_(varsigma_b) f(x) omega = 0`. Correct and complete.
* Consequences: iterated residues at branch points commute (no diagonal poles); rays and circles can be exchanged with the slot-1
  residue by Fubini (compact contours; `omega` is holomorphic at `0`, `O(dp/p^2)` at `oo` in each slot). Correct.
* The EO pins, checked in the TeX source of v4: Thm 4.2 (poles only at branch points, `(k,g) != (1,0)`), Thm 4.4 = (4-23), (4-24),
  Thm 4.5 = (4-25), Thm 4.6 (symmetry), Cor. 4.1 = (4-26)-(4-29), (4-7) `W_1^(0) = 0`; App. A: (A-1), (A-5) as above, (A-25) is
  "for any rational `R(x)` with no pole at `x(a_i)`: `Res_(a_i) R(x(p)) W = (1/2) Res_(a_i) R(x(p))(W(p) + W(pbar)) = 0` due to theorem
  4.4". G4's remark that EO tacitly need the other points over `x(a_i)` to be regular is right. The kernel translation (EO's
  `dE_q(p) = (1/2) Int_q^qbar B`, `omega(q) = (y(q) - y(qbar)) dx(q)`: minus the engine; `W^EO_(g,n) = (-1)^n omega_(g,n)`) re-derived.
* F2 is right and matters for the wording only: `Res_(p=0) x^5 omega_(1,1) = -15/4` (by hand: `x^5 = p^-5 (1 + p^5/4)^5`,
  `omega_(1,1) = -(15/4) p^4 (1 + O(p^5)) dp`); `rf4_check_b` B3: the 35 residues `Res_(varsigma_b) x^k omega_(1,1)` vanish, 40 of 45
  control residues of `p^k omega_(1,1)` do not (a consistency check, as G4 says; LR F4 applies).
* Remark (1) of G4_3 s.2 (coinciding critical values) and remark (2) (the [DW] route) are correct and not needed.
