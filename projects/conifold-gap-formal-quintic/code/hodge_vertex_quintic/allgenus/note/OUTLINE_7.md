# Outline, part 7: s07 flat versus hatted; s08 the connection problem (the saddle series, Theorem P, the glue)

## s07_comparison — target 4 pp

Purpose: show, graph by graph and in every genus, that under (F_g) the outer term of the natural form equals the flat term of s05
(the comparison theorem), so that s05-s06 apply to the object of Theorem A.

Covers: DC2 s.1 D12-D15, s.2.4 (`lem:flh-*`, `thm:flh-main`, `rmk:flh-fail`), s.3 part 5, s.4 C9-C11, s.5 E12-E15, s.7 G4, W3, W4; DA A7,
C18; S15 s.0-s.5; FR s.0-s.5; OU2 s.1; ST2 s.2 step 6; ST2b K8, K9, O14-O20; GS s.2.3. Subsections `sec:flh-lemmas`, `sec:flh-thm`.

* `def:flh-outer` (defn; DC2 D12): a CohFT `Omega` on `(\mathbb V, eta, 1)` with dilaton-shifted ancestor potentials; `O(z) = Sum O_k z^k :
  \mathbb V -> \mathbb V'` symplectic, `O_0` invertible, `O_0^-1 1' = \varkappa^-1 1` with a scalar unit `\varkappa` (S15's "(U)" is printed
  "(unit)"; REV R8 (i): `\varkappa`, not `kappa`); edge `E_O`, leaf point `q_O = O^-1(z)(-z 1')`, `Cont_Gamma[E, q; Omega]` (including
  `|Aut Gamma|^-1`); TFT rescaling `(\varkappa.Omega)_(g,n) = \varkappa^(2g-2+n) Omega_(g,n)`.
* `lem:flh-frame` (lem; S15 Lemma 1): constant isometries move between the CohFT and the data. Proof [W] S15 s.2; FR s.1.
* `lem:flh-dil` (lem; S15 Lemma 2): `Cont_Gamma[E, q; \varkappa.Omega] = Cont_Gamma[E, q; Omega]`, `q` read as absolute (each potential at
  `q + z 1` of its own unit, FR s.1), from `lem:set-dilaton` with `\varkappa = d`; state `d(0) = 1` (`lem:gau-diag`; GS s.2.3; DEC Q5). Proof
  [W] S15 s.2; FR s.1.
* `lem:flh-hat` (lem; S15 Lemma 3 with the count of DA A7): for every leg-free genus-`g` graph, the hatted term (`def:set-outer`) equals
  `u^(2g-2)` times the outer term: `G_Gamma = s^(10(g-1)) T_Gamma`. Proof DC2 s.3 (i)-(iii) [W] S15 s.2 (edges unchanged since `nu_a + nu_b = 2`
  on eta-pairs; leaves get `s`; `z = s^4 zhat` gives `s^(4(3g-3))`; `xhat` gives `s^-(2g-2)`), cross-checked with OU2 s.1's count (DA A7
  steps 2-3 [W]). Correction DC2 G4: the sources prove it for "the description of the hatted code"; here it is a lemma about the
  definition `def:set-hat` and refers to no implementation. Correction FR C5/O18: `x = s varpi_0` must be defined.
* `lem:flh-comm` (lem; S15 Lemma 4): `R.(\varkappa.\TFT) = \varkappa.(R.\TFT)` for `R(0) = 1` (translation `\varkappa^-1 \ttr_R`, [PPZ15] Def. 2.5).
  Proof [W] S15 s.2; FR s.1.
* `lem:flh-inputs` (lem; S15 s.3 (I1)-(I4); FR s.1): (I1) `Omega^mod_(v(u)) = (Psi^m)_*(M.\TFT^m)`, `M = M_0(tau zhat)` (because `z/varsigma^4 =
  tau zhat` and `def:set-modcohft`); (I2) R0 in every colour: `thm:gau-R0` (iv) with `cor:gau-colours`; (I3) `prop:gau-symp`: `E_P` is a power
  series and `R' = P Psi^m` symplectic (under (I2) with `D = 1` in the window); (I4) `R'_0 = P_0 Psi^m`, `P_0^-1 1_Q = d^-1 1_m`, a TFT diagonal
  in an orthonormal basis is fixed by its unit, so `\TFT^Q = d.\TFT^m`, and `d(0) = 1` (`lem:gau-diag`; the formula `d = tau^(1/4) varpi_0` is
  not used). Proof [W] S15 s.3; (I1) is definitional here (OUTLINE_2 s.0.4 item 18).
* `thm:flh-main` (thm; "the comparison theorem"): let `g >= 2`, `Gamma` a leg-free stable graph of genus `g`, `F_Gamma` its term in the flat
  sum of `P_can.Omega^mod_(v(u))` (`def:red-flatsum`), `T_Gamma` its outer term (`def:set-outer`, `M = M_0(tau(u) zhat)` exact). Assume (F_g).
  Then `T_Gamma = F_Gamma`, hence `G_Gamma = u^(2g-2) F_Gamma`. Window: at genus `g` only `P` through `z^(3g-3)` enters (the one-vertex leaf
  `z^(3g-2)` needs `P^-1_(3g-3)`), so (F) is needed only as (F_g); genus 3 sees only the `z^5` coefficient (`[z^10] D_0 = ([z^5] D_0)^2/2` is
  forced by `prop:gau-symp`; REV R8 (b)). Proof DC2 s.3 steps 1-6 [W] (S15 s.4; FR s.1); step 5 [W-] (FR (I2); ST2 s.2); step 6 [W-]:
  `d(0) = 1` gives the u-adic convergence of the dilaton re-expansions (FR C5; `lem:tra-conv` (f)). Corrections FR C1-C6 (= ST2b O14-O20):
  leaf `T' = z y + q_(R')`; (F_g) for (F) (K9, W3); (I1) referenced; graph-local, no total potential, no TR (S15 s.4).
* `rmk:flh-fail` (rmk; S15 s.5, FR): where the comparison can fail ((1)-(6)); item (3): the composition of R-actions ([PPZ15] Prop. 2.4;
  for the unit-preserving action Def. 2.13 and Prop. 2.14) is needed only to identify the TOTAL with `F_g` (`cor:mai-gap`). Numbers only in
  appendix C.
* `rmk:flh-dilfactor` (rmk; S15 s.4; FR): `prod_v d^(2g_v - 2 + n_v) = d^(2g-2)`; only symplecticity of `P` beyond (F_g) is used.

Length 4 pp.

## s08_connection — target 10 pp

Purpose: define the full-Stirling saddle series `c(z; L)` and prove its rationality (DEC Q6), reduce `D_0 = 1` to one value at `q = 0`:
Theorem P identifies the canonical colour-0 WKB function with the full-Stirling one (an exact identity at the conifold end), and the glue
gives `D_0(z) = Phi^SP(z,1)/h(z)` using (QRR) once.

Covers: DD1 s.1 D1.1-D1.16, s.2.1-s.2.5, s.3, s.4, s.6 (addendum A-D), s.7; DD2 D6, GAM-1, GAM-2, s.3A (GAM-1), s.4; POLY s.1-s.6; PR;
ER s.2-s.4, s.7; EVID s.2.1-s.2.4, s.8; GAM s.0-s.3; GP s.4-s.5, s.8; GMR s.2; PB s.0-s.3; PS s.1-s.4; PS2 s.6-s.7; RS1 A4-A6; RS2 items
17-26; ST4 s.1; ST4b K1-K11, O1-O18; ST5b K3-K5, K8, K11-K13, O8-O14; FSL4 s.10; CP s.1.6; CP4 s.4; RA s.1.3; G5 s.1, G5_3 s.4 with G5R
s.1, G5R_2 s.4. Subsections `sec:con-problem`, `sec:con-form`, `sec:con-P`, `sec:con-glue`, `sec:con-finite`.

### sec:con-problem — the normalisation problem and the polynomial structure
* `rmk:con-problem` (rmk): the canonical gauge is normalised at the conifold point, the actual rows at `q = 0`; Property N(g) needs (F_g);
  `D` is a connection constant (a series built at `u = 0` compared at `u = 1`), not a statement at `q = 0` alone (PR C4; ST4b K5).
* `rmk:con-genus` (rmk; DD1 rmk:nrm-genus): `D_0 = 1 mod z^(N+1)` iff the coefficients through `z^N` vanish; with `log D_0 = Sum a_j z^(5j)` only
  odd `j` are independent (`D_0(z) D_0(-z) = 1`, from `prop:gau-symp` for the actual and the canonical gauge, PR C2), `a_j` first matters at
  genus `ceil((5j+3)/3)`. [W] EVID s.2.4; PR C2; GP s.8. Finite ranges only in appendix C.
* `lem:con-u0` (lem; POLY Lemma 1): for every `n >= 1` and every choice of the free coefficients, `ph_n(0) = Fs_(0,n) := [t^n] Fs_0`. Proof DD1
  s.3.1 steps 1-4 [W] (POLY s.2; GMR s.2.1); step 3 [W] (ER s.7 item 5: the `u^0` recursion is that of the model, `prop:set-modrows`).
  `rmk:con-model`: the conifold model does not fix the normalisation; the free coefficient `frak f_m` is the TOP coefficient.
* `prop:con-deg` (prop; degree bound, ER result B; not used by the proof of the main theorem): every `ph_n` (every colour, every choice of
  the free coefficients) has degree `<= floor(4n/5)` and every resonance is trivial; equivalently `A_n` lies in `L Q[L]` with degree `<= 4n +
  1`, and `frak f_m = [u^(4m)] ph_(5m)` is the value at `L = 0`. Proof G5_3 s.4.2-s.4.3 [G]: step 3 with the degree bound `deg RHS_n <= 4n + 5`
  of G5R_2 s.4.2 (F4.1: the weight argument `x : 1`, `z : -4`, `theta : 0`); step 5 as G5_3 s.4.3 (the left factor `1 - x^5`, the right side
  `(1 - x^5) rho_n` with `deg rho_n <= 4n`, resonance coefficient `[x^0] rho_n = [u^(4m)] r_(5m) = 0`; ER s.4's valuation lemma with the left
  factor of POLY s.1); step 4 no longer needs `lem:gau-poly`; it also gives existence for every choice of the free coefficients. F4.2: this
  completes ER s.4 (ER s.7 item 3 relabelled the bound). F4.3: write `L` (not `x`), local letters for G5's `M`, `P`, `Pi`, `E`, `rho_n`; the tags
  `(L0)`-`(L4)` are not printed. Equality is NOT claimed (computed only). Delete the resonance `\draftnote`.
* `rmk:con-points` (rmk): the three points `q = 0` (value `h`), `q = oo` (value `Ginf`, s09), `u = 0` (top coefficients) (GAM s.0; ER s.4).
  Dropped (D8): `rmk:con-ZZ`, `lem:con-p5k`.

### sec:con-form — the full-Stirling series, its rationality, the 1-form and `Phi^SP`
* `def:gam-c` (defn; DD2 D6, CP s.1.6, the definition of record; placed here by DEC Q6, label kept): `\Cstar := C \ {0, 1, zeta, .., zeta^4}`,
  `q(L) := (1 - L^-5)/5^5` (no fifth root); for `L` in `\Cstar`, a disc around `L` and branches of `log(5 sigma)`, `log(zeta^j - sigma)`:
  `lambda := -5 log(5L) + Sum_j log(zeta^j - L)` (`e^lambda = -q(L)`; local to s08-s09); `SS_0(sigma) := lambda sigma + 5 sigma log(5 sigma) -
  5 sigma + Sum_j [(zeta^j - sigma) log(zeta^j - sigma) - (zeta^j - sigma)]`, `SS_1 := (1/2) log(5 sigma) - (1/2) Sum_j log(zeta^j - sigma)`,
  `SS_2k := bk_k [(5 sigma)^(1-2k) + Sum_j (zeta^j - sigma)^(1-2k)]`, `A_St := exp(SS_1 + Sum_k z^(2k-1) SS_2k)`; `c(z; L)` := the normalised
  stationary-phase series (constant term 1) of `e^(SS_0/z) A_St d sigma` at `sigma = L`. It is the normalised saddle series of `Gamma(1 + 5s)
  prod_(j=0..4) Gamma(zeta^j/z - s) (-q)^s d sigma` (`sigma = zs`) with every Gamma factor replaced by `Gamma_St` (`def:fsp-stirling`): the
  `(log z)/z` coefficient is `-Sum_j zeta^j = 0`, the `log z` coefficient `2` (CP s.1.6; RF s.1). Taylor data with the local power sums
  `p_k(L) := Sum_j (zeta^j - L)^-k` (REV R8 (f)), `-1/SS_0''(L) = -L(1 - L^5)/5`. Colour `i`: the same at `sigma = zeta^i L`. Measure remark:
  `Gamma(1 + 5s)` carries `dy`, `Gamma(w)` carries `dx/x` (DD2 s.4; s.7 item 10).
* `prop:gam-rational` (prop; GAM-1; placed here by DEC Q6, label kept): `c(z; L)` does not depend on the disc or the branches; every `c_n`
  lies in `\AL = Q[L, L^-1, (L^5 - 1)^-1]`; `lambda` and the sign convention of `q` only fix the dictionary `L <-> q`. Proof DD2 s.3A steps 1-5
  [W] (CP s.1.6; CP4 s.4 (B1)-(B10); RA s.1.3), from `lem:fsp-univ`, `lem:fsp-const`.
* `lem:con-form` (lem; DD1 lem:thp-form (a)-(d)): (a) with Stirling for `Gamma(1 + 5s)` and `1/Gamma(1 + s - zeta^j/z)`, `sigma = zs`, `x = 1/sigma`,
  the exponent is `-WW(x)/z`, `WW = Xi/x + GG(x)`, `Xi := -log(1-u)`, `GG(x) = Sum_(m>=1) x^(5m-1)/(m(5m-1))`; `dWW/dx = -(Xi + log(1 - x^5))/x^2`;
  (b) the amplitude `frak a(x,z)`, `log frak a = -(1/2) log(1 - x^5) + Sum_k bk_k (xz)^(2k-1)[5^(1-2k) - 5 Sum_P binom(2k-2+5P, 5P) x^(5P)]`,
  independent of `u`, `frak a(0,z) = 1`; (c) hatted data `x = sX`, `z = s^4 zhat`: `WWhat(X; u) = Sum_(m>=1) (u^(m-1)/m)[1/X + X^(5m-1)/(5m-1)]`,
  critical point `X = 1` for every `u`, `WWhat''(1) = 5/(1-u)`, `WWhat(1) = (5/4)/tau`; (d) the series `c(z; L)` of `def:gam-c` (first Mellin-Barnes
  summand) and this form (the second, `q^s Gamma(1+5s)/prod_j Gamma(1 + s - zeta^j/z)`) give the same normalised series (DD2 GAM-2). Proof [W]
  POLY s.5; PB s.2.1; PR s.2; GMR s.2.2; (d) [W] Stage-4 K6, PB P9, CP s.1.7 (reflection `lem:fsp-refl` plus branch bookkeeping; RF-C5's direct
  route lacks the `(log z)/z` sentence for the second summand: use the reflection route). Correction ST4b O7/GMR C6: this phase is NOT UR
  s.8's; they share only the critical data.
* `def:con-PhiSP` (defn; DD1 D1.13): `Phi^SP(zhat, u)` := the normalised stationary-phase series (constant term 1) of `e^(-WWhat/zhat) frak ahat(X)
  dX` at `X = 1`, `frak ahat(X) := frak a(sX, s^4 zhat)`; `ph^SP_0 = 1`; the residue `VV := Res_(x=0) e^(-WW/z) frak a dx` in `Q[z^-1][[u]]`. The
  normalisation: division by the Gaussian factor times the amplitude value (DD1 s.4.5; GMR C2: the `z`-dependent factor is common to both
  sides, never "absorbed by `D(0) = 1`").

### sec:con-P — Theorem P
* `thm:con-Pa` (thm; Theorem P(a)): statement G5 s.1.2 (i)-(iv) [G]: with `\AL` (derivation `D`, `D L = (L^6 - L)/5`) and the embedding `\jmath :
  \AL -> Q((s))`, `L -> s^-1` (G5's `iota`, renamed): (i) `PF5^L(L c(z; L)) = 0` in `\AL[[z]]` (`PF5^L`: `zD` replaced by `zD + L`); (ii) `\jmath(c)(z) =
  Phi^SP(s^-4 z; s^5)` in `Q((s))[[z]]`, i.e. `\jmath(c_n) = s^(-4n) ph^SP_n`; (iii) `Phi^SP` solves HPF in `Q[[u]][[zhat]]`; (iv) every `ph^SP_n` is
  a polynomial in `u`, hence `c_n(L) = L^(4n) ph^SP_n(L^-5)` in `\AL` and `c_n` lies in `Q[L, L^-1]`. Proofs G5 s.1.3-s.1.5 [G] with G5R s.1:
  F1.1 (the identity (1.1) of G5 s.1.3 holds for amplitudes in `Q[sigma, z]`; Step 3 applies it only to `sigma^m` and `prod_k (5 sigma + kz)`),
  F1.3 (Step 4 over `\AL(zeta)`, then Galois invariance), F1.4 (local letters for G5's `c(eps)`, `Delta`, `g`, `psi_eps`, `e(eps)`); Steps 1-2 of
  G5 s.1.3 are `lem:fsp-generic` (appendix A). Inputs: `def:gam-c`, `prop:gam-rational`, `lem:con-form` (a)-(c), `lem:gau-HPF`, `lem:gau-poly`,
  `lem:fsp-univ`, `lem:fsp-ift`, `lem:fsp-dq` (iterated clause, needs `D_q frak a = 0`, RL-C2), `lem:fsp-transl` (b), `lem:fsp-stirfe`,
  `lem:fsp-cov` with `lem:fsp-homog`, `lem:fsp-generic`. Not used: Theorem U, (QRR), (CYC), `prop:con-deg`, `thm:gam-main`. Corrections PB P5
  (an infinite convergent sum of exact forms, `lem:fsp-ibp`), ST5b O12. Delete the domain `\draftnote`.
* `lem:con-cfu` (lem; the uniformising map, source "Lemma CFU"): unique `Y = X(1 + Sum_(k>=1) u^k Y_k(X^5))`, `Y_k` in `Q[T]`, `deg Y_k <= k`,
  and `V(u)` with `WWhat(X; u) = V/Y + Y^4/4 = W_V(Y)`; `V = tau^(-5/4)`, `Y_1 = (1 + T)/18`; `X = 1` goes to `Y_c = V^(1/5) = tau^(-1/4)`; the
  inverse is of the same form; commutes with `X -> zeta X`. Proof [W] PB s.2.4; POLY s.5; PS s.2.3. Name it "the uniformising map", with
  no attribution (the attribution question goes to `DRAFT_NOTES.md`; no `\draftnote` in the chain).
* `lem:con-red` (lem): with local letters for the sources' `beta(Y)` and the exact amplitude `Gamma` (REV R8 (e)): the 1-form `frak ahat(X(Y))
  X'(Y) dY`; `Y^n == (n-3) zhat Y^(n-4) + V Y^(n-5)` (`n >= 4`) modulo exact forms; reduction to `Sum_(m=-1..3) pi_m Y^m e^(-W_V/zhat) dY` plus an
  exact form, `pi_m` in `Q[zhat][[u]]`; Z/5 grading. Proof [W] PB s.2.5; POLY s.6 step 2; PR s.3; RL s.3.4.
* `lem:con-S2` (lem; source "(S2)"): `SP_(eta=1)[eta^(nu_b) e^(-(1/eta + eta^4/4)/t) d eta]`, normalised by `sqrt(2 pi t/5) e^(-5/(4t))`, equals
  `eps_b Fs_b(t)`. Proof [W] PS s.1 (i)-(vi) (PB s.2.6); inputs `lem:fsp-dq` with a moving critical point (PS-C1), `lem:fsp-ibp`, `lem:fsp-cov`
  and `lem:fsp-homog` (PS-C3), `prop:set-modrows` (uniqueness half).
* `lem:con-id` (lem; signs, weights, identification (a)-(d)): `SP_(Y_c)[Y^(nu_b) e^(-W_V/zhat) dY] = N_zhat eps_b tau^((1-nu_b)/4) Fs_b(tau zhat)`;
  `Phi^SP = Sum_b eps_b pi_(nu_b) tau^((1-nu_b)/4) Fs_b(tau zhat)`; for every regular gauge row 0 reads `Phi/varpi_0 = Sum_b tau^((1-nu_b)/4) G_b
  Fs_b(tau zhat)`, `G_b = s^(nu_b) P_0b`; triangular uniqueness; hence `G^SP_b = eps_b pi_(nu_b)/varpi_0`. Proof [W] PB s.2.6; PS s.1; PR s.3.
* `lem:con-relab` (lem; relabelling): `VV = s Res_X[e^(-WWhat/zhat) frak ahat dX]` coefficientwise, `[u^N zhat^j] \hat\cV = [u^(N + (1-4j)/5) z^j] VV`
  (REV R8 (e): `\hat\cV`, not `HH`); one summable family summed in two orders with finite fibres. Proof [W] PS s.2.2 (with PS-C4); PB s.2.7.
* `lem:con-key` (lem; key lemma): in `s Sum_m pi_m Res_Y[Y^m e^(-W_V/zhat) dY] = Sum_b pi_(nu_b) s^(-nu_b) I^C_(nu_b)(v(u), z)` the key `(u^E, z^beta)`
  receives `pi` only from `u`-orders `<= E + (4 beta - 1)/5` through finitely many terms; only `beta = 4 mod 5`. Proof [W] PB s.2.7; PS s.2.3.
  Correction PS-C5/RS1 A6: finite fibres give a formal sum; membership in `Q((z))[[u]]` follows from the coefficientwise identity with `VV` in
  `Q[z^-1][[u]]`, or from `lem:con-weight` (kept under D8 as one sentence with its one-line proof, PS s.2.3: the weight `5N - n - 4j >= 0`).
* `prop:con-res` (prop): `VV = varpi_0 (P^SP Y^C)_0` in `Q((z))[[u]]` in every conifold frame (REV R8 (e): `varpi_0`, not `C_0`; `P^SP := Q^SP
  (Q^m)^-1`), from `lem:con-relab`, `lem:con-key`, `lem:fsp-res` (ring PS-C2), the model residues `Res_(Y=0)[Y^m e^(-W_V/zhat) dY] = s^(-(m+1))
  I^C_m(v(u), z)`, `lem:gau-vanmod`, `lem:con-id`; signs cancel (`eps_b^2 = 1`), powers of `s` cancel. Proof [W] PB s.2.7; PR s.3; PS s.1-s.2.
* `thm:con-Pb` (thm; Theorem P(b)): `Phi^SP = Phi^can`. Hypotheses (complete list, DD1 thm:thp-b): `thm:gau-U` for `(A, A~)` with the hypotheses
  of `prop:gau-quintic`, `prop:gau-equalprod`; `thm:gau-G` with `lem:gau-uniq`, `lem:gau-R1`, `lem:gau-vanmod`, `lem:gau-vanq` (a); `cor:gau-cycsys`;
  `prop:set-modrows` (uniqueness half); `lem:gau-poly`, `thm:con-Pa`, `lem:con-form`..`prop:con-res`; appendix A (`lem:fsp-ibp`, `lem:fsp-cov`,
  `lem:fsp-cfu`, `lem:fsp-dq`, `lem:fsp-res`, `lem:fsp-stirfe`, `lem:fsp-transl` (b), `lem:fsp-homog`); the conifold frames (`varpi_0(0) = 1`,
  `t_c'(0) = 1`); `ph^can_0 = 1` (`lem:gau-diag`). Not used: (QRR), the finite route, the degree bound, Theorem Gamma. Proof DD1 s.3.3 steps 1-7
  [W] (PB s.2; PS s.1-s.3): step 5 both solve HPF, `Phi^SP = Phi^can \mathsf s(z)` (`\mathsf s` a local letter, in `1 + z^5 Q[[z^5]]`), so
  `P^SP = \mathsf s(z) P_can` and `(P^SP Y^C)_0 = \mathsf s(z)(Y^C_Q)_0` by `thm:gau-G` for `P_can` (PS-C6); step 6 `prop:con-res`; step 7 at `u^1`:
  `VV_1 = -1/z`, `[u^1](Y^C_Q)_0 = -1/z`, `[u^0](Y^C_Q)_0 = 0` (PS-C7), `varpi_0(0) = 1`: `\mathsf s = 1`. Corrections ST5b O10, PR C6, PB P1, PS-C1..C7.
* `cor:con-frame` (cor): `Phi^can` does not depend on the frame (PB s.0; `Phi^SP` is frame-free). `rmk:con-710` (rmk): at `u^2` the residue
  produces `-(1/2 + 1/5)/z = -(7/10)/z`, the vanishing period's `7/10` (POLY s.5; GP s.5). Omitted: GAM's "Theorem A(i)" (DD1 rmk:thp-Ai).

### sec:con-glue — the glue
* `thm:con-glue` (thm): `h(z) D_0(z) = Phi^can(z, 1)` in `Q[[z]]`; with `thm:con-Pb`, `D_0(z) = Phi^SP(z, 1)/h(z)` exactly. Proof: `thm:gau-R0` (iv) at
  `u = 1` (`s = 1`, `zhat = z`, `L = 1`), `Phi^act(z, 1) = Sum_n A^act_n(q=0) z^n = h(z)` by (QRR) in the form of `thm:set-qrr` (G1_4 s.3.1); [W] PS
  s.4.2 (e); PB s.3. Not used: the finite route, the degree bound. No `\draftnote`.
* `cor:con-A2` (cor): `D = 1` iff `Phi^SP(z, 1) = h(z)` (POLY s.6 Corollary; PR s.2; with `cor:gau-colours`); the form `c(z; 1) = h(z)` is
  `cor:gam-equiv` (s09).
* `rmk:con-why` (rmk; DD1 rmk:glu-why): every property used so far (PF5, symplecticity, Z/5, polynomiality) is invariant under `D -> D exp(odd
  series in z^5)`; an exact identity is needed. [W] GP s.8; GMR s.2.8; ST4b O10. That only the symplecticity of `G` detects `D = 1` among the
  outputs of the natural form is computed (G5R_2 s.3.3; appendix C).

### sec:con-finite — the finite route (kept short under D8; not used by the main proof)
* `lem:con-N` (Lemma N, with (CYC) among its inputs, PS-C9) and `thm:con-conn` (connection form (i)-(iii)): statements and proofs DD1 s.2.3,
  s.3.2 [W]; `lem:con-P` only if the proof of `thm:con-conn` needs it. They explain the computations of appendix C; further conifold frames give
  no independent evidence (ER s.3; ST4b K7). Dropped (D8): `prop:con-pfree`, `prop:con-row`.

Length 10 pp.
