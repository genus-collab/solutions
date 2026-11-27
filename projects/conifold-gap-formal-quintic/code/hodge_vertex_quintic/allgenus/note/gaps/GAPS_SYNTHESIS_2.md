# Gap phase synthesis (key gaps/critic), part 2 of 3: conditionality, what is open, changes to the OUTLINE files (s01-s03)

Continuation of `GAPS_SYNTHESIS.md` (verdict s.0, status table s.1, adjudications s.2). Same conventions.

## 3. Exact conditionality of the main theorem after this phase

Theorem A (`thm:mai-main` = `thm:int-main`) is PROVED modulo exactly the following inputs; every other step of the chain C1-C13 (with the new
row C1') is written out in the note's sources or in the gap notes and refereed. Numbers are those of the arXiv versions named.

**3.1 (Q) quoted theorems** (statements read in the full text in this phase or by the designer; TR only on the model curve):
* [Eynard:2007kz] math-ph/0702045v4: Def. 4.2 (4-9) with the kernel (4-3) (`W^EO_(h,n) = (-1)^n omega_(h,n)`); Thm 4.2 (poles only at the
  branch points); Thm 4.6 (symmetry); Thm 5.1 (5-14)-(5-15) with s.5.3 (5-27)-(5-30) and Thm 4.7 (4-30), as pinned by the designer (not
  touched in this phase). Thm 4.4 (4-23)-(4-24), Cor. 4.1 (4-26)-(4-29) and (A-25) become optional pins: `lem:mod-xfun` is proved without them.
* [DuninBarkowski:2012bw] (DOSS, arXiv:1211.4021) Thm 4.1, with its graph weights (s.2.3), as recalled in [Dunin-Barkowski:2015caa]
  (1509.06954v2) s.2.6 (2-30)-(2-35), with the action of its s.2.3 (dilaton leaves with the sign of [PPZ15] Def. 2.13; s.2.6), the
  compatibility of its Prop. 4.3 (1) and the sentence before its Thm 8.1; at `f = 4 W_v`, i.e. (8-1) with `n = 4`, `a = (0,0,0,0,4v)`
  (critical points `p^5 = v`, critical values `5 p_c^4`, pairwise distinct). This replaces "DNOPS Thm 8.1" (G4R R-3; R3).
* [ABDKS25] (2212.00320v3) Thm 7.1: (7.6), (7.7) (with its factor `1/hbar^2`), (7.9), (7.10); the kernel (1.1) (EO's printed kernel).
* [PPZ15] (1303.1043v3): Def. 0.3 (iii); Rem. 2.1; Def. 2.2 with the convention sentence after it (p.20); Prop. 2.4; Def. 2.5; Prop. 2.9;
  Prop. 2.12; Def. 2.13; Prop. 2.14.
* [Hock25] (2502.18115v2) Cor. 1.2, (1.3): only while route 1 is the proof of record of the Gap Theorem (D4). Route 2 is now complete (G4
  item 3: ABDKS (7.10) and the exact identity `S_h = -4^h B_2h/(2h)`), so Hock can leave the conditionality (s.7 Q4).

**3.2 (QRR)**, unchanged as a quoted input: `R(q = 0) = diag(h_(w_0), .., h_(w_4))` in Givental's convention (the columns `R e~_i e^(u_i/z)`
solve the QDE; legs `R^-1(psi)`), now used for **every colour**: row 0 of colour `i` at `q = 0` (G1 s.3.1) enters `def:set-rows`,
`lem:gau-colrows` Step 3 (s03) and `thm:con-glue` (s08). Sources as in the design: [GJR18] (17), [GJR17] (19)-(20), [Giv01] s.2.2 Thm 1
with the Remark, [MR2276766] Thm 1, Cor. 1. Its derivation from the CohFT (G1 s.3.2) is written and refereed, with two pins missing (O1).

**3.3 (M)**, replacing (S1): [MR1408320] (alg-geom/9603021v2) Theorem 11.6, Calabi-Yau case `n = 4`, `r = 1`, `l_1 = 5`, `lambda' = 0`,
rewritten for `s_i` (s.2.1), with the set-up of s.9 (p.23-27) and Prop. 11.4; specialised to `lambda_i = zeta^i`: `J(T(q) H, z) = I(q, z)/C_0(q)`,
`T = I_1/C_0`. Hypotheses: genus-0 convexity of `O(5)`, polynomiality in `lambda` (G1 s.1.4). `\draftnote{IMRN numbering not compared}`.

**3.4 (G0)**, new name for what `prop:set-cyclicLR` uses besides (M): the genus-0 twisted classes form a CohFT with unit for `eta` ([GJR18]
s.4.1; [PPZ15] Def. 0.3), polynomial and `S_5`-symmetric in `lambda` (proper moduli spaces; `E_(0,n,d)` a vector bundle of rank `5d + 1`),
with the string, divisor and splitting axioms, the genus-0 TRR, and the quantum differential equation ([GJR18] s.5.2 eq. (14)) (G1 s.2.1).

**3.5 (T)**, replacing (S2): [MR2917177] Theorem (2) in paragraph (1.15) with [GJR18] s.4.1, Thm 4.1: the CohFT of the formal quintic is the
graph sum of `def:set-graphsum` with the R-matrix of `thm:set-qrr`. In Theorem A it identifies that R-matrix (C1); in Corollary B it identifies
the total with `F_g`.

**3.6 (IFT)**, replacing (S3): the holomorphic implicit function theorem in LMI step 2, unless the note uses `lem:fsp-ift` (formal, proved
in appendix A), which suffices (LA s.2.1). Recommended: use the formal lemma and drop (IFT) (s.7 Q2).

**3.7 Standard results cited** (read in the full text): Arbarello-Cornalba alg-geom/9406008v1, formula (1.7) and `kappa_0 = 2g - 2 + n` (p.5),
for `lem:set-dilaton` (i). Definitional: [LPfq] s.0.2 (5) (the formal quintic), (10) (fibre weight `5 lambda_i`).

**3.8 Not used by the chain** (may be cited elsewhere): `Int_(Mbar_(1,1)) psi_1 = 1/24` (only the genus-1 log term of `lem:set-dilaton` (ii));
Mumford's relation, [MR1728879] (only `rmk:set-qrrvertex` and the derivation of (QRR)); [MR2276766] Thm 2 (modern form of (M)); [GJR18]
(13), Guo-Ross 1703.06955 Props. 6.1-6.2, [LPfq] (17), (19)-(20), [ZZ08] (corroboration of `prop:set-cyclicLR` only); [CGYZ26] (remark to
`cor:mai-gap`); DNOPS Prop. 5.2, the uniqueness of homogeneous R-matrices, TRM's sign-of-`kappa` argument, the second half of DNOPS
Thm 8.1 (all dropped, G4 F7).

**3.9 The corollaries.** Corollary B (`cor:mai-gap`): Theorem A, (T), [PPZ15] Def. 2.13, Prop. 2.14, and proved statements only
(`lem:gau-confexp`, `cor:gau-rowsNN`, `thm:gau-R0` (iv), `thm:gam-D`, `prop:gau-symp`, `lem:set-M0sympl`, `lem:flh-frame`). `cor:mai-natfact`:
no input beyond those of Theorem A (`thm:gam-D`, `prop:set-cyclicLR`, `thm:gau-R0`, `prop:gau-symp`, `lem:gau-R1`, `lem:set-M0sympl`).
`thm:con-Pa`, `prop:con-deg`: no quoted input (appendix A, `lem:con-form`, the s03 lemmas named in G5).

**3.10 The conditionality sentence** (for `sec:int-inputs` and the closing paragraph of `thm:mai-main`; wording for the writer): Theorem A
holds modulo the quoted theorems (Q) of [Eynard:2007kz], [DuninBarkowski:2012bw] (as recalled in [Dunin-Barkowski:2015caa]), [ABDKS25] and
[PPZ15] with the statement numbers above, the value (QRR) of the R-matrix at the large-radius point, the genus-0 mirror theorem (M) and the
genus-0 axioms (G0) of the twisted theory, and Teleman's reconstruction (T); every other step is proved in the note. R7's `\draftnote` listing
open chain steps is deleted (there are none); if (QRR) is promoted (s.7 Q1), one `\draftnote` names O1 instead.

## 4. What is still open, by label

* **O1 `thm:set-qrr` (b)**, the derivation of (QRR) (G1 s.3.2): two classical inputs, not read in the full text by anyone: (i) Mumford's
  relation `c(E) c(E^vee) = 1`, i.e. `ch_(2k)(E) = 0` for `k >= 1` (C3b-2); (ii) for every `k >= 1`, `ch_(2k-1)(E) != 0` in `H^*(Mbar_(g,n))` for
  some `(g, n)` (C3b-4). Route: s.2.2. Outside the conditionality while (QRR) is quoted.
* **O2 the remarks** (R7, last entry; not assigned): every remark kept under D8 must carry its proof from the sources or say "not used; not
  proved here" (rule 10). Planned `\draftnote`s inside remarks, none used by a proof: `rmk:red-leaf` (O-9), `rmk:gam-blind` (the evaluation
  at `L = 0`), `rmk:gam-motivation` (GEN's lemmas), `rmk:mai-tau` (uniqueness of `tau`), `rmk:mai-degree` (only G5's a priori bound
  `floor((5n+4)/4)` is proved), `thm:mod-gap` route 1's remark on Hock s.3.2.
* **O3 bibliographic, not mathematical:** Givental 1996 IMRN numbering; Teleman journal numbering (existing `\draftnote`); Arbarello-Cornalba
  journal data (G1R C4-4, from memory); CGYZ26 re-read in the PDF (Defs. 1.1-1.2, Ex. 1.3, Rem. 1.4, Conj. 1.5, Prop.-Def. 4.5 (1); G5R F2.4);
  LLY97 only if added as a second source for (M) (not read).
* **O4 not claimed, not needed** (state so where relevant): injectivity of `ev` on `NN` (G5 s.2.6 item 3); `deg ph_n = floor(4n/5)` as an
  equality; `P_can(0, z)` not polynomial (VERIFIED to `z^18` only; Lemma E does not need it); regular semisimplicity of `A~` (true; not used).

## 5. Changes the design editor must make (OUTLINE files; s01-s03 here, s04-App. C in part 3)

**OUTLINE.md.**
* Rule 8: "(S1)-(S3)" -> "(M), (G0), (T) (and (IFT) if kept)" (R2); conditionality as in s.3.10.
* s.0.1: s02 uses "(M), (G0), (T), (QRR), (PPZ)"; "Keep visible": (QRR) enters s02 (rows), s03 (`lem:gau-colrows`) and `thm:con-glue`.
* s.0.2, chain table: C1 inputs (M), (G0), (T), (QRR), (PPZ); new row **C1'** `prop:set-cyclicLR` with its corollary (inputs (M), (G0));
  C6 adds `lem:gau-colmod`, `lem:gau-colrows`, `cor:gau-rowsNN` (inputs C1', (QRR) at `q = 0`); C9 inputs (DOSS via DNOPS s.2.6), (EO) Thms
  4.2, 4.6, (ABDKS), optionally (Hock).
* "Exact final list of inputs": replace by s.3.1-s.3.7; add the R7 table "steps of the chain left as `\draftnote`" with the entry "none".
* Main-theorem statement: unchanged by this phase (R6 is separate); keep "(alpha, beta) in Q^2" (D14), the scope of (N0).

**OUTLINE_2.md.** s.0.5: mark D2, D3, D5, D6, D7, D15 closed, with pointers (G1, G5, G4, G4, G2, G2). s.0.4: add the settled conflicts
"DNOPS s.2.3 dilaton-leaf sign" (PPZ's sign), "Givental Cor. 11.7" (pin Thm 11.6), "trivial frames: ten".

**OUTLINE_3.md (s01, s02).**
* s01 `sec:int-inputs`: the paragraph of s.3.10; R2 renaming; no list of open chain steps.
* `def:set-genus0`: `Delta_i = C_0^2 zeta^(3i)/L^2` and `D u_i = zeta^i L` become consequences (`cor` of `prop:set-cyclicLR`, G1 s.2.10 (a)).
* `thm:set-mirror`: (M) as in G1 s.1.1 with the pin of s.3.3 and the rewriting of s.2.1 (one paragraph: G1 s.1.3 items 1-5, s.1.4);
  delete the `\draftnote` (S1); state (G0) next to it; the Yukawa sentence stays evidence (appendix C).
* New **`prop:set-cyclicLR`** after `thm:set-mirror`: statement G1 s.2.0; Lemmas A-F (G1 s.2.2-s.2.7) as local lemmas; proof G1 s.2.8 with
  one sentence "only derivatives along the line `t = T H` are used" (G1R s.2.1); new corollary (G1 s.2.10 (a)-(b)): `Delta_i`, `D u_i`, and
  "`C_0 v_0` solves PF5 for every solution of the QDE". Local letters (G1R s.2.3): `s_a`, `g`, `W`, `w`, `c_k`, `phi_k`, `Mop`, `X_a`, `Fhat`,
  `ctilde_m`, `d_T`, `ell` -> declared letters (s.6).
* `def:set-rows`: define `Phi_i := C_0 Q^(i)_0/(zeta^i L)`, `A_n := [z^n](C_0 Q^(0)_0)`; the WKB form of colour `i` is the corollary above;
  both covariance clauses (`Q^(i)_(a,n) = zeta^(i(a+1-n)) Q^(0)_(a,n)`, `Phi_i(zhat,u) = Phi_0(zeta^-i zhat,u)`) cite `lem:gau-colrows` (s03),
  not CS s.1.
* `thm:set-qrr`: delete `\draftnote` (a) (proof G1 s.3.1). For (b): either (i) keep (QRR) quoted and print Steps 1-5 (G1 s.3.2, with
  C3b-1 and C3b-3 as in G1R s.4.2, s.4.4) as a remark "how (QRR) follows from [list]" with one `\draftnote` naming O1; or (ii) promote it
  (s.7 Q1). Drop "Prop. 8.5(ii) (s.8.2)" (ar5iv label); cite "8.7 Proposition (ii)" (v3).
* `rmk:set-qrrvertex`: subsumed by Steps 2-3 of (b) (keep one of the two).
* `def:set-graphsum`: replace the `\draftnote` (DA C6) by the naturality sentence (s.1 item (b)).
* `def:set-M0`: add **`lem:set-M0sympl`** (one proof: G4 s.5 or G1 s.5.1) and delete the `\draftnote`.
* `def:set-modcohft`: define `omega^m_v` by `Deltahat_i^(1/2) := zeta^(4i) varsigma` (s.2.10); "Hurwitz Frobenius manifold of DNOPS (8-1)"
  only as a remark, at `a = (0,0,0,0,4v)` (R3); add **`lem:set-modhom`** (G1 s.5.2 with the naturality sentence of G1R s.7); delete the
  `\draftnote`.
* `lem:set-dilaton`: statement with the domain hypothesis (`y = q + z 1` formal, `\varkappa` in `1 + m`; G1R C4-1, C4-2); proof G1 s.4 with
  PPZ15 Def. 0.3 (iii), Prop. 2.12, Def. 2.13 and Arbarello-Cornalba (1.7); one sentence that the Taylor-series form needs `\varkappa(0) = 1`
  (G3R K4); the genus-1 log term either with a pinned `1/24` or dropped (not used); delete the `\draftnote`.
* `def:set-confexp`: `NN` as presented in G5 s.2.2 (generators, `Lam` with `D Lam = 1`, `D i_2` from the concomitant); `Cf` with
  `Lam -> log(1-u)` and the automorphism remark of G5R F2.1; "rows ... lie in `NN` (FQ2 ...)" -> forward reference to `cor:gau-rowsNN` (R1 item 3).

**OUTLINE_4.md (s03).**
* `prop:gau-quintic`: add the step of s.2.9 (R1 item 4) with inputs `prop:set-cyclicLR`, `lem:gau-confexp`.
* `thm:gau-U`: the statement of G2 s.5.1 (a)-(d) over any field of characteristic 0, (e) as a remark under `prop:gau-equalprod`; proofs
  G2 s.3-s.6; delete the `\draftnote`s G1-G5, G7 and the one on `C[[u]]`; add **`lem:gau-lead`** before the uniqueness clause and
  **`lem:gau-degree`** before `prop:gau-graded`. Renames (G2R C5): the free constant `gamma`, `beta = u Bprod`, the field `k`, `S`, `V(A, A~)`.
  Optional: index shifts by `m = 1..5` (G2R R1).
* `lem:gau-comm`, `lem:gau-lift`, `lem:gau-scal`: statements and proofs G2 s.3, s.4.1, s.4.2.
* `prop:gau-graded`: step 1 by `lem:gau-degree`; one sentence on `zeta^(-gr)`-conjugation and averaging (G2 s.6.4); delete `\draftnote{G5}`.
* `lem:gau-diag`: scope (N0) printed; step 4 = G2 s.7 (general value `p_(0,1)(0)^5 prod rho^U_a`, value 1 under (N0)); the second proof
  from `lem:gau-lead` may replace steps 1-5; delete `\draftnote{G8}`.
* `thm:gau-G` Cor. (a): proof G2 s.8; delete `\draftnote{G9}`. `prop:gau-symp` step 3: cite `thm:gau-U` (b) for the pair `(A~, A~)`.
* `lem:gau-wkb`: every colour, over any field containing `Q(zeta)` (G2R C2); proof: colour 0 as planned, colour `i` by the substitution
  `z -> zeta^-i z` (s.2.5) or the transport operator `5 zeta^(4i) D`.
* `thm:gau-R0`: (iii) adds the input "WKB form of the actual colour-0 row" (corollary of `prop:set-cyclicLR`); (iv) `D_0` in `Q[[z]]`, `D_i`
  in `Q(zeta)[[z]]`, comparison ring of colour `i` `Q(zeta)[L, 1/L][[z]]` (G2 s.11 item 3); prove (iv) for colour 0 only and get the colours
  `i != 0` from the three-line computation (G2R C4); the last sentence: (QRR) enters s03 only through `lem:gau-colrows`.
* New **`cor:gau-rowsNN`** after `thm:gau-R0` (G1 s.2.10 (c)); new **`lem:gau-colmod`**, **`lem:gau-colrows`** (G2 s.9.3, s.9.4 with C3 and
  `psi_i = L Phi_i`) in `sec:gau-R0` after `lem:gau-wkb`.
* `cor:gau-colours`: proof G2 s.9.2; delete `\draftnote{G6}`.
