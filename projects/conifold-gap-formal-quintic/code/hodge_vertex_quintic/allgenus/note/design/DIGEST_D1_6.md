# Digest D1, part 6 (addendum): dependency graph, correction register, two further remarks

Addendum to `DIGEST_D1.md` .. `_5.md` (the piece map of `_4.md` lists `_5.md` as the last part; this file was added after it). Same
abbreviations and labels.

## A. Dependency graph of the statements of s.2 (arrows: "is used by")

* `lem:pol-HPF` -> `lem:pol-u0`, `lem:pol-poly`, `prop:pol-deg`, `thm:thp-a` (step 6), `lem:fin-N` (step 3).
* `lem:pol-poly` -> `thm:thp-a` (step 7), `thm:thp-b` (step 5), `thm:glu-main` (a), `lem:pol-Finf`.
* `lem:thp-form` -> `thm:thp-a`; `lem:thp-cfu` -> `lem:thp-red` -> `lem:thp-id`; [PropS] + [FSP] -> `lem:thp-S2` -> `lem:thp-id`.
* `lem:thp-relab` + `lem:thp-key` + `lem:thp-id` + [GI] Lemma V -> `prop:thp-res`; (`lem:thp-weight` optional, used by nothing).
* `thm:thp-a` + [CYC] + [U] + [GI] (Theorem G, Lemma Q) + `prop:thp-res` -> `thm:thp-b` -> `cor:thp-A2`, `cor:thp-frame`.
* [R0] ([U], [CYC]) + `lem:pol-poly` + `lem:glu-wkb` + [QRR] -> `thm:glu-main`; with `thm:thp-b` and [HF] -> `D_0 = c(z; 1)/h(z)`.
* `thm:glu-main` + [GAM] + `cor:glu-colours` -> `thm:glu-D` -> (F_g) for every `g` -> the Stage-2 chain (other digests).
* `thm:thp-a` -> [HF] (CP3 s.3.10 uses P(a) for polynomiality; P(a) does not use [HF]).
* [LMI0] + [HF] + `thm:thp-a` + `thm:thp-b` -> `prop:pol-Ginf`.
* Finite route: `lem:glu-wkb` + [U] + (R1) + [CYC] -> `lem:fin-N`; `lem:fin-P`; `lem:fin-N` + [QRR] -> `thm:fin-conn`; `lem:fin-N` +
  `lem:fin-P` + (R1) -> `prop:fin-pfree`; (R1) -> `prop:fin-row`; these feed only the evidence (s.5 E1-E5, E16).
* Not used by the all-orders chain: `lem:fin-P`, the formula for `K_q`, `prop:pol-deg`, `lem:pol-Finf`, `lem:pol-p5k`, `rmk:thp-Ai`
  (RF2 s.5; PB s.1 H9 with O13; K3).

## B. Correction register (where each correction lands in the digest)

| correction | content | lands in |
|---|---|---|
| PS-C1 | PB's FL3 is also used in (S2): cite FSL FL4 (s.4.5) | `lem:thp-S2` |
| PS-C2 | PB's residue ring needs `ord_(1/w) c_n -> +oo`, or use FSL's ring `L` | D1.15; `prop:thp-res` |
| PS-C3 | cite L1a and Lemma H for `x = sigma x'`, `z -> sigma^4 z'` and for `Y = Y_c eta`, `w -> tau w` | `lem:thp-S2`; `lem:thp-id` |
| PS-C4 | name the family; two sentences (finite hatted fibre; `H` is the hatted-order sum) | `lem:thp-relab` |
| PS-C5 / A6 | finite fibres give a formal sum; weight lemma, or the coefficientwise identity with `V` in `Q[z^-1][[u]]` | `lem:thp-key` |
| PS-C6 | Theorem G for `P^SP` is redundant (regular gauge, or multiply by `s(z)`) | `thm:thp-b` step 5 |
| PS-C7 | state `[u^0](Y^C_Q)_0 = 0` | `thm:thp-b` step 7 |
| PS-C8 | the ring of comparison `Q[L, 1/L][[z]]` in `Q((s))[[w]]`; strip the same exponential | `thm:glu-main`; D1.8; s.4.6 |
| PS-C9 | Lemma N's reading; add (CYC) to its inputs | `lem:fin-N`; `prop:glu-same` |
| PS-C10 | formal WKB uniqueness PROVED in one line | `lem:glu-wkb` |
| RS2 items 13, 27, 28 | (CYC)'s role in R0; normalisations refer to the fq_conifold family | [R0]; `thm:thp-b`; `prop:fin-row` |
| RS2 item 40 | (Gamma) blind to the twist measure | s.4.4 |
| SYN5b K3, P7 | Lemma N not needed for all orders; the glue is shorter | `thm:glu-main`; `thm:glu-D` |
| SYN5b K4, PB P1 | (CYC) is an input | `thm:thp-b`; `lem:fin-N`; [R0] |
| SYN5b K5 | the relabelling is complete | `lem:thp-relab` |
| SYN5b K8, PB P8 | regularity at `L = 1` is automatic given P(a); (Gamma) is the value | `rmk:glu-K8` below |
| SYN5b K11, K12 | FL dictionary; (S2) PROVED | s.4.12; `lem:thp-S2` |
| SYN5b O8, O10, O12, O13, O14 | lists of hypotheses; underclaim on FL3, FL6; "not used anywhere" rephrased; `Ginf` needs less | `thm:glu-D`; `thm:thp-b`; `thm:thp-a`; `prop:pol-Ginf` |
| Stage-4 O1, O2, O9, O13, O14 | "VERIFIED given (U)"; frames not independent; superseded test claim; proof/computation mixed; 28 and 7 independent identities | s.5 |
| Stage-4 O3, O4 | degree bound PROVED; interval `[-4n/5, 0]` | `prop:pol-deg` |
| Stage-4 O5, O12 | still a connection statement | `cor:thp-A2` |
| Stage-4 O6 | GAM's Theorem A: (i) without (U); (ii) given (U), Theorem G, (FL), with gaps G2-G3 closed by POLY's route | `rmk:thp-Ai`; use the POLY/PB route for `thm:thp-b` |
| Stage-4 O7 | not UR s.8's phase | `lem:thp-form` |
| Stage-4 O8, O11 | 1702.06096 also defines the formal quintic; 1803.01409's details unconfirmed | s.6 |
| Stage-4 O10 | "no ODE-formal argument" is true for the listed properties only | `rmk:glu-why` below |
| Stage-4 O15 | divisibility necessary, not sufficient | `prop:fin-row` |
| Stage-4 O17 | the model fixes nothing at `u = 0` | `rmk:nrm-model` |
| GR C1, C2, C6 | (H) not presupposed; the common factor `(2 pi z/5)^(1/2)`; different phase | `prop:pol-deg`; s.4.5; `lem:thp-form` |
| PR C1, C4, C6, C7, C8 | Lemma 3 PROVED by quotation; wording; the formal lemmas listed; divisibility; cite math/0108100 s.10 | `prop:pol-deg`; `cor:thp-A2`; `thm:thp-b`; `prop:fin-row`; s.6 |
| PR C3 | "c of `Phi^(0)` = `1/F_inf^QRR`" is equivalent to the main result, not a control | s.5 (do not list it as a control) |
| ER s.7 items 1-7 | as Stage-4 O1-O5; the QRR convention; still a connection problem | s.5; `prop:pol-deg`; `lem:pol-u0`; `thm:fin-conn` |
| RF2 C2, C3, C4 | the Jacobian belongs to the amplitude; `c_n` in `Q[L]` needs P(a), [HF] and the bound (or LMI at 0); the u = 1 convention is (QRR) | [HF]; `prop:pol-Ginf`; `thm:fin-conn` |

## C. Two further remarks the writers may want (sources and status)

* **`rmk:glu-why`** (rmk) *Why an exact input is needed.* Any two colour-0 WKB series differ by a constant factor `D(z)`; PF5, symplecticity
  (`D(z) D(-z) = 1` only), Z/5, (H) and the polynomial structure in `L` are all invariant under `D -> D exp(odd series in z^5)`. So no argument
  using only these properties decides `D`; Theorem P removes the conifold end by an exact identity (the residue at the pole), and (Gamma)
  supplies the input at `q = 0`. Src. GP s.8; GR s.2.8; Stage-4 O10 (correct for the properties listed; Theorem P and the chart are formal,
  their input is the integral). St. PROVED (as an invariance statement).
* **`rmk:glu-K8`** (rmk) *Regularity at `L = 1` is automatic.* By P(a) and [HF], `c_n(L) = L^(4n) ph^SP_n(L^-5)` is a Laurent polynomial, so
  "regular at `L = 1`" in (Gamma) is automatic and its content is the value `ph^SP_n(1) = h_n`; the chart proof of Theorem Gamma proves the
  regularity independently (it does not need polynomiality). Src. PB P8; SYN5b K8; PB s.3 last bullets. St. PROVED given P(a).

## D. What the note must contain from this scope (check-list for the writers)

1. D1.1-D1.14 definitions in the note's notation, with the clash list of s.4.13 resolved.
2. `defn:nrm-D`, `rmk:nrm-model`, `rmk:nrm-genus` (the last without the finite ranges, which go to the appendix).
3. `lem:pol-HPF`, `lem:pol-u0`, `lem:pol-poly`, `prop:pol-deg` (ER's proof; the ZZ route as a remark), `rmk:pol-points`; `prop:pol-Ginf`
   after Theorem Gamma.
4. `thm:thp-a`, `lem:thp-form`, `lem:thp-cfu`, `lem:thp-red`, `lem:thp-S2`, `lem:thp-id`, `lem:thp-relab`, `lem:thp-key`, `prop:thp-res`,
   `thm:thp-b` with the full hypothesis list, `cor:thp-A2`; optional `lem:thp-weight`, `rmk:thp-710`, `cor:thp-frame`.
5. `lem:glu-wkb`, `thm:glu-main` (with PS-C8's sentence), `cor:glu-Laurent`, `cor:glu-colours`, `thm:glu-D`, `rmk:glu-why`, `rmk:glu-K8`.
6. The finite route (`fin`) as a subsection or appendix, flagged as unused by the main proof; `prop:glu-same` there.
7. Evidence (s.5) only in the evidence appendix, with the independence caveats; no script, file or log names.
8. `\draftnote`s: (QRR) in the H9 form (s.7 item 2); ER s.4 step (5) if it cannot be expanded (s.7 item 6); any statement of `rmk:thp-Ai`
   (s.7 item 8).
