# Referee report on the theory worker's reduction (REDUCTION_DRAFT.md, _2.md, _3.md, th_core.py, th_mla.py)

Date 2026-11-13. Folder `allgenus/theory/referee/` (referee key `theory`). AI-written; no human has read it. Labels as in the worker's
note (PROVED / VERIFIED / CONJECTURAL / REFUTED); for the verdict on each claim: CONFIRMED (re-derived here by another route),
WEAKENED (true with an added hypothesis or a narrower scope than stated), REFUTED, UNVERIFIABLE (not checkable here). No file of the
worker or of the project was modified; all new files are in this folder. Primes: 2^89-1 (project), a third prime
p3 = 1208925819614629174706261 (least prime = 1 mod 10 above 2^80) and a 128-bit prime (least prime = 1 mod 10 above 2^127).

## 0. Verdict: sound, with corrections (no error in the reduction; the evidence for the Model Lemma is narrower than stated)

1. **The reduction is correct.** Lemma E (edge contraction), Corollary E' (edges in `T_v (x) T_v`), Lemma T (leaf point on the slice,
   `V^(1)(q_T) = T_P`, `s = d T_P - v = O(u^2)`), Lemma F (finiteness) and the proof of Theorem R(a), (b) were checked line by line and
   hold as written, conditional on ML_Lambda and on the frame/hatting identification of s.1.5. Prop. S for `r != rt` (Corollary U) is
   right. The deformation route, isotropy of `Sigma_C`, Lemma C and the Fock-sheaf reformulation are correct as reformulations.
2. **Independent confirmations (new routes, exact over Q where stated).** (a) Lemma E: random symplectic `R` (generic and Z/5-graded),
   random `f`, both slots, exact Q; non-symplectic and wrong-sign controls fail (`rf_lemE.log`). (b) Corollary E' and Lemma T in the
   FLAT frame on the actual canonical gauge (the worker and the gauge/model workers checked hatted tables mod p): `z^16`, `u < 16`, frames
   (0,0) and (1/3,-2/7), exact Q; the stored normalisation fails as it must (`rf_flat.log`). (c) **The genus-2 graph sum of
   `P.Omega^mod_{v(u)}` computed in the flat frame** (never computed before; the project computes the hatted sum): all six edge graphs
   regular, vertex `= (1/1200) t_c^-2 + regular`, and the seven per-graph values EQUAL the hatted values on record
   (`1/1200 - 7u/6000`; edge `[u^2]` = 479/756000, 1/2250, -127/2016000, 0, 0, 0): this confirms Theorem R at genus 2 directly in
   the form in which it is stated, and the s.1.5 identification graph by graph at genus 2; dressings violating (A) give singular
   edge graphs (`rf_g2flat.log`). (d) The ML_Lambda test re-derived (own frames at general v with homogeneity checked, own M_0 from
   the LG recursion fed to MVertex6, own contraction, third prime): identical counts 43/418/425/127, all vanish (`rf_mla_*.log`).
   (e) Annihilating NON-gauge dressings `exp(u^m z N) P` in the flat frame at genus 2: regular, vertex `c_2 T_P^-2` with `T_P != t_c`
   (Theorem R needs no flatness, as claimed). (f) ML_Lambda on the deep part of `T_v` (never reached by th_mla): 15 lists at genus 1-4, all vanish.
3. **Corrections / weakenings.**
   (i) "1013 singular lists" overstates what is tested: only **326** of them (5, 65, 168, 88 at genus 0-3) vanish by a cancellation
   between non-zero flat correlators, i.e. test the corrections that define `T_v`; in **642** (22, 328, 253, 39) every flat correlator
   of the expansion vanishes although it is dimensionally allowed (non-trivial vanishings of single flat correlators, blind to the
   frame corrections); **45** (16, 25, 4, 0) are trivial (dimension or the string equation alone). (ii) **No list tested by th_mla
   probes the part of `T_v` that depends on the `i >= 1` terms of the vanishing series** (`z^-5` and deeper corrections, first in
   `frame(6,0) = z^6 phi_0 + v z^5 H - (v^5/120) H`): every singular list containing such a frame is removed by the psi-deficit cut
   <= 6 at genus 2-4, and none exists at genus 0, 1 with <= 5 insertions (`rf_mla_scope.log`). The referee tested 15 such lists
   directly (s.3a: `<frame(6,0) H^5>_1`, 3 at genus 2, 6 at genus 3, 5 at genus 4, using the gap values for the pure
   `tau_0(H)^m` terms): all vanish, and the frames built from a wrong series do not. (iii) The table's "102 skipped", "262 skipped"
   include REGULAR lists; the singular lists skipped are 60 and 164. (iv) Corollary N / N_9 are stated as PROVED (conditional)
   without listing the s.1.5 identification (flat P-action sum = the project's hatted outer sum, graph by graph) among their
   hypotheses, although s.7 lists it as not re-derived; after (c) above it is VERIFIED at genus 2 only. (v) "each of the
   three hypotheses on P is needed in examples (s.3.4)": for (Sp) there is no example (without (Sp) the action is undefined).
4. **Open, as the worker says.** ML_Lambda has no proof for any h (the `H^2`-tower via arXiv:2211.12259 is the model worker's);
   (H), (N) only in a window; Theorem U AI-audited only; Prop. I only model-outwards.

## 1. What was read

Worker: REDUCTION_DRAFT.md, _2.md, _3.md, th_core.py, th_mla.py, th_mla_89.log, th_mla_61.log. Background: NATURAL_FORM.md s.0-s.3,
s.5, s.6; natural_form/theory/UNIFORM_REDUCTION.md; natural_form/lg/LG_MODEL.md; the plan; ../gauge/GAUGE_IDENTITY.md, _HAT.md;
../univ/UNIVERSALITY.md; ../model/MODEL_SINGULAR*.md, ms_ktest.py/.log; ../tr/TR_MODEL.md; ../xyswap/XYSWAP_MODEL.md. Abstracts of
arXiv:1509.06954 and arXiv:2212.00320 (WebFetch): the worker's descriptions of both match the abstracts.

## 2. Pen-and-paper check

**Lemma E (CONFIRMED).** With `alpha(x) = eta f(-x)`, the contraction is `Res_x alpha^T E_R(x, psi1) dx`. Expanding `1/(x + psi1) =
Sum_j (-psi1)^j x^(-j-1)` is legitimate in `Q((x^-1))[[psi1]]` because the numerator is divisible by `x + psi1` (this IS (Sp)); the
`eta^-1` term has no residue (`alpha x^(-j-1)` has powers `<= -2`); the other gives `-beta_{>=0}(-psi1)^T eta^-1 R^-1(psi1)^T`,
`beta = R^-1(x)^T alpha(x) = eta (R f)(-x)` by `R^-1(x) = eta^-1 R(-x)^T eta`. Signs, transposes and the substitution `x -> -psi1` are right.
**Corollary E' (CONFIRMED).** `(P z^-k J^C)_{>=0} = z^-k (P J^C)_{>=k}`; `z^-k J^C` is in `Hc_-` for `k >= 1` (`J^C` has a `z^0` term);
`E_P = E_R` for `P = P_0 R` because `P_0^T eta P_0 = eta`. Equivalence noted by the worker and correct: (A) <=> the edges annihilate
`V^(1)` <=> they annihilate the whole flag.
**Lemma T (CONFIRMED).** `V^(k)(q_T) = Omega(P z^-k J^C, -z phi_0) = -[z^(k-2)] Y'_0` (the worker's `(-1)^(k+1) [z^(k-2)] Y'_0(-z)` is the
same); `k >= 3` by (A), `k = 2` by the grading (`Y'_0` has only `z^-n`, `n = 1 mod 5`); `(P_0^vec)_00 = (P_0)_33 = d`, `d (P_0)_00 = 1` by
(Sp) at `z^0`; `[u^1 z^-1] Y'_0 = -v'(0) p_00(0)` gives `s = O(u^2)`. Lowered versus vector form: `eta = 5 J`, `J` the involution
`a <-> 3 - a`, so `P^vec = J P J` and (Sp) holds in either form; the worker's identification is harmless. New remark: `s` is frame
independent, `s = tau^(1/4) varpi_1 - v` (`varpi_0 t_c = varpi_1`); both frames give `s = u^2/30 + (127/3150) u^3 + ...`.
**Lemma F (CONFIRMED).** Translation insertions have psi-degree `>= 2` and cost one dimension each: at most `3h - 3 + n` of them.
**Theorem R (CONFIRMED as a conditional proof).** Points checked: every vertex of a connected graph with an edge has `n >= 1` and
genus `<= g - 1` (a genus-g vertex would need genus-0 leaves of a tree); genus-0 vertices with `<= 2` half-edges vanish by dimension
(`T = O(z^2)`); the frame expansion of an edge converges u-adically (corrections of `frame(k, a)` carry `v^((4(k-j) + c_a + 1)/5)`,
exponent `>= 1`) and agrees with the flat expansion because the double sum converges u-adically and only finitely many flat slots
are non-zero; for the vertex graph the sum over the number of leaves may be extended to all `m` (the added terms are zero) and
re-ordered u-adically, which justifies `Phi_g(v + s)`; `d^(2g-2) Phi_g(d T_P) = c_g T_P^(2-2g)`. Step (2) needs ML_Lambda(h) only for
lists with at least one `T_v` insertion, as stated. No gap found. The proof uses neither the u-ODE nor the normalisation beyond (A).
Part (c) (legs `w` of the target enter as `P^-1 w`, and `V^(k)(P^-1 w) = V'^(k)(w)`) is a correct one-line outline.
**The identification of s.1.5 (WEAKENED: argued, not derived by the worker; VERIFIED at genus 2 here).** A constant (in z) change of
frame of `H`, a rescaling `z -> lambda z` together with the degree grading, and the dilaton rescaling of the TFT are symmetries of
Givental's formula graph by graph, which is the worker's reason; the referee agrees but notes that Corollaries N, N_9 depend on it
and that only genus 2 has been checked (s.3: the flat per-graph values equal the hatted ones).
**Corollary U (CONFIRMED).** (R1) with `r != rt`: `r P_22(0,z) = rt P_11(0,z)`; Prop. S: `t(z) = P_22(0,-z) P_11(0,z) = (rt/r) c(z) c(-z)`;
with `r = -1/lambda^2`, `rt = -1`, `c = 1/lambda` gives `t = 1`; `T_P = T_A/lambda`.
**Corollary L (CONFIRMED given ML_Lambda(g)).** First-order term `c_g (2 - 2g) vhat^(1-2g)(t_c - vhat)` with `t_c - vhat = u^2/30` gives
`-(g-1) c_g u^(3-2g)/15`; at genus 2 confirmed numerically (s.3, control X1).
**Deformation route, Lemma C, s.5.3 (CONFIRMED as reformulations).** (D0)-(D3) follow by differentiation; the infinitesimal Lemma E
is `-(X f)_{>=0}`, zero on the flag by (D2). Lemma C: `[E(d,d)/2, L_f] = E(dL_f, d) = 0` by Lemma E. s.5.4 (bracket not closed) right.
**Isotropy of Sigma_C (CONFIRMED).** `eta^-1` is proportional to `eta`, so self-duality gives constancy in u; `(eta^-1)_11 = 0`.
**Consistency checks of s.2.1 (CONFIRMED).** `l_{-1} V^(k) = -V^(k+1)` (the `z^0` part of `q` drops out for `k >= 1`); `<frame(1,0)>_h = 0`
for `h >= 2`, `= 1/8` at `h = 1`; `<frame(1,0) H^m>_1 = 0` for `m >= 1` (dilaton + `<H^n>_1`), as ML_Lambda(1) requires.

## 3. Independent computations (all scripts in this folder, own code; project modules read-only as stated)

| check | route (different from the worker's) | result |
|---|---|---|
| Lemma E (`rf_lemE.py`) | exact Q; random symplectic `R = exp(X)` (generic and Z/5-graded, `z^8..z^10`), random `f` with 3-6 terms; both slots; controls: non-symplectic `R`, flipped sign | formula holds on every complete coefficient (6 cases, both slots); `E_R` divisible by `psi + psi1` iff (Sp); both controls fail. CONFIRMED |
| Corollary E', Lemma T, `d`, `s` (`rf_flat.py`, `rf_flat_core.py`) | FLAT frame, exact Q; canonical gauge from `ur_gauge` (read-only) normalised here; own `tau`, `v`, `Y^C` (model system checked on monomials); `V^(k)(q_T)` by direct pairing (the symplectic identity of Lemma T is NOT used) | frames (0,0) and (1/3,-2/7), `z^0..z^16`, `u < 16`: (Sp) exact; (G) on 575 complete coefficients; `T_P = t_c` through `u^14`; `V^(1)(q_T) = t_c`, `V^(2..8)(q_T) = 0` on all complete orders; `d (P_0)_00 = 1`, `d = tau^(1/4) varpi_0`; `s = u^2/30 + (127/3150) u^3 + ..` in BOTH frames; flat edge divisible by `psi + psi1` and its contraction with `V^(1..8)` zero on 2985 complete coefficients; control without `(-1)^n`: 92 slots non-zero. Stored normalisation: (G) fails (77), `V^(6)(q_T) != 0`, edge contraction fails (92), symplectic defect at `z^10`. CONFIRMED |
| Theorem R at genus 2, flat frame (`rf_g2flat.py`) | the seven genus-2 graphs of `P.Omega^mod_{v(u)}` in the FLAT frame: flat edges `E_P`, leaf point `q_T`, dilaton factor `d^(2h-2+n)`, model correlators `C v(u)^e` with `C` from MVertex6 fed with the referee's own LG `M_0` (128-bit prime, rationally reconstructed, re-checked mod 2^89-1); graph sums exact over Q | frame (0,0): vertex `u^-2(1/1200 - 7u/6000 + 287u^2/864000 ..)`; edge graphs `[u^0] = [u^1] = 0`, `[u^2]` = 479/756000 (loop on g=1), -127/2016000 (bridge g1-g1), 1/2250 (two loops on g=0), 0, 0, 0: EXACTLY the hatted values on record (UNIFORM_REDUCTION s.6.2). Frame (1/3,-2/7): edge graphs regular, vertex `= (1/1200) t_c^-2 + regular`. Controls: leaves removed: vertex `-1/900 u^-1` = the leafless value `-(3g-3)(4/9) c_g`, so the leaves supply `-1/18000 = -(g-1) c_g/15` (Corollary L); `(1 + u^m z E_21) P` (symplectic, flat-regular, violates (A)): edge graphs singular from `[u^m]`, m = 0, 1. CONFIRMED (genus 2) |
| annihilation without flatness (E3 analogue), flat frame (`rf_g2flat.py`) | `P -> exp(u^m z N) P`, `N = E_04 + E_43` (Pi-type, infinitesimally symplectic, annihilating; the dressed connection `exp(zN) A exp(-zN)` depends on z, so `P` is not a Theorem-U gauge to any cyclic system) | m = 0: all edge graphs regular; vertex `u^-2(1/1200 - u/3000 + ..) = (1/1200) T_P^-2 + regular` with `T_P = u + u^2/5 + ..` (NOT `t_c`); m = 1: regular, `T_P = t_c + O(u^3)`. Exactly the prediction of Theorem R(a), (b) for a non-gauge `P`. CONFIRMED (genus 2) |
| ML_Lambda test (`rf_mla.py`, `rf_mla_core.py`) | own frames at GENERAL v (homogeneity and `e_j in Z_{>=0}` asserted), own list enumeration, own contraction, own `M_0` (closed LG recursion; regressed against the LG_MODEL table) in MVertex6, third prime p3 | regression (`rf_mla_regress.log`): `F_2 = 1/1200`, `F_3 = 1/25200`, `<HHH>_0 = -5`, `<H>_1 = -1/12`, `<tau_1(1)>_1 = 5/24`, `<frame(1,0)>_1 = 1/8`, `<frame(3,3)>_2 = 3/3200`, `<frame(4,2)>_2 = -41/9600` (the worker's numbers), string identity. Test: genus 0/1/2/3 with <= 5/5/4/3 insertions: 43/418/425/127 singular lists, all vanish, same skip pattern. CONFIRMED (counts and values) |

**3a. Scope of the ML_Lambda evidence (`rf_mla_g*_n*.log`, second block of each log = refined statistics; `rf_mla_scope.log`;
`rf_deep.log`; `rf_deep2.log`).**

| genus, insertions | singular lists | cancellation between non-zero flat terms | every flat term 0, some dimensionally allowed | trivial (dimension / string) | C1 | C3 non-zero | C5 non-zero |
|---|---|---|---|---|---|---|---|
| 0, <= 5 | 43 | 5 | 22 | 16 | 5 | 4 of 10 | 4 of 9 |
| 1, <= 5 | 418 | 65 | 328 | 25 | 64 | 46 of 212 | 50 of 170 |
| 2, <= 4 | 425 | 168 | 253 | 4 | 165 | 123 of 274 | 120 of 209 |
| 3, <= 3 | 127 | 88 | 39 | 0 | 88 | 76 of 95 | 62 of 69 |

C3 = one correction coefficient of one frame doubled (the element leaves `T_v`); C5 = frames built with `Y^C(z)` in place of `Y^C(-z)`.
Both make the test fail on a large fraction of the affected lists, so the 326 cancellation lists are a genuine test of `T_v`.
**Deep part of `T_v`.** C4 (frames from a vanishing series with the sign of `(-1/(4z))^i` flipped, i.e. wrong from the `z^-5` terms on)
affects no list tested by th_mla: the frames that depend on the `i >= 1` terms are `frame(6,0)`, `(7,0)`, `(7,4)`, `(8,0)`, `(8,3)`,
`(8,4)`, ..., and every singular list containing one has a flat term of psi-deficit > 6. The referee evaluated those lists whose only
high-deficit flat terms are pure `tau_0(H)^m` correlators, using the gap `<H^m>_h = Phi_h^(m)(1)` for them (input: `F_h = c_h v^(2-2h)`,
VERIFIED `h <= 5`) and MVertex6 for the rest:

| genus, insertions | singular lists with a deep frame | evaluable | vanish | C4 (wrong series) non-zero |
|---|---|---|---|---|
| 1, <= 6 | 1 (`<frame(6,0) H^5>_1`: `<tau_6(1) H^5>_1 = -5/24`, `<tau_5(H) H^5>_1 = 7/24`, `(1/120) <H^6>_1 = 1/12`) | 1 | 1 | 1 (value 1/6) |
| 2, <= 5 | 10 | 3 | 3 | 3 |
| 3, <= 4 | 45 | 6 | 6 | 6 |
| 4, <= 2 | 12 | 5 | 5 | 5 |
| 4, <= 3 | 72 | 5 (the same five) | 5 | 5 |

So ML_Lambda also holds on the deep part where it can be reached (VERIFIED, p3, 15 lists: 1, 3, 6, 5 at genus 1-4), e.g. the genus-4 one-point function
`<frame(6,0)>_4 = 0`: the off-K remainders `+-19/100800000` of `<tau_5(H)>_4` and `<tau_6(1)>_4` (model worker, s.3 there) cancel.
**Extensions attempted and aborted** (MVertex6 too slow with five or six legs): genus 1 with six and genus 2 with five insertions
(`rf_mla_g1_n6.log`, `rf_mla_g2_n5.log`, marked ABORTED; no result). **Hessian data point of the refutation** (`rf_deep.log`):
`<tau_1(H) tau_2(H^3)>_h / F_h'(1) = -23/12` (h = 2) and `-13/48` (h = 3), as the model worker states (genus dependent).

## 4. Attempts to refute (what was tried, what happened)

* *Conventions tuned to pass?* No. `T_v` is fixed by `Y^C` alone (the referee rebuilt the frames at general `v` from the residue
  formula, with homogeneity asserted, and they agree with th_core at `v = 1`); the exponent `e = 2 - 2h + Sum w(k, a)` is the Euler
  grading (consistent with string, dilaton, `<HHH>_0 = -5/v`); the pairing sign is forced by `V^(k) = Omega(z^-k J^C, .)`, and the
  opposite sign (C5) fails. In the flat genus-2 sum nothing was adjusted: the per-graph values came out equal to the record.
* *A test that cannot fail?* The ML_Lambda test can fail (C3, C5); the flat graph sum can fail (X3, X4 singular; X1 changes the vertex
  by exactly the leaf term); Lemma E and the flat Lemma T/E' checks fail for non-symplectic `R`, for the stored normalisation and
  without the sign `(-1)^n`. But 687 of the 1013 singular lists involve no cancellation (45 of them are trivial), and C4 is vacuous
  on the worker's range (s.3a).
* *"Regular" confused with "zero"?* No: ML_Lambda demands ZERO for negative exponent (regular = monomial `C v^e`, `e >= 0`), and
  Theorem R concludes membership in `Q[[u]]`; both are what the proofs use. The genus-2 edge graphs have non-zero `[u^2]` (regular).
* *Low (g, n) presented as all-order?* The reduction (Theorem R) is all-genus and proved; ML_Lambda is labelled CONJECTURAL and its
  test ranges are stated. Overstatement only in "1013 lists" (s.3) and in Corollary N_9's "finite statement" (true, but the finite
  list was not written down; it includes the deep frames that th_mla never tested).
* *Non-gauge, annihilating P (the real content of "flatness is not used").* Tried in the flat frame (`exp(u^m z N) P`): regular, vertex
  `c_2 T_P^-2` with `T_P != t_c` (m = 0). Consistent with Theorem R; no counterexample found.
* *Quoted theorems.* arXiv:2212.00320 and arXiv:1509.06954 are cited only for CONJECTURAL routes; the abstracts say what the worker
  says ("unramified y, arbitrary rational x"; "under some conditions"); the worker states it did not read the hypotheses. The
  `H^2`-tower via arXiv:2211.12259 is the model worker's claim, correctly attributed. arXiv:1303.1043 (edge formula) and
  arXiv:math/0108100 (Fock-space formula, used only in Lemma C, which is not used) are standard.
* *Refutation of the naive singular structure (model worker's data).* The logic is right for the stated ansatz (on singular pairs
  the Hessian `(<XY>_h - F_h'' PiC(X) PiC(Y))/F_h'` must be genus independent). One data point re-computed in s.3a.

## 5. Hard rules

Worker: files only in `allgenus/theory/`; all ASCII; every file under 250 lines (140, 127, 116, 128, 91, 23, 31); no author names
(only terms such as Toda, Virasoro, Fock, Bernoulli, Hurwitz); papers cited by arXiv number; one process at a time (th_mla runs at
19:05, 19:08-19:12:37, 19:12:59, sequential); labels used throughout, with the weakenings of s.0 item 3. Referee: own folder and own
scratch only; no worker module imported (project modules `ur_gauge`, `fq_conifold` (through `ur_gauge`, frame patched in memory
only), `fq_fp.FpRing`, `fq3_mvertex.MVertex6/MConst` read-only); `py -I -B`, stdlib only; no run above 30 minutes (two extension runs were stopped by the referee before their time-out,
logs marked ABORTED). Process budget (2): exceeded briefly - one `rf_g2flat.py` run (about 5 s) and a few sub-second helper
processes (file edits, the ASCII check) ran while two background runs were active.

## 6. Files (this folder; all `py -I -B`, stdlib only)

| file | content | log |
|---|---|---|
| `rf_lemE.py` | Lemma E, exact Q, random symplectic `R`, controls | `rf_lemE.log` |
| `rf_flat_core.py`, `rf_flat.py` | flat-frame (G), Lemma T, `d`, `s`, Corollary E' on the canonical gauge, stored-normalisation control; args `K NU a b` | `rf_flat.log` (frames (0,0), (1/3,-2/7), `K = 16`, `NU = 16`; first block `K = 8`, `NU = 10`) |
| `rf_g2flat.py` | genus-2 graph sum of `P.Omega^mod_{v(u)}` in the flat frame, controls X1-X4, E3-type dressings | `rf_g2flat.log` |
| `rf_mla_core.py`, `rf_mla.py` | ML_Lambda test re-derived (own frames, own `M_0`, third prime), cancellation statistics, controls C1, C3, C4, C5 | `rf_mla_regress.log`, `rf_mla_g<h>_n<NT>.log` |
| `rf_mla_scope.py` | which singular lists the deficit cut removes; frames depending on the `i >= 1` terms of `Y^C` | `rf_mla_scope.log` |
| `rf_deep.py` | `<frame(6,0) H^5>_1` (the one deep list in reach) and its C4 control; the model worker's Hessian data point | `rf_deep.log` |
| `rf_deep2.py` | singular lists with a deep frame at genus 2-4, pure `tau_0(H)^m` terms from the gap, C4 control | `rf_deep2.log` |
