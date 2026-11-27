# Reduction of Conjecture N (part 3 of 3): deformation route, other candidate statements, Stage-2 routes, gaps

Continuation of `REDUCTION_DRAFT.md` (s.0-s.2) and `REDUCTION_DRAFT_2.md` (s.3: Lemma 1, Theorem R). Same labels and conventions.

## 4. The deformation route (task item 2)

**4.1 Set-up (PROVED).** Join the pulled-back model `A_0 = v' Amod(v)` to a system `A_1` (quintic, or any self-dual cyclic one-pole system
with the same product) by a path `A_s` of such systems with the same product and the same residue (e.g. interpolate the units `f, g` of
`b = (f, B/(f^2 g^2), f, g, g)`). Let `P_s` be the canonical gauges (`P_0 = 1`) and `X_s := (d_s P_s) P_s^-1`. Then
 (D0) `z X_s' = [A_s, X_s] + d_s A_s` (differentiate `z P_s' = A_s P_s - P_s At`);
 (D1) `X_s` is flat-regular, graded and infinitesimally symplectic, `X_s(-z)^T eta + eta X_s(z) = 0` (differentiate Prop. S);
 (D2) `X_s Y^C_{A_s} = d_s Y^C_{A_s}`, which has no positive power of `z` (differentiate Theorem G): **X_s annihilates the whole flag**
      `V_s^(k)`, `k >= 1`, by Lemma E; its translation keeps the base point on `Lambda_s` and moves `V_s^(1)` by `d_s T_{A_s}` (Lemma T).
 (D3) `d_s Omega^s = X_s.Omega^s` (left action), i.e. for the potentials at the base point
      `d_s F_g = (1/2) Sum_{g1+g2=g} B_s(dF_g1, dF_g2) + (1/2) B_s(d^2 F_(g-1)) + L_s F_g`, with the edge `B_s(psi, psi') = (X_s(psi) eta^-1 +
      eta^-1 X_s(psi')^T)/(psi + psi')` (the `epsilon`-derivative of `E_R` at `R = exp(epsilon X_s)`) and `L_s` the linear vector field
      `q -> -X_s q` (legs and translation).

**4.2 What must hold.** Suppose one only knew a NAIVE singular structure on all of `Hc_+`, `F_h = Phi_h(Vn(q)) + regular` with some non-linear
genus-0 function `Vn` (`dVn = V^(1)` at the base point). Inserting it into (D3) gives, besides regular terms,
`Phi_g'(Vn) [L_s Vn + B_s(dF_0, dVn)]` (the genus-0 transport of `Vn`; equals `d_s Phi_g(Vn)` if `Vn` is transported), terms containing
`B_s(dVn, .)` (zero by annihilation at the base point), and `(1/2) Phi'_(g-1)(Vn) B_s(Hess Vn)`, of order `Vn^(3-2g)`, which cannot cancel
against the transport (other power of `Vn`; universal constants `c_h`). So the plan's worry is correct for the naive form: one needs
`B_s(Hess Vn) = 0`. In fact the naive form itself is false: the model worker found singular two-point functions off `Lambda`
(`<t1(H) t2(H^3)>_h`, `<t2(H) t1(H^3)>_h`, `<t3(H) t0(H^3)>_h`) whose "Hessian" `(<XY>_h - F_h'' PiC(X) PiC(Y))/F_h'` depends on `h`, so no
genus-independent `Vn` exists (REFUTED, `../model/MODEL_SINGULAR.md` s.4). *Resolution (PROVED from (D2) and Lemma 0).* `B_s` lies in `Sym^2 T_v`, so only `Hess Vn` restricted to `T_v (x) T_v` at points
of `Lambda` enters, and there the singular part has no Hessian if ML_Lambda holds (`V^(1)` is linear). Example: for the string-invariant
completion `Vn = V^(1) + Sum_{n>=1} x^n V^(n+1)/n!` (`x` a genus-0 coordinate with `l_{-1} x = 1`), `Hess Vn = dx (.) dV^(2)` at the base
point, and `B_s(dx, dV^(2)) = 0` because `B_s` annihilates `V^(2)`. The string counterexample `<tau_1(H) tau_0(1)>_h = <H>_h` sits in a direction
(`tau_1(H)` alone) that is not in `T_v (+) Q phi_1`; inside `T_v` the frame `frame(2,0) = tau_2(1) + tau_1(H)` compensates it exactly
(`<frame(2,0) tau_0(1)>_h = <frame(1,0)>_h = 0`, string equation).
**Answers.** (i) The infinitesimal gauge equation does not force `B_s(Hess Vn) = 0` for a Hessian taken off `Lambda`; it forces
`B_s in Sym^2 T_v` (D2), and that suffices. Note `(P z^-1 J^C)_{>=0} = 0` already implies `(P z^-k J^C)_{>=0} = 0` for all `k`: annihilation
of `V^(1)` by the edges (the plan's "(G)") IS annihilation of the whole flag. (ii) Comparison with E3: annihilation alone (plus (Sp) and
flat regularity) suffices, at first order (rank 0 on ker Ann, genus 2, 3) and non-linearly (A1-A4 at genus 2, 3 and one dressing at genus 4,
all 378 graphs) - exactly as Theorem R predicts; flatness (D0) is never used. (iii) What the model must satisfy: ML_Lambda(h) for all `h`.
**4.3 Status of the route.** The deformation route is the infinitesimal form of Theorem R(c); it gives nothing beyond it, but one useful
corollary: **Proposition I.** ML_Lambda(h <= g) is an invariant of the gauge class: if it holds for one system of the shape (with its
canonical normalisation), it holds for all - apply Theorem R(c) to `P_{A'} P_A^-1`, which is flat-regular and satisfies (Sp) and (A)
(`P_A^-1 Y^C_A = Y^C`). Status: PROVED from the model to any system (Theorem R(c) as written); from a system to another system or back
to the model SKETCH: the proof of Theorem R uses the homogeneity of the source only for the bookkeeping of u-orders (regular monomials,
u-adic convergence of frame expansions), which for a non-homogeneous source must be replaced by the hatted u-filtration (NATURAL_FORM
s.6.1), uniform over the class. So ML_Lambda may be proved on any convenient representative; the model (homogeneous, TR allowed) is the
natural one.

## 5. Other candidate statements (task item 3)

**5.1 Toda / second difference (CONJECTURAL).** With `Z^(c=1)(x) = exp Sum_h hbar^(h-1) Phi_h(x)` (genus 0, 1 terms included),
`Phi = log` of a double-gamma (G-type) function of `x/hbar`, the second difference `Phi(x + hbar) + Phi(x - hbar) - 2 Phi(x)` is `log x` + const
(normalisations to be fixed). Statement (ML-Toda): on `Lambda_v` the second difference in the `phi_1` (= residue) direction of
`log A^mod` equals `log V^(1)` + regular. ML_Lambda implies ML-Toda; conversely ML-Toda plus homogeneity fixes the `Phi_h` up to the kernel of
the second difference (affine functions), so the Bernoulli constants `c_h` would FOLLOW. Natural source: the extended rationally 4-constrained
KP description (arXiv:2211.12259); shifting the residue `a_5` by `hbar` should be a dressing (gauge) transformation of the Lax operator
`L = d^4 + .. + phi d^-1 psi`, whose tau functions form a lattice obeying a discrete (Toda-type) bilinear equation; for `x + v/x` (rank 2) this
is the Toda lattice of the GUE.
**5.2 TR on the model and the x-y swap (allowed by the user's ruling).** Spectral curve `x = p^4 + v/p`, `y = p`, `omega_{0,2}` = the standard
bidifferential of P^1. (a) arXiv:1509.06954 (abstract read): TR on the superpotential curve of a semisimple point of a conformal Frobenius manifold
reproduces, under conditions, the CohFT at that point; its examples include Frobenius manifolds for hypermaps (search summary; to be
checked that `H_{0;1,4}` with form `dp` is covered). (b) arXiv:2212.00320 (abstract read): the universal x-y swap formula, and **an explicit
closed formula for the TR differentials of any spectral curve with unramified `y` and arbitrary rational `x`** - our curve exactly (`y = p`
is unramified, `x` is rational; the precise hypotheses of that theorem were not read). (c) arXiv:2312.16950 (search summary read): log-TR,
compatible with the x-y swap, efficient when the dual side is trivial - as here, since `X = p` has no ramification point. Expected
mechanism for ML_Lambda (CONJECTURAL, nothing checked): the closed formula expresses every `omega_{h,n}` through local data of `x` and
`y dx` at the poles of `x` (`p = 0`, encircled by the vanishing cycle, and `p = oo`); the `v`-singular terms come from `p = 0`, where the
deformations seen by `T_v`-insertions are invisible, so that on `Lambda` only the residue `v` survives, giving `Phi_h(V^(1))`.
**5.3 Two polarisations / Fock sheaf (PROVED reformulation).** In the language of arXiv:1411.7039, `P^` changes the opposite polarisation
`Hc_-` to `P^-1 Hc_-` with `Hc_+` fixed; Lemma E says the propagator between them vanishes exactly on `Hc_- intersect P^-1 Hc_-`, and (A) says
the negative half of the isotropic c = 1 sector `Sigma_C` (s.2.2) lies there. So the c = 1 part of the potential is polarisation
independent along `Sigma_C` (Lemma C), while the rest of `Hc_-` is not; ML_Lambda is the statement that the singular part lives on `Sigma_C`.
**5.4 Lie algebra of first-order operators.** The infinitesimal dressings that annihilate the flag, `a_C = {X : (X J^C)_{>0} = 0}` (flat-regular,
graded, infinitesimally symplectic), are not obviously closed under brackets (`[X, X'] J^C = X (X' J^C) - X' (X J^C)`, and `X' J^C` is an
arbitrary vector without positive powers, on which `X` creates positive powers in general), so "invariance under a Lie algebra" is not the
natural formulation; the right one is the groupoid of pairs (theory, flag) of Prop. I, on which the singular part is covariant: for
`X` in `a_C`, `X^` acts on `Phi_h(V^(1))|_Lambda` as the first-order operator `L_X` (transport of `V^(1)`), plus regular terms (Theorem R(c),
infinitesimally). Virasoro: the model is conformal, so its descendant potential satisfies the Virasoro constraints; the generators should
preserve `Sigma_C` (it is spanned by a homogeneous solution), and their restriction to it would be the c = 1 constraints (SKETCH, unchecked).

## 6. Ranked Stage-2 routes (task item 5)

What is to be proved is ML_Lambda(h) for all `h` (with `n = 0`: the gap `F_h = c_h v^(2-2h)`). By homogeneity it is a statement about the
model at one point; by Prop. I it may be proved on any representative of the gauge class.

| rank | route | difficulty | first concrete lemma |
|---|---|---|---|
| 1 | TR on `x = p^4 + v/p`, `y = p` + closed formula for unramified `y` (arXiv:2212.00320), identification arXiv:1509.06954 | moderate: literature exists; hypotheses and the extraction must be checked | L1: the model's ancestor correlators at `v`, in the basis `x^(c_a) dx` of LG_MODEL s.2, are the expansions of TR `omega_{h,n}` at the five branch points (normalisations fixed against `th_mla` data at `h <= 2`); L1': in the closed formula, the polar part in `v` of `omega_{h,n}` contracted with `T_v`-vectors vanishes, and `F_h = c_h v^(2-2h)` |
| 2 | combinatorics of the descendant towers: the `H^2`-tower is a hypermap generating function (arXiv:2211.12259 Thm 2), whose `v`-exponent is the number of vertices `>= 1` (model worker s.8, PROVED given that theorem) | moderate: needs a combinatorial model for the towers of `1`, `H^3`, `H^4` and for MIXED towers | L2: an analogue of arXiv:2211.12259 Thm 2 for the full slice `K` (all four towers and `tau_0(H)`), with the `v`-exponent counted by a non-negative statistic |
| 3 | genus-by-genus for `h = 0, 1`: genus-0 reconstruction (the only log of `F_0` is in the residue) and the genus-1 formula of the Frobenius manifold | moderate, but stops at `h = 1` | ML_Lambda(0): genus-0 ancestor correlators with a `T_v` insertion and negative exponent vanish (from the S-matrix and `c_555 = 1/a_5`) |
| 4 | extension across codimension 2 on the 5-dimensional Hurwitz space: singular only along `{a_5 = 0}`; at a generic point there two critical points meet the pole (local model `x + v/x`, GUE) | high: needs ML_Lambda for the GUE in families and extension across the caustic | the singular part at a generic point of `{a_5 = 0}` is the GUE c = 1 term in the local vanishing coordinate |
| 5 | KP / Toda (s.5.1): dressing shift of the residue, second-difference equation (would give the Bernoulli constants) | high | the dressing transformation shifting `a_5` by `hbar` maps the tau function of arXiv:2211.12259 to itself up to a gauge factor |
| 6 | another representative of the gauge class (Prop. I) with a known CohFT | unknown; no candidate | none |

Independent of all routes: (H) and (N) to all orders (connection problem between `q = 0` and the conifold; gauge worker s.6), Theorem U
(AI-audited only), and a human check of Lemma E and Theorem R.

## 7. What is not established (gaps, named)

* **ML_Lambda(h)**: VERIFIED on the lists of s.2.4 (genus 0-3, at most 5/5/4/3 insertions, psi-deficit `<= 6`, two primes) and on the
  model worker's lists (genus 0-4, exact); PROVED only on the `H^2`-tower (model worker, given arXiv:2211.12259 Thm 2); the gap
  (`n = 0`) VERIFIED for `h <= 5`. Corollary N_9 needs it for `h <= 9`, on the finite lists that occur in polar orders.
* **(H), (N)**: verified through `z^24`/`z^25`; not proved (gauge worker).
* **s.1.5 identification** of the project's hatted outer graph sum with the flat P-action graph sum: standard (frame changes and the exact
  dilaton factor), supported by the univ engine's regressions against the stored genus-3/4 tables; not re-derived line by line here.
* Lemma C uses the Fock-space formula of arXiv:math/0108100 (conventions not re-checked); Theorem R does not use Lemma C.
* Everything is AI-written and AI-checked; Lemma E, Lemma T and Theorem R have not been refereed. The tests of s.2.4 use the project's
  symbol `M_0` mod `zhat^7` through `th_core` (MVertex6), i.e. the same vertex code as the production runs; an independent engine (the
  model worker's exact `ms_engine`) agrees with it on 6182 flat correlators (`../model/ms_regress.log`, `ms_regress2.log`).
* Theorem R(c) (legs) is proved by the same argument as (a), (b), written in one line; Prop. I is PROVED only from the model outwards.
* The univ experiments quoted in s.3.4 were run by another worker and not re-run here; they are cited as evidence, not as input.

## 8. Files (this folder)

| file | content |
|---|---|
| `REDUCTION_DRAFT.md`, `_2.md`, `_3.md` | this note |
| `th_core.py` | model correlators in the flat frame (MVertex6 on the project symbol), the frames of `T_v` (inherited from the first attempt, reviewed: pairing, frames, index contraction) |
| `th_mla.py`, `th_mla_89.log`, `th_mla_61.log` | the ML_Lambda test of s.2.4 with controls C1, C2 |
