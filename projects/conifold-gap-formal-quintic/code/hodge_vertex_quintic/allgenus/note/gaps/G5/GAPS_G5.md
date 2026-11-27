# Gap closing, cluster G5: the connection problem and the main theorem (s08, s10)

Author-facing; not part of the note. Key `gaps/G5`. Pieces: this file (s.0 summary, s.1 item 1); `GAPS_G5_2.md` (s.2 item 2,
Corollary B and the constant of CGYZ26); `GAPS_G5_3.md` (s.3 item 3, s.4 item 4); `GAPS_G5_4.md` (s.5 evidence, s.6 suggestions
for the editor, s.7 dependencies on other clusters, s.8 sources read). Labels: PROVED (argument written here, every step), VERIFIED
(exact computation; range and field stated), CONJECTURAL, SKETCH, OPEN. Finite checks are evidence, never proof. Notation of
`NOTATION.md` (+ `_2`), labels of `LABELS.md` (+ `_2`); letters introduced here are local and marked so. Nothing outside
`note/gaps/G5/` was modified. AI-written; no human has read it.

## 0. Summary

| # | item | resolution | status |
|---|---|---|---|
| 1 | `thm:con-Pa` step 5: the domain of Theorem P(a) | **(b) PROVED** (s.1). The flat identity lives in `A_L := Q[L, L^-1, (L^5-1)^-1]`; the flat-hatted comparison is an identity in `Q((s))`, via the embedding `L -> 1/s`; HPF then holds in `Q[[u]][[zhat]]`. No real interval, no analytic continuation, no circularity; `prop:gam-hatflat` becomes part (iv) of P(a) | PROVED (formal); VERIFIED `n <= 9` at five points |
| 2 | `cor:mai-gap` (Corollary B) and the constant of [CGYZ26, Conj. 1.5] | **(b) PROVED given Theorem A**, with **(c) CITED** [PPZ15] Props. 2.4, 2.14, Def. 2.13 (read) and **(a) FOUND** a written proof of the regrouping in `natural_form/genus4_referee/review/REVIEW.md` s.2.3. Exact definition of the conifold expansion of `F_g`; frame scope `(alpha, beta)` in `Q^2`, extended to `C^2` and to every scaling of `varpi_0`, `varpi_1` (s.2.5). `c_g` IS the constant of [CGYZ26, Conj. 1.5], in the same normalisation, and the two conifold expansions agree on the period ring (s.2.6) | PROVED mod Theorem A and the inputs of `def:set-outer` |
| 3 | `cor:mai-natfact`: what is stated without proof | **(b) PROVED** (s.3): `Psi(u) := Qhat(u, zhat = 0)`; (a)-(c) follow from `thm:gau-R0` (iv), `D = 1`, (R1), `lem:gau-diag` (whose missing line G8 is written in s.3.3), symplecticity. OU3 s.3's relation `Psi eta^-1 Psi^T = eta^-1` is wrong in the note's convention; the correct one is `Psi^T eta^-1 Psi = 1` | PROVED; VERIFIED checks |
| 4 | `prop:con-deg` step 5: the triviality of the resonance | **(b) PROVED** (s.4): `P_j = (1 - x^5) o P~_j` for `j >= 1`; the right-hand side of the transport equation is `(1 - x^5) rho_n` with `rho_n` in `x Q[x]`, and the resonance coefficient is `[x^0] rho_n = [u^(4m)] r_(5m) = 0`. The two ingredients are in the sources (ER s.4's valuation lemma; the left factor `(1-u)` of POLY s.1, EVID s.2.1 (b)); their combination is written nowhere. It also gives existence for every choice of the free coefficients, without P(a) or Theorem U, and proves EVID's VERIFIED-only (b') | PROVED; VERIFIED `n <= 30` |

No item is left OPEN in this cluster. Items 2 and 3 are conditional on statements of other clusters (s.7), notably R1
(`prop:set-cyclicLR`, `cor:gau-rowsNN`, cluster G1) for the membership of the rows in the differential ring, and the symplecticity of
`M_0` (`def:set-M0`, R7's list).

## 1. Item 1: `thm:con-Pa` step 5 (the domain of Theorem P(a))

### 1.1 What was missing
DD1 s.3.3 step 5 (from PB s.2.2): "On `0 < u < 1` the analytic change of variables `sigma -> X` identifies the flat SP with the
hatted SP." The outline (`OUTLINE_7.md`, `thm:con-Pa`) asks to name the domain; it suggests `lem:fsp-param` (i), "as rational
functions" (CP3 s.3.10). Two problems. (1) "As rational functions" presupposes that `ph^SP_n` is rational, which is the conclusion
(iv) of P(a): circular inside P(a) (CP3 s.3.10 uses P(a) for polynomiality, DD1 s.7 item 7). (2) "On `0 < u < 1`" needs the hatted
`u`-series to converge on `(0,1)` and their sums to be the analytic stationary-phase coefficients there (true, `lem:fsp-param` (i)-(ii),
but never written), and it gives no ring for the HPF equation. Searching HVQ: no source names the ring for step 5. The device used
below is the one PS s.4.2 (d), s.4.5 uses for the glue ("`Q[L, 1/L][[z]]` embeds in `Q((s))[[w]]`"; "two rational functions that
agree as expansions at `L = oo` are equal"), so the resolution is in the spirit of the sources, but it is written nowhere for P(a).

### 1.2 Statement (proposed precise form of `thm:con-Pa`; colour 0, colour `i` by covariance)
Let `A_L := Q[L, L^-1, (L^5 - 1)^-1]` with the derivation `D`, `D L = (L^6 - L)/5`; it contains `q = (1 - L^-5)/5^5`, and `D = q d/dq`.
Let `c(z; L) = Sum_n c_n(L) z^n` be the full-Stirling saddle series of `def:gam-c`; by `prop:gam-rational` (`lem:fsp-univ`) every
`c_n(L)` lies in `A_L`. Let `iota : A_L -> Q((s))`, `L -> s^-1`. Write `PF5^L` for PF5 with `zD` replaced by `zD + L`.
* (i) (flat) `PF5^L (L c(z; L)) = 0` in `A_L[[z]]`.
* (ii) (comparison) In `Q((s))[[z]]`: `iota(c)(z) = Phi^SP(s^-4 z; s^5)`, i.e. `iota(c_n) = s^(-4n) ph^SP_n` for every `n`.
* (iii) (HPF) `Phi^SP` solves HPF in `Q[[u]][[zhat]]`.
* (iv) (polynomiality) every `ph^SP_n` is a polynomial in `u`; hence `c_n(L) = L^(4n) ph^SP_n(L^-5)` in `A_L`, so `c_n` lies in `Q[L, L^-1]`.
The last clause is the first statement of `prop:gam-hatflat`, which needs no separate proof.

### 1.3 Proof of (i): the flat equation, in `A_L`
Notation local to this proof: `S` the phase `SS_0` of `def:gam-c` written as a formal phase at `sigma = L`, i.e. by its Taylor data
`S^(k)(L)`, `k >= 2` (`S'(L) = 0`); `S''(sigma) = 5/sigma - Sum_j 1/(sigma - zeta^j) = -5/(sigma(sigma^5 - 1))`, so every `S^(k)(L)`
lies in `A_L`; `Delta := -S''(L) = 5/(L(L^5 - 1))`, a unit of `A_L`. For an amplitude `g(sigma, z)`, polynomial in `sigma` with
coefficients in `A_L[z]`, put `c^[g](z; L) := SPn_L[S; g A_St]/A_St(L, 0)` (`def:fsp-wick`), so `c = c^[1]`. The Taylor data at `L`
of `A_St/A_St(L,0)` lie in `A_L[[t]][[z]]` (binomial series for the half-density, `SS_2k` rational in `sigma`).

*Step 1 (the generic point).* Let `eps` be a formal parameter and `F(sigma; eps) := S(sigma) + eps sigma` (`eps` = the change of `lambda`;
`D = d/d lambda`, since `e^lambda = -q`). By `lem:fsp-ift` (over `A_L[[eps]]`) there is a unique critical point
`c(eps) in L + eps A_L[[eps]]`; differentiating `S'(c(eps)) + eps = 0` gives `c' = -1/S''(c) = c (c^5 - 1)/5`. Let
`psi_eps : A_L -> A_L[[eps]]` be the ring homomorphism `L -> c(eps)` (well defined: `c`, `c^5 - 1` are units). The map
`exp(eps D) = Sum_m eps^m D^m/m!` is also a ring homomorphism (Leibniz), and `e(eps) := exp(eps D)(L)` satisfies
`e' = exp(eps D)(D L) = e (e^5 - 1)/5`, `e(0) = L`. A formal ODE `y' = y (y^5 - 1)/5` has exactly one solution with `y(0) = L`
(its coefficients are fixed recursively), so `c(eps) = e(eps)` and `psi_eps = exp(eps D)` on `A_L`.
*Step 2 (`D` on saddle series).* The Taylor data at `c(eps)` of `F(.; eps)` (orders `>= 2`) and of the normalised amplitude
`g A_St/A_St(c(eps), 0)` are the images under `psi_eps` of the Taylor data at `L` of `S` and of `g A_St/A_St(L, 0)`: they are values at
`c(eps)` of rational functions of `sigma` (and of the binomial series of the half-density), and the Taylor map of a rational function
commutes with composition. By `lem:fsp-univ` (U1), `SPn_(c(eps))[F(eps); g A_St]/A_St(c(eps), 0) = psi_eps(c^[g])`. Put
`r := A_St(L,0)^2 = 5L/(L^5 - 1)` (in `A_L`); then `A_St(c(eps),0)/A_St(L,0) = (psi_eps(r)/r)^(1/2)` (binomial series). Apply `lem:fsp-dq` (FL4)
with base ring `A_L` at `eps = 0`, to the amplitude `g A_St/A_St(L,0)`, which is free of `eps`:
`d/d eps|_0 [ e^((F(c(eps);eps) - F(L;0))/z) (Delta(eps)/Delta)^(-1/2) SPn_(c(eps))[F(eps); g A_St/A_St(L,0)] ] = SPn_L[S; (sigma/z) g A_St]/A_St(L,0)`.
The right side is `(1/z) c^[sigma g]`. On the left, the bracket is `e^((F(c(eps);eps) - F(L;0))/z) (psi_eps(L)/L) psi_eps(c^[g])`, because
`(Delta(eps)/Delta)^(-1/2) (psi_eps(r)/r)^(1/2) = (psi_eps(r/Delta)/(r/Delta))^(1/2)`, `r/Delta = L^2`, and a series with constant term 1 has
one square root with constant term 1; moreover `d/d eps F(c(eps); eps) = c(eps)` (criticality) and `d/d eps|_0 psi_eps = D` (Step 1).
Hence `(L/z) c^[g] + (D L/L) c^[g] + D c^[g] = (1/z) c^[sigma g]`, i.e.
    `(zD + L)(L c^[g]) = L c^[sigma g]`   in `A_L((z))`.                                                     (1.1)
*Step 3 (PF5).* `PF5 = (zD)^5 - 1 - q prod_(k=1..5)(5 zD + kz)` with `q` acting by multiplication on the left. Iterating (1.1) (the
amplitudes `sigma^m g` stay free of `eps`; this is the iterated clause of `lem:fsp-dq`, RL-C2) and using `A_L`-linearity of `SPn`:
`PF5^L (L c) = L ( c^[sigma^5 - 1] - c^[q prod_k (5 sigma + kz)] )`.
*Step 4 (the shift).* Work with the Mellin-Barnes form `q^s Gamma(1+5s)/prod_j Gamma(1+s-zeta^j/z)` of POLY s.5-s.6: its phase has the
same Taylor data `S^(k)(L)`, `k >= 2`, and its normalised amplitude the same Taylor data as those of `def:gam-c` (`lem:con-form` (d);
the half-densities differ by a constant phase), so it has the same `c^[g]`.
With `a := (sigma^5 - 1) A_St`, the translate of `a` by `p = z` (`lem:fsp-transl` (b), over `A_L`) is
`b_z(t) = exp((S(L+t+z) - S(L+t))/z) a(L+t+z)`. In `A_L[[t]][[z]]`: `exp(S'(L+t)) = q (5 sigma)^5/(sigma^5 - 1)` at `sigma = L + t`
(both sides are 1 at `t = 0`, since `q (5L)^5 = L^5 - 1`, and have the logarithmic derivative `S''`), and the formal Stirling
functional equation `lem:fsp-stirfe` for `Gamma(1 + 5s)` (five shifts) and for each `1/Gamma(1 + s - zeta^j/z)` gives
`b_z = q prod_(k=1..5)(5 sigma + kz) A_St`. This is DD1 s.3.3 steps 2-3 [W] (POLY s.6 (a); PB s.2.2; FSL4 s.9.2), now read in
`A_L[[t]][[z]]` (the arguments `1/y = z/(5 sigma)`, `z/(sigma - zeta^j)` lie in `z A_L[[t]]`). By `lem:fsp-transl` (b),
`SPn_L[S; b_z] = SPn_L[S; a]`, so the bracket of Step 3 vanishes. QED (i).
*Alternative for Step 1-2 (holomorphic).* For every `L_0` in `CC_*` the coefficients are holomorphic near `L_0` (`lem:fsp-param` (i)) and
`lem:fsp-dq` holds holomorphically; each coefficient of (1.1) is then an identity of rational functions on an open set, hence in `A_L`.
Either version is complete; the formal one keeps the note free of analysis.

### 1.4 Proof of (ii): the comparison, in `Q((s))`
*Step 1.* `iota` is injective (`Q(L) -> Q((s))` is a field embedding) and `D`-equivariant: on `Q((s))`, `D = -(1-u) d/du =
-((1 - s^5)/(5 s^4)) d/ds`, and `D(s^-1) = (1 - s^5)/(5 s^6) = iota((L^6 - L)/5)`; a derivation-compatible ring map on generators is so
on the localisation.
*Step 2 (change of variables).* Over `A := Q((s))` (a field), the formal map `phi(X) := L/X`, `sigma - L = L Sum_(m>=1) (1 - X)^m`, has
`phi(1) = L`, `J := phi'(1) = -L`, a unit. `lem:fsp-cov` (L1a) over `A` gives `SPn_L[S; a] = SPn_(X=1)[S o phi; (a o phi) phi'/J]`.
*Step 3 (the phases agree).* In `A[[X - 1]]`: `S(phi(X)) - S(L) = -s^4 (WWhat(X; u) - WWhat(1; u))` (`lem:con-form` (c)). Proof: both
sides vanish at `X = 1`; their `X`-derivatives are `-(L/X^2) S'(L/X)` and `(L/X^2) log((1 - uX^5)/(1 - u))` (from `lem:con-form` (a),
(c): `WWhat'(X) = Sum_m u^(m-1)(X^(5m-2) - X^-2)/m`); so it suffices that `S'(L/X) = -log((1 - uX^5)/(1 - u))`. Both vanish at `X = 1`,
and both `X`-derivatives are images of the rational function `5 u X^4/(1 - u X^5)` (`u = s^5`): on the left `S''(L/X) (-L/X^2)` with
`S'' = -5/(sigma(sigma^5 - 1))`, on the right by differentiation. Two ring homomorphisms `Q[s^(+-1), X^(+-1), (1 - s^5 X^5)^-1] -> A[[X-1]]`
(expand at `X = 1` first, then in `s`; or expand `(1 - s^5 X^5)^-1 = Sum_P u^P X^(5P)` in `Q[X][[u]]`, then at `X = 1`) agree on the
generators, hence everywhere; so the two derivatives are equal. Characteristic 0: a series in `X - 1` is fixed by its derivative and
its value at `X = 1`.
*Step 4 (the amplitudes agree).* `(A_St o phi) phi'/(J A_St(L,0)) = frak ahat/frak ahat(1, 0)` in `A[[X-1]][[z]]`, `z = s^4 zhat`. The
half-density: at `sigma = L/X` (`L^5 = 1/u`), `5 sigma/(sigma^5 - 1) = 5 L u X^4/(1 - uX^5)` and `phi'/J = 1/X^2`, so the normalised left half-density
is the square root with constant term 1 of `(1 - u)/(1 - uX^5)`, i.e. `(1 - uX^5)^(-1/2)/(1 - u)^(-1/2)`, that of `frak ahat`. The Stirling parts: `z^(2k-1)(5 sigma)^(1-2k) =
5^(1-2k)(u X zhat)^(2k-1)` and `z^(2k-1) Sum_j (sigma - zeta^j)^(1-2k) = (u X zhat)^(2k-1) Sum_j (1 - zeta^j s X)^(1-2k) =
(u X zhat)^(2k-1) 5 Sum_P binom(2k-2+5P, 5P) u^P X^(5P)`, the exponent of `frak ahat` (`lem:con-form` (b)-(c)); equal by the same
two-homomorphism argument.
*Step 5 (the rescaling).* Write the amplitude of Step 4 as `a(X; z) := frak ahat(X; z/s^4)/frak ahat(1, 0)`. `lem:fsp-homog` (Lemma H) with
the unit `mu = s^4` of `A`: `SPn_(X=1)[s^4 (-WWhat); a(.; z)](z) = SPn_(X=1)[-WWhat; a(.; s^4 zhat)](zhat)` at `zhat = z/s^4`, and
`a(.; s^4 zhat) = frak ahat(.; zhat)/frak ahat(1,0)`. With Steps 2-4 and the normalisations (`c_0 = ph_0 = 1`):
`iota(c)(z) = Phi^SP(z/s^4; u)` in `A[[z]]`, i.e. `iota(c_n) = s^(-4n) ph^SP_n`. The hatted data lie in `Q[[u]]` (`[t^k] WWhat(1+t) =
Sum_m u^(m-1)[(-1)^k + binom(5m-1,k)/(5m-1)]/m`, `WWhat''(1) = 5/(1-u)`, `frak ahat` in `Q[X, zhat][[u]]`), so `Phi^SP` is computed in
`Q[[u]][[zhat]]` and, by (U1) for `Q[[u]] -> Q((s))`, equals the same functional over `A`. QED (ii).

### 1.5 Proof of (iii) and (iv)
(iii) Apply `iota` to (i): `PF5^L(s^-1 iota(c)) = 0` in `Q((s))[[z]]` (Step 1 of 1.4). The substitution `z = s^4 zhat` is a ring
isomorphism `Q((s))[[z]] -> Q((s))[[zhat]]` (`s^4` is a unit), under which `D` (at fixed `z`) becomes `D + (D zhat) d/dzhat`,
`D zhat = (4/5) zhat (1 - u)/u`. The chain-rule identity `(zD + L)(s^e F) = s^(e-1) A_e F` of `lem:gau-HPF` (POLY s.1) holds in
`Q((s))[[zhat]]`; with `e = -1` it turns `PF5^L(s^-1 Phi^SP)` into `s^-1 HPF(Phi^SP)` (POLY s.1; DD1 `lem:pol-HPF`). So `HPF(Phi^SP) = 0`
in `Q((s))[[zhat]]`; HPF has coefficients in `Q[u, zhat]` and `Phi^SP` lies in `Q[[u]][[zhat]]`, so the identity holds there.
(iv) `lem:gau-poly` (`ph^SP_0 = 1`). Then `iota(L^(4n) ph^SP_n(L^-5)) = s^(-4n) ph^SP_n(s^5) = iota(c_n)`, and `iota` is injective. QED.

### 1.6 For the writer of s08
* Replace the `\draftnote` of `thm:con-Pa` by the statement 1.2 and the proofs 1.3-1.5. Inputs: `lem:fsp-univ`, `lem:fsp-ift`,
  `lem:fsp-dq` (iterated clause), `lem:fsp-transl` (b), `lem:fsp-stirfe`, `lem:fsp-cov`, `lem:fsp-homog`, `lem:con-form` (a)-(c),
  `prop:gam-rational`, `lem:gau-HPF`, `lem:gau-poly`. Not used: Theorem U, (QRR), (CYC), `prop:con-deg`, `thm:gam-main`.
* `prop:gam-hatflat` then cites `thm:con-Pa` (iv) for `c_n(L) = L^(4n) ph^SP_n(L^-5)`; its separate proof (CP3 s.3.10 with RF2 C2) may
  stay as a remark. The order `thm:con-Pa` -> `prop:gam-hatflat` has no circularity: (ii) does not use (iv).
* Step 1 of 1.3 is a small formal lemma (FL4 at the generic point); the editor may place it in appendix A after `lem:fsp-dq`
  (suggested local name only; no label is registered by this cluster).
* Evidence (appendix C only): `g5_flat_hatted.py` (s.5): `c_n(L_0) = L_0^(4n) ph^SP_n(L_0^-5)` exactly for `n <= 9` at
  `L_0 = 2, 3, -2, 1/2, 5/3`, i.e. also at `u = 32` and `u = -1/32`, outside `0 < u < 1`.
