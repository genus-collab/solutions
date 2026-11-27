# Stage 4, Gamma route: the connection problem D_0 = 1 (2026-11-14)

Folder `allgenus/stage4/gamma/` (task key `gamma`). New files only; project modules are imported read-only (named in s.7). AI-written and
AI-checked only: no human and no referee has read this. Labels: **PROVED** (argument written here or in `GAMMA_PROOFS.md`, or a quoted
theorem), **VERIFIED** (exact computation, range stated), **CONJECTURAL**, **SKETCH**, **REFUTED**. Conventions: `UNIFORM_REDUCTION.md`
(UR) s.1-s.5, `allgenus/gauge/GAUGE_IDENTITY.md` (GAUGE) s.1-s.6, `natural_form/referee/ref_rows.py` (hatted WKB). `u = 1 - 5^5 q`,
`L = u^(-1/5) = (1 - 5^5 q)^(-1/5)`, `s = 1/L`, `z = s^4 zhat`, `zeta = e^(2 pi i/5)`, colour 0 unless stated (all others by Z/5).
Proofs: `GAMMA_PROOFS.md`. Primes p61 = 2^61-1, p89 = 2^89-1.

## 0. Bottom line

1. **The precise reduction (s.1, PROVED).** The colour-0 WKB coefficients `R_n(L) := L^(4n) ph_n(L^-5)` are polynomials in `L`
   (arXiv:0710.0889 Thm 4(ii)). `ph_n(0) = f_n` (the model's stationary-phase coefficients) holds AUTOMATICALLY for every normalisation, so
   the model fixes nothing at `u = 0`. The free constant at `n = 5m` sits at `u^(4m)`; it equals `R_5m(L = 0)`, the value at **q = infinity**.
   (H) is the identity `ph_n(1) = h_n` for `5 not | n` (h_auto); (N) says that the free constants equal their canonical values `kappa_m`.
2. **A conifold-natural integral representation (s.2).** Replace every Gamma factor of the Mellin-Barnes integrand of the formal quintic,
   `Gamma(5s+1) prod_j Gamma(zeta^j/z - s) (-q)^s`, by its full Stirling series. In `x = 1/(z s)` this is ONE 1-form
   `Omega = exp(W(x;u)/z) A(x;z) dx` whose phase `W = log(1-u)/x + Int_0^x log(1-t^5) t^-2 dt = -(u/x + x^4/4) + ...` is a deformation of the
   model `-(x^4/4 + v/x)`, with five saddles `x^5 = u` and a simple pole at `x = 0` (the oscillatory-integral picture of UR s.8, now derived).
3. **Theorem A (PROVED given Theorem U; s.3).** (i) `Res_{x=0}` of the rows of `Omega` is EXACTLY the quintic's vanishing solution `Y^C_Q`;
   (ii) the formal stationary-phase expansion of `Omega` ("full-Stirling WKB") is `P_can Q^m`: **the canonical normalisation is the
   full-Stirling Mellin-Barnes normalisation.** Hence `D(z) = c(z;L)/R(z;L)` for any `L`, `c` = full-Stirling saddle series (colour 0).
4. **Identity (Gamma) (s.4).** `c(z; L) = R(z; L)`: the QRR-normalised WKB of the formal quintic IS the steepest-descent series of its own
   Mellin-Barnes integral with every Gamma Stirling-expanded. Equivalently `c(z; 1) = h(z)`. By Theorem A: **(H) + (N) to all orders <=>
   (Gamma).** VERIFIED: exactly over Q for `n <= 30` (three values of `L`), mod p61 and p89 for `n <= 200` (two values of `L`), with
   controls. Proof: SKETCH (a single-Gamma connection at `q = 0`; two routes, s.4.3).
5. **Consequence.** `D = 1 + O(z^201)`: `a_1 = a_3 = ... = a_39 = 0` (E4 had `z^24`). So **(F_g) holds for every g <= 67**, modulo Theorem U
   (VERIFIED computation + Theorem A). For all genera (F) is reduced to the identity (Gamma), which is purely about Gamma functions.
6. **Closed forms (VERIFIED / PROVED).** Free constants: `Sum_m R_5m(0) z^(5m) = Ginf(z) := exp(Sum_l (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1)))`,
   the Stirling factor of `prod_j Gamma(zeta^j/z + 1/5)` (q = infinity Gamma normalisation; VERIFIED m <= 40; called `g(z)` in the
   `gm_compare` logs). Vanishing period (PROVED from Theorem A(i), VERIFIED to u^26, two frames):
   `varpi_0 t_c = -Sum_(k>=1) alpha_(k-1) (log(1-u))^k/k!`, `alpha_m = [y^m] exp(Sum_k b_k (5^(1-2k) - 5) y^(2k-1))` (Stirling factor of
   `Gamma(5w)/Gamma(w)^5` at `w = 1/y`; the sum is its integrated Borel transform).

| point | `L` | WKB normalisation | Gamma factors Stirling-expanded | status |
|---|---|---|---|---|
| q = 0 | 1 | `R(1) = h` (QRR) | `Gamma(1+5/z) prod_(j!=0) Gamma((zeta^j - 1)/z)` | definition |
| q = infinity | 0 | `R(0) = Ginf` | `prod_j Gamma(zeta^j/z + 1/5)` | VERIFIED m <= 40 |
| conifold | infinity | top coefficient `f_n`; canonical gauge | all six (saddle of `Omega`) | Theorem A (PROVED mod U) |

## 1. The WKB polynomials and what the model fixes

**Set-up (ref_rows).** `psi_0 = L e^(u_0/z) Phi(zhat, u)`, `Phi = Sum_n zhat^n ph_n(u)`, `ph_0 = 1`; the `zhat^(n+1)` coefficient of (HPF) is
`(1-u)[(4n - 5 theta_u) ph_n + rho_n(u)] = 0` (`rho_n` from `ph_0..ph_(n-1)`; the factor `1-u` is manifest, `gm_qrr.py`). So `ph_n` is a
polynomial, unique for `5 not | n`; for `n = 5m` the coefficient of `u^(4m)` is free (the residue `[u^(4m)] rho_n = 0` is needed and holds).

**Flat form.** `R_n(L) := L^(4n) ph_n(L^-5) = Sum_k ph_(n,k) L^(4n-5k)` is the coefficient of `z^n` in `psi_0/(L e^(u_0/z))`.
* (PROVED) `R_n in Q[L]`, i.e. `deg ph_n <= floor(4n/5)`: our `psi_0/(L e^(u_0/z))` is `h(z)` times the series `Sum_s Phi_s z^s/L` of
  arXiv:0710.0889 (5) for `n = 5` (their `F(w, x)` is our colour-0 Frobenius series with `w = 1/z`, `x = q`), and Thm 4(ii) there gives
  `Phi_s in L Q[L]`. VERIFIED here for `n <= 30` (degrees `0,0,0,2,3,4,4,5,..`; `= floor(4n/5)` except `n = 2`).
* (PROVED) `ph_n(0) = f_n` for EVERY choice of the free constants: the `u^0` part of (HPF) is the model's recursion (pivot `4n`), and the
  free constants enter only at `u^(4m)`, `m >= 1`. So "ph_n(0) as fixed by the model" is not a condition; UR s.2 used it via (H).
* The free constant `x_m := [u^(4m)] ph_5m = [L^0] R_5m = R_5m(0)`: the colour-0 WKB at `q = infinity` (`L = 0`) is the scalar series
  `x(z) = Sum_m x_m z^(5m)`. Rescaling the WKB by `(1 + d z^5)` shifts `x_1` by `d` and nothing below `z^5` (control C in `gm_compare`).

**The reduction (PROVED).** Let `kappa(z) = Sum_m kappa_m z^(5m)` be the `L^0`-values of the canonically normalised colour-0 WKB
(`Q D = P_can Q^m` with `D = 1`). Then
* (H) holds through `z^N` iff `ph_n(1) = h_n` for all `5 not | n <= N` with the polynomial solution (VERIFIED `n <= 200`, p61, p89);
* given (H), (N) holds through `z^N` iff **(PC_m)** `x_m = kappa_m` for `5m <= N`, where `x_m = h_5m - Sum_(k != 4m) ph_(5m,k)` is fixed by `q = 0`.

`kappa_m` is identified in s.3: `kappa_m = [L^0] c_5m(L)` (full-Stirling saddle series). Empirically `x_m = kappa_m = [z^(5m)] Ginf`.
**Example n = 5** (exact): `ph_5 = 47735919/128000000 - (203511/250000) u + (548547/1000000) u^2 - (6789/62500) u^3 + (4537/3937500) u^4`;
`ph_5(0) = f_5` (model), `ph_5(1) = h_5 = -257071/8064000000` (QRR), free constant `x_1 = 4537/3937500 = B_6(1/5)/6 = [z^5] Ginf`.
**Checks.** E4 (GAUGE s.6) is (PC_m) for `m <= 4` plus (H) to `z^24`; our runs give `m <= 40`. The stored normalisation
`P_stored = P_can (1 + m5 z^5)`, `m5 = -43397/807975000`, has free constant `x_1 + m5 = 4437977/4039875000`, whereas
`kappa_1 = x_1 = 4537/3937500 = B_6(1/5)/6`: the stored normalisation is excluded at `z^5`, as E4 found.

## 2. The full-Stirling Mellin-Barnes 1-form

**MB integrand.** `a(s) = Gamma(5s+1) prod_(j=0..4) Gamma(zeta^j/z - s)`, `sigma = z s`: `(1 - sigma^5) a(s) = prod_(k=0..4)(5 sigma - k z) a(s-1)`, so
`Int a(s)(-q)^s ds` solves the Picard-Fuchs equation `[prod_j (zD - zeta^j) - q prod_(k=1..5)(5zD + kz)] psi = 0` (PROVED, standard; the
residues at `s = zeta^i/z + d` give the Gamma-normalised Frobenius solutions `Z_i`, those at `s = -k/5` the solutions at `q = infinity`).
**Stirling form.** `Gamma_St(y) := exp((y - 1/2) log y - y + (1/2) log 2pi + Sum_k b_k y^(1-2k))`, `b_k = B_2k/(2k(2k-1))`, satisfies
`Gamma_St(y+1) = y Gamma_St(y)` as formal series in `1/y` (classical). Replacing all six Gammas, the `log z` terms cancel and
`a_St(s)(-q)^s dsigma = (elementary) exp(S_0/z + S_1 + Sum_k z^(2k-1) S_2k) dsigma` with
`S_0' = 5 log(5 sigma) - log(1 - sigma^5) + log(-q)`, `S_1 = (1/2) log(5 sigma/(1 - sigma^5))`, `S_2k = b_k[(5 sigma)^(1-2k) + Sum_j (zeta^j - sigma)^(1-2k)]`.
The saddles are `u sigma^5 = 1`; colour `i` is `sigma = zeta^i L`, and `S_0(zeta^i L) = u_i + const` (`D S_0 = sigma`).
**One-dimensional form (PROVED, direct computation).** In `x = 1/sigma` (logarithms combine; the branch is fixed by the colour-0 saddle):

    Omega = exp(W(x; u)/z) A(x; z) dx,     W = log(1-u)/x + G(x),    G(x) = Int_0^x log(1 - t^5) t^-2 dt = -Sum_l x^(5l+4)/((l+1)(5l+4)),
    A(x; z) = (1 - x^5)^(-1/2) exp( Sum_k b_k z^(2k-1) x^(2k-1) [ 5^(1-2k) + Sum_j (zeta^j x - 1)^(1-2k) ] ),     A(0; z) = 1.

`W' = x^-2 log((1 - x^5)/(1 - u))` vanishes exactly at `x^5 = u`; `W(zeta x) = zeta^-1 W(x)`; `zD Omega = x^-1 Omega` (`D log(1-u) = 1`).
The `j = 0` Stirling term `(z x/(1 - x))^(2k-1)` is singular at `x = 1` (`q = 0`, colour 0) and `(z x/5)^(2k-1)` comes from `Gamma(5s+1)`.
**Hatted form.** `x = s xi`: `W/z = What(xi; u)/zhat`, `What = (log(1-u)/u) xi^-1 - Sum_l u^l xi^(5l+4)/((l+1)(5l+4))`, and `A` is a power
series in `u` with polynomial coefficients (`z^(2k-1) x^(2k-1) = (u zhat xi)^(2k-1)`); the colour-0 saddle is `xi = 1` for all `u`.

## 3. Theorem A: the canonical normalisation is the full-Stirling normalisation

Rows: `Omega_0 = Omega/C_0`, `Omega_(a+1) = (x^-1 Omega_a + z (D rho_a) Omega)/c_a` (`rho_a` the coefficient of `Omega`), with the flat-frame
`c = (C_1, L^5/(C_0^2 C_1^2), C_1, C_0, C_0)` of UR s.1 (`C_0 = varpi_0`, `C_1 = D t_c`), any conifold frame.

**Theorem A (PROVED given Theorem U, UR s.3; proof `GAMMA_PROOFS.md` s.1-s.5).** Let `Y^MB_a := Res_(x=0) Omega_a` (order by order in `u`)
and `Q^fS_(a,i) :=` the formal stationary-phase expansion of `Omega_a` at the colour-`i` saddle (natural Gaussian normalisation). Then
* (i) `Y^MB = Y^C_Q`, the vanishing solution of GAUGE s.2 (value `-e_1` at `u = 0`; in particular the normalisation constant is `K = 1`);
* (ii) `Q^fS = P_can Q^m_nat`, `Q^m_nat` the model's stationary-phase WKB of `Int exp(-(X^4/4 + v/X)/z) X^m dX` at `v(u) = u tau^(-5/4)`.

Ingredients: (a) PF for both functionals (difference equation + formal translation invariance); (b) (H) for `Q^fS` (hatted structure),
so `Q^fS (Q^m)^-1 = P_can / c^fS` by Theorem U; (c) **Lemma CFU**: a unique hatted map `xi = Xi psi(Xi^5; u)`, `psi = 1 + Sum_k u^k phi_k`,
`phi_k in Q[y]`, `phi_k(1) = 0`, with `What(xi; u) = -(Xi^4/4 + Xi^-1)/tau(u)`; solvable at each order because the critical values agree:
`What(1; u) = -Sum_n 5u^n/(5n+4) = -(5/4)/tau(u)` (exact); (d) Brieskorn reduction of the pulled-back rows to the model rows `X^m dX`,
`m in {0,-1,3,2,1}`: ONE matrix `P^MB` intertwines stationary phase and residue; (e) `Y^MB_1(0) = -A(0; z) = -1`, so Theorem G gives `c^fS = 1`.
**Corollary A.** `D(z) = c(z; L)/R(z; L)` for every `L`, `c(z; L) := Sum_n c_n(L) z^n` the colour-0 full-Stirling saddle series (`c_0 = 1`).
(H) + (N) to all orders <=> `c = R` <=> `c(z; 1) = h(z)`. The vanishing period formula of s.0 item 6 is `[z^-1]` of (i).
**VERIFIED.** (i): `gm_residue.py 26 20 [1/3 -2/7]`: `Res Omega = C_0 Y^C_Q,0` for `z^-1..z^-20` mod `u^26` and row 1 for `z^0..z^-19`,
exact over Q, frames (0,0) and (1/3,-2/7); control without Stirling corrections fails. (c): `gm_cfu.py 12`: all orders `k <= 12`, `deg phi_k = k`
(so the flat map is even jointly analytic), control `tau = 1` obstructed at `k = 1`. (ii) through `z^24`: E4 + (Gamma).

## 4. Identity (Gamma)

### 4.1 Statement
**(Gamma)** For every `n >= 0`: `c_n(L) = R_n(L)` in `Q[L]`. Equivalently (both sides are colour-0 WKB series, so their ratio is a constant
series): `c(z; L0) = R(z; L0)` at ONE value `L0`, or `c(z; 1) = h(z)`: the full-Stirling saddle series, continued to `q = 0` (where the
saddle meets the poles of `Gamma((1 - sigma)/z)`), equals the Stirling factor of the OTHER five Gammas at `sigma = 1`. In words: the singular
Gamma contributes exactly nothing. The analogous statement at `q = infinity` (`Gamma(5s+1)` singular, first pole `s = -1/5`) is `c(z; 0) = Ginf(z)`.

### 4.2 Evidence (VERIFIED; `gm_compare.py`, logs `gm_compare_30_125.log`, `gm_compare_20_200.log`)
* exact over Q: `c_n(L) = R_n(L)` for `n <= 30` at `L = 2, 1/2, 7/3` (`L = 1/2` is `q < 0`); (H) h_auto and resonances for `n <= 30`;
* mod p61 and p89: `D(z) = c(z;L)/R(z;L) = 1 + O(z^201)` at `L = 2, 3/5`; h_auto for `n <= 200`; free constants `= [z^(5m)] Ginf`, `m <= 40`;
* controls: shifting `x_1` by 1 gives `D - 1` first non-zero at `z^5`; shifting `x_40` by 3 gives it at `z^200` (as predicted).
Cost: `n <= 125` in 61 s per prime; the comparison is `O(N^4)` and needs no gauge, no model and no periods.

### 4.3 Proof status: SKETCH (two routes, details in `GAMMA_PROOFS.md` s.6-s.7)
* **Formal route.** Near `sigma = 1` write `sigma = 1 - z eta`; the integrand is `Gamma_St(eta) Y^-eta exp(V)`, `V = Sum_m z^m P_m(eta)`.
  *Single-Gamma lemma* (PROVED): the formal saddle series of `Gamma_St(eta) Y^-eta eta^r` is EXACTLY `(-Y d/dY)^r e^(-Y)` (it obeys
  `-Y dF/dY = Y F` by `Gamma_St(eta+1) = eta Gamma_St(eta)` and translation invariance), i.e. the residue sum of the exact `Gamma`. Matching the
  outer series (fixed `L`) with the perturbative inner one inside the overlap ring `Q[[z/eps, eps^2/z]]` (`eps = L - 1`) gives `c(z;1) = h(z)`;
  the reorganisation lemma for formal Laplace expansions in that ring is not written out (SKETCH).
* **Analytic route.** (A1) The Gamma-normalised Frobenius solution is a lattice sum `Z_0 = Sum_(d>=0) E(1 + dz)` of a positive, log-concave
  analytic function (reflection formula: `1/Gamma` is entire); Euler-Maclaurin + Laplace give `Z_0 ~ (elementary) e^(u_0/z) L c(z; L)` at
  fixed `q in (0, q_c)` (PROVED modulo routine estimates). (A2) `Z_0 ~ (elementary) e^(u_0/z) L R(z; L)` at fixed `q`: uniform validity of
  the Frobenius-normalised WKB down to the regular singular point `q = 0` (classical uniform simplification for leading matrices with
  distinct eigenvalues `zeta^j`; precise citation not checked). (A1) + (A2) give (Gamma). No purely ODE-formal argument can decide (Gamma):
  the ratio of two WKB series is an arbitrary constant series; an exact (analytic or overlap) input is needed.

## 5. What this means for Stage 4 and Conjecture N

* (F_g) (STAGE2_SYNTHESIS s.2) is `D == 1 mod z^(3g-2)`. By Theorem A and 4.2 it holds for `g <= 67` (VERIFIED, modulo Theorem U and the
  AI-written proof of Theorem A). Conjecture N(g), `g <= 67`, then rests on (U), the quoted theorems of STAGE2_SYNTHESIS s.2 and Theorem A.
* For all `g`: Conjecture N <= (U) + quoted theorems + Theorem A + (Gamma). (Gamma) is a statement about Stirling series of six Gamma
  functions; it no longer involves the gauge, the model, the periods or graph sums.
* The total-potential obstruction (MECHANISM s.2) is respected: Theorem A is graph-free and fixes the gauge itself.
* GAUGE s.6 called (N) "a connection problem between q = 0 and the conifold". It splits: conifold <-> full-Stirling MB (Theorem A, formal,
  proved via the vanishing cycle = residue at the pole of `Omega`), and full-Stirling MB <-> q = 0 (identity (Gamma), a one-Gamma connection).

## 6. Literature (read for this note)

* arXiv:1803.01409 (Lho-Pandharipande, formal quintic) s.2.4.4 eq. (18) and Prop. 5: the I-function at `H = lambda_i` is
  `e^(mu lambda_i/z) Sum_k R_k (z/lambda_i)^k`, `R_0 = L`, `R_1 = (3/20)(L - L^5)`, `R_k in C[L^(+-1)]` (via arXiv:0710.0889 Thm 4); formal.
  (arXiv:1702.06096, cited in the brief, is the local P^2 paper; the formal quintic paper is arXiv:1803.01409.)
* arXiv:0710.0889 (Zagier-Zinger) Thms 1-4: for `F(w,x) = Sum_d x^d prod_(r=1..5d)(5w+r)/prod_(r=1..d)((w+r)^5 - w^5)` (our `f_0`, `w = 1/z`)
  the expansion `F ~ e^(mu w) Sum_s Phi_s w^-s` with `Phi_0 = L`, `Phi_s in L Q[L]`, `Phi_s(0) = delta_(s,0)`, as formal series (coefficientwise in `x`).
* arXiv:math/0110142 (Coates-Givental, QRR): source of the q = 0 factors `h_w`; arXiv:0903.1463 (Iritani, Gamma integral structure) and
  arXiv:1404.6407 (Galkin-Golyshev-Iritani): Gamma-normalised solutions; we use only the elementary Mellin-Barnes form of s.2.
* Uniform simplification at a regular singular point with a large parameter: the Iwano-Sibuya-Nishimoto line and Wasow's monograph
  (no arXiv numbers; statement not checked here). Second-order analogue: Debye expansions of `J_nu(nu x)`, uniform as `x -> 0`.

## 7. Files and runs (all `py -I -B`, stdlib, exact; working directory outside the project)

| file | content | run / log | result |
|---|---|---|---|
| `gm_series.py` | truncated series over Q / F_p, Bernoulli, Lagrange inversion | library | - |
| `gm_qrr.py` | independent colour-0 WKB polynomials (HPF + QRR), `h`, `R_n(L)` | library | - |
| `gm_fs.py` | full-Stirling MB saddle series `c(z; L)` at rational `L` | library | - |
| `gm_compare.py` | (Gamma), (H), free constants `= Ginf`, controls | `30 125`, `20 200` | ALL PASS |
| `gm_residue.py` | Theorem A(i) on the actual vanishing solution, period formula, control | `26 20`, `26 20 1/3 -2/7` | ALL PASS |
| `gm_cfu.py` | Lemma CFU to `u^12`, control | `12` | ALL PASS |

Read-only imports (only `gm_residue.py`): `natural_form/theory/ur_gauge.coefficients`, `allgenus/gauge/gi_core.quintic_vanishing`.
Timings: `gm_compare 30 125` 141 s, `gm_compare 20 200` 823 s (that log was produced by a byte-identical copy in the session scratch
folder and copied here), `gm_residue` and `gm_cfu` a few seconds.
Nothing else in the project was read by code or modified.

## 8. Caveats

* Theorem A is proved here (AI-written, unrefereed); it uses Theorem U (AI-audited only) and GAUGE's Theorem G and Lemma 1.
* (Gamma) is VERIFIED to `z^200`, not proved; both proof routes are sketches. The q = infinity identity `x = Ginf` is VERIFIED only.
* The formal stationary-phase and residue functionals on `Omega` are defined order by order (in `z`, resp. `u`); `GAMMA_PROOFS.md` s.1
  lists the properties used (killing exact forms, invariance under analytic changes of variable, translation invariance).
* Frames: Theorem A holds in every conifold frame `(alpha, beta)` (`C_1(0) = -1` in all of them); the residue check ran in two frames.
