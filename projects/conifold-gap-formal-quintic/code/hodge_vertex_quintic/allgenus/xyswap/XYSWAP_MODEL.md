# The LG model x^4 + v/x through the x-y swap: closed formulas, Bernoulli free energies, regular legs (2026-11-13)

Folder `allgenus/xyswap/`. AI-derived and AI-checked only; no human has read it. Labels: PROVED (argument written here, or a quoted
theorem with arXiv number, statement and checked hypotheses), VERIFIED (exact computation; orders, field/primes stated), CONJECTURAL,
REFUTED. Python 3.14, stdlib only, `py -I -B`, at most two processes, every run < 1 minute.

**Curve and conventions.** `x = z^4 + v/z`, `y = z` on `P^1`, `B = dz dz'/(z - z')^2` (`z` = the task's `p`; on the line of the project
`v = a_5`, the residue coordinate). Branch points `4 z^5 = v`. Code convention (`xy_tr.py`): `omega_{0,1} = y dx`,
`K(z0, z) = (1/2) Int_{sigma z}^{z} B(z0, .) / ((y(z) - y(sigma z)) dx(z))`. The recursion (1.1) of arXiv:2212.00320 (ABDKS) is written with
`Int_z^{sigma z}` and `omega_{1,0} = -y dx`; VERIFIED: `omega^code_{g,n} = (-1)^n omega^ABDKS_{g,n}` (all cases below). `F_g` (g >= 2) is
*defined* by the dilaton formula `F_g = (1/(2-2g)) Sum_a Res_a Phi omega_{g,1}`, `dPhi = y dx`, `Phi = (4/5) z^5 - v log z`. With the
opposite kernel orientation and the same `Phi`, every `F_g` (g >= 2) changes sign; the code convention is the one in which the
project's normalisation of the frozen theory (LG_MODEL.md s.5.4: `c_2 = -1/240`, `c_3 = 1/1008`) comes out.
Leg pairing: a leg is a function `Phi` on the curve (its differential `dPhi` is the leg differential);
`<omega_{g,n}, Phi_1 .. Phi_n> := Sum Res_{a_1} .. Res_{a_n} prod_j Phi_j(z_j) omega_{g,n}`. It vanishes on `Phi = f(x)` for any `f` holomorphic
near the branch points (`Res_a f(x) omega = 0` by the linear loop equation), so legs matter modulo functions of `x`.
Dictionary used for legs (fixed by the TR string and dilaton equations, VERIFIED below): `tau_0(1) <-> -z`, `tau_0(H) <-> -log z`
(`H = d/dv`), `tau_k(a) <-> Int tau_{k-1}(a) dx`; so `tau_1(1) <-> -(4/5) z^5 + v log z`, `tau_1(H) <-> -(x log z - z^4/4 + v/z)`.

## 0. Results in brief

1. **Which swap theorem applies (s.1).** arXiv:2212.00320, **Theorem 7.1**: for `S = P^1`, `y = z` and ANY rational `x` (with simple zeros
   of `dx`) the TR differentials are given by an explicit closed formula. Our curve is literally their example (7.2), "enumeration of
   hypermaps", `x = z^(r-1) + z^(-1)`, with `r = 5` (rescale `z`). `x`, `y` are meromorphic, `dx`, `dy` have no residues: neither log TR
   (arXiv:2312.16950) nor the generalized TR of arXiv:2408.02608 is needed. The dual side (`x^vee = z` unramified) is trivial; the
   simple pole of `x` (residue `v` of `x dz`) enters only through the vertex factor
   `G = exp(-w (S(w hbar d_z) - 1)(v/z)) = ((1 - w hbar/2z)/(1 + w hbar/2z))^(v/hbar) e^(v w/z)`.
   The swap theorems cover `n >= 1` only; the free energies are not swap invariant here (`F^vee_g = 0`, `F_g != 0`).
2. **Closed formulas = own direct TR (VERIFIED).** Exactly over `Q(zeta_5)`: `(0,3), (1,1), (0,4), (1,2), (2,1)`; mod `2^89-1`: also
   `(3,1), (2,2), (1,3)`. Explicit formulas for `(0,3), (1,1), (1,2), (2,1)` in s.2 (checked exactly against TR).
3. **Free energies (PROVED mod Thm 7.1; VERIFIED).** On the whole family `x = z^4 + c3 z^3 + c2 z^2 + c1 z + c0 + V/z` (the Hurwitz space
   `H_{0;1,4}` with the pole moved to 0): **`F_g = B_2g/(2g(2g-2)) V^(2-2g)` exactly, for all `g >= 2`** - no dependence on `c3, c2, c1, c0`.
   So `F_g` of the Hurwitz Frobenius manifold depends on the residue only (the "primary-sector lemma" of UNIFORM_REDUCTION s.6.5 item 3,
   for `g >= 2`, on the whole manifold). Mechanism: the polynomial part of `Phi` contributes 0 (degree), the `-V log z` part gives `V`
   times the vanishing period `Int_inf^0 omega_{g,1} = (B_2g/2g) V^(1-2g)`, which is a Borel-summed Gaussian vertex,
   `Sum_g [eps^2g] Int_0^inf e^(-u) eps^2/sinh^2(eps u) du`. Also `F_1 = -(1/12) log v` on the line (EO variational formula).
   VERIFIED: `F_2`, `F_3` by direct TR on the line (`v = 4`; `Q(zeta_5)` and `2^89-1`) AND at generic points of the 4-parameter family
   (own generic TR code, branch points chosen in `F_p`, no symmetry, both primes), where the polynomial part of `Phi` contributes
   exactly 0 (control); the period for `g <= 10` from the closed formula.
4. **Regular legs (s.3).** Class `L = Q[params][z] + Q[params][x]` (single valued, no `log z`, regular at the pole modulo functions of
   `x`, coefficients polynomial in the parameters). **Model Lemma part (ML1), PROVED mod Thm 7.1:** for `2g-2+n > 0`, every correlator
   with ALL `n >= 1` legs in `L` is a polynomial on the Hurwitz space (in particular regular at `V = 0`). VERIFIED at `(1,1), (2,1), (1,2),
   (2,2)` exactly and mod two primes: 28 non-zero `L`-correlators all with exponent `>= 0`, and 11 `L`-correlators that homogeneity would
   allow to be singular vanish; every bad leg (`log z`, the dilaton, `tau_1(H)`, `z^-5`) gives singular values.
5. **String equation reconciled.** `tau_1(H)` is a logarithmic leg (monodromy `x` around the vanishing cycle `|z| = eps`); the unit leg
   (in `L`) acts on the other legs by `-d/dx`, which maps `tau_1(H)` to `log z`. VERIFIED: `<omega_{g,2}, tau_0(1), tau_1(H)> = <omega_{g,1},
   tau_0(H)>` = `-(1/12) v^-1` (g = 1), `(1/120) v^-3` (g = 2). Mixed correlators (one `L` leg, one bad leg) are singular in general.
6. **Euler-dilaton = log-free leg.** `e = tau_1(1) + v tau_0(H) <-> -(4/5) z^5` (the polynomial part of `-Phi`), in `L`; TR identity
   `Sum Res (4/5) z^5 omega_{g,n+1}(z, J) = (2-2g-n - v d_v|_x) omega_{g,n}(J)`; VERIFIED `<e, H>_1 = <e, H>_2 = 0`. The edge of the Euler-dilaton dressing that the project
   computed (UNIFORM_REDUCTION s.6.3(ii)) has all its legs in `L` at the model point. Whether ALL edges of the regular gauge are in `L`
   was not tested (s.4): it is the remaining step of the mechanism, together with a statement on the leaves.

## 1. The x-y swap theorems and which one applies

* **arXiv:2212.00320** (Alexandrov, Bychkov, Dunin-Barkowski, Kazarian, Shadrin), read: s.1 and s.6-7. *Theorem 1.8/1.14*: if the
  `omega^{(g)}_{m,0}` satisfy TR for `(S, x, y, B)` with meromorphic `x`, `y` (standing assumptions s.1.1: critical points of `x` simple,
  `y` regular with `dy != 0` there), all zeros of `dy` simple and `x` regular at them, then the `omega^{(g)}_{0,n}` given by the universal
  graph formula (1.13) (`n >= 1`) satisfy TR for `(S, y, x, B)`. *Theorem 7.1*: for `S = CP^1`, `y = z`, `x` an arbitrary rational
  function (s.7: "for all possible rational x such that dx has only simple zeros"), with initial form `-y dx`,
  `omega^{(g)}_{m,0} = [hbar^2g] Sum_{r_i>=0} prod_i D_{x_i}^{r_i} [w_i^{r_i}] Sum_Gamma hbar^{2g(Gamma)} prod_{(ij)} W_ij prod_i W_i`,
  `W_i = -e^{-w_i (S(w_i hbar d_{z_i}) - 1) x(z_i)} dz_i/w_i`, `W_ij = w_i w_j/((z_i - z_j)^2 - (hbar^2/4)(w_i + w_j)^2)`, `S(u) = sinh(u/2)/(u/2)`,
  `D_x omega = d(omega/dx)`, `Gamma` connected simple graphs on `m` numbered vertices. *Hypotheses for our curve:* `dx = (4z^5 - v) dz/z^2`
  has five simple zeros (`v != 0`), `y = z` is regular with `dy = dz != 0` everywhere, `x` rational with poles at `0` (simple) and `inf`
  (order 4): satisfied. Their Remark 7.3: the formulas extend naturally to the discriminant (here `v = 0`).
* **arXiv:2211.08917** (Hock), abstract: a simpler formula relating `W_{g,n}` and `W^vee_{g,n}` (genus 0 in that paper; the general-genus
  version is the formula (1.13) above, also proved by Hock, arXiv:2304.03032 per ABDKS s.1.5), for meromorphic `x`, `y` with simple and
  distinct ramification points. Consistent with the above; not needed separately.
* **arXiv:2312.16950** (ABDKS, log TR), abstract: a patch of TR in the presence of logarithmic singularities (`dx` or `dy` with simple
  poles with residues), for which the x-y swap relation is proved. Not needed: here `x`, `y` are meromorphic. The pole of `x` produces a
  logarithm only in `Int x dz = z^5/5 + v log z` on the *dual* side, and Theorem 7.1 handles it through `S(w hbar d)` acting on `v/z`.
* **arXiv:2408.02608** (ABDKS, generalized TR for degenerate/irregular curves), abstract only: agrees with ordinary TR in the regular
  case, which is ours.
* **Free energies.** None of the above covers `n = 0`. Naive symplectic invariance fails: `F^vee_g = 0`, `F_g = c_g v^(2-2g)`. The whole
  discrepancy is the Gaussian term of the pole (s.2.3). (Literature on x-y swap of `F_g` not searched beyond the four abstracts.)

**Decision.** Theorem 7.1 of arXiv:2212.00320 applies verbatim; all statements below are conditional on it (and on the identification
"model CohFT = TR on this curve", assumed by the task and checked by the `tr` worker; here only consistency: `F_1`, `F_2`, `F_3`,
`<tau_1(1)>_1 = 5/24`, `<e>_1 = 1/8` agree with the project's numbers up to the sign conventions of UNIFORM_REDUCTION s.6.3).

## 2. Closed formulas

### 2.1 General structure (from Theorem 7.1)
`omega_{g,n} = (-1)^n [hbar^2g] Sum_r prod_i D_{x_i}^{r_i}[w_i^{r_i}] Sum_Gamma hbar^{2g(Gamma)} prod W_ij prod_i (-dz_i/w_i) A(w_i, z_i) G(w_i, z_i)`,
* `A(w, z) = exp(-(1/2) w^3 hbar^2 z^2 - (1/80) w^5 hbar^4)`: the vertex of `x = z^4` (the `v -> 0` limit, the 4-spin type curve);
* `G(w, z) = exp(-(v/hbar)[log((1 + a)/(1 - a)) - 2a])`, `a = w hbar/(2z)`, `= exp(-v Sum_{k>=1} w^(2k+1) hbar^(2k)/(4^k (2k+1) z^(2k+1)))`:
  the Gaussian factor of the simple pole; `v` enters only here and in `D_x = d(z^2 . /((4 z^5 - v) dz))`.
So `omega_{g,n} = N_{g,n}(z; v) prod dz_i / prod (4 z_i^5 - v)^(6g-4+2n)` with `N` a polynomial in `z` and `v` (PROVED: each term of (7.9) is a
Laurent polynomial in the `z_i` over powers of `4z_i^5 - v` and `(z_i - z_j)`; the total has no diagonal poles and is regular at `z_i = 0`
for `v != 0`, where `4z^5 - v` is a unit; the same for the whole family with `q(z) = z^2 x'(z)` in place of `4 z^5 - v`).

### 2.2 The four requested cases (code convention; `h = z^2/2 + v/(12 z^3)`, `k = 1/80 + v/(80 z^5)`) - VERIFIED exactly against TR
* `omega_{0,3} = Sum_i D_{x_i}[ dz_1 dz_2 dz_3 / prod_{j != i} (z_i - z_j)^2 ]`  (2 points, `Q(zeta_5)`)
* `omega_{1,1} = -D_x^2(h dz) = -5 z^4 (3 v^2 + 36 v z^5 + 8 z^10) / (4 z^5 - v)^4 dz`  (3 points)
* `omega_{1,2} = -D_1^3[h(z_1) dz_1dz_2/(z_1-z_2)^2] - D_2^3[h(z_2) ...] + (1/4)(D_1^2 + 2 D_1 D_2 + D_2^2)[dz_1dz_2/(z_1-z_2)^4]`  (2 points)
* `omega_{2,1} = -D_x^4(k dz) + (1/2) D_x^5(h^2 dz) = [266112 z^34 - 10384448 v z^29 - 20583920 v^2 z^24 - 7658960 v^3 z^19 - 690770 v^4 z^14
  - 12053 v^5 z^9 - 8 v^6 z^4] / (4 z^5 - v)^10 dz`  (3 points)
In `h` and `k` the first term is the 4-spin vertex (`A`), the `v`-term the Gaussian vertex (`G`). Files: `xy_explicit.py`, `.log`.

### 2.3 The Gaussian part: the vanishing period and `F_g` (PROVED mod Thm 7.1)
For `n = 1` only the single vertex occurs and `omega^ABDKS_{g,1} = d Phi_g`, `Phi_g = Sum_{r>=1} T_r`, `T_r = D_x^(r-1)([w^r] W_1)/dx` (`[w^0] W_1 = 0`).
*Monomial count* (any `P = z^4 + c3 z^3 + c2 z^2 + c1 z + c0`, pole `V/z`): a factor of `[hbar^2g] e^{-w(S-1)x}` is either Gaussian
(`V w^(2k+1) hbar^(2k) z^-(2k+1)`) or non-Gaussian (`w^3 hbar^2 z^e`, `e = 2, 1, 0`, or `w^5 hbar^4`), and `1/x' = -(z^2/V)(1 - z^2 P'(z)/V)^-1`.
Near `z = 0` a product of Gaussian factors gives `z^(m + r + 1)` with `m + r + 1 = 0` exactly, and every non-Gaussian factor or correction
of `1/x'` raises the power by at least 2. Near `z = inf`, `1/x' ~ z^-3/4` and `T_r = O(z^(5-10g))`. Hence each `T_r` is regular at 0,
`Phi_g(inf) = 0`, and `Phi_g(0)` is the constant term of the pure pole `x = V/z`: in `xi = V/z`, `T_r = V (2V)^-2g (r-1)! [eps^2g][t^(r+1)] e^(E_0)`,
`E_0 = -Sum eps^2k t^(2k+1)/(2k+1)`, `eps = hbar/(2V)`; Borel summation (`(r-1)! = Int_0^inf e^-s s^(r-1) ds`) and `eps s = tanh(eps u)` give
`Sum_r (r-1)! [t^(r+1)] e^(E_0) = Int_0^inf e^-u eps^2/sinh^2(eps u) du`, and `(x/sinh x)^2 = Sum (1-2g) 2^2g B_2g x^2g/(2g)!`:
**`Phi^ABDKS_g(0) - Phi^ABDKS_g(inf) = -(B_2g/(2g)) V^(1-2g)`**, i.e. `Int_inf^0 omega^code_{g,1} = (B_2g/2g) V^(1-2g)` (g >= 1).
*Free energy.* `Phi = Phi_pol - V log z`, `Phi_pol = z P(z) - Int P + V` (polynomial). `Sum_a Res_a Phi_pol omega_{g,1} = -Res_inf(Phi_pol omega_{g,1})`
is a polynomial in `(c, V)` (s.2.1) of weighted degree `5(2-2g) < 0` (weights `c_k: 4-k`, `V: 5`), hence 0. For the log: the contour
around a cut from 0 to inf gives `Sum_a Res_a(log z omega) = Phi(0) - Phi(inf)` for the primitive `Phi` of `omega` (residue-free,
holomorphic at 0 and inf). Therefore **`F_g = (1/(2-2g)) (-V) Int_inf^0 omega^code_{g,1} = B_2g/(2g(2g-2)) V^(2-2g)`** for `g >= 2`, on the
whole family. `F_1`: the EO variational formula with `d_v(y dx)|_x = -dz/z` gives `d_v F_1 = Int_0^inf omega_{1,1} = -1/(12 v)` (the sign of
that formula fixed by consistency with the dilaton definition at `g >= 2`), i.e. `F_1 = -(1/12) log v` on the line.
*Verification.* `xy_run1_*.log`: closed-formula `Phi_g(0) = -B_2g/(2g) v^(1-2g)` and `T_r` regular at 0 and vanishing at inf, `g = 1..10`
(exact, symbolic `v`); `F_2 = -1/3840`, `F_3 = 1/258048` at `v = 4` by own direct TR (`Q(zeta_5)` for `g = 2`, `2^89-1` for `g = 3`).
`xy_trgen_*.log`: own TR at generic points of the family (branch points chosen in `F_p`, `c3 c2 c1 V != 0`, no symmetry used):
`F_2 = c_2 V^-2` (2^89-1 and 2^61-1, two seeds), `F_3 = c_3 V^-4` (both primes); control: the polynomial part of `Phi` contributes
exactly 0 and dropping `-V log z` destroys the identity. `xy_legs_*.log`: `Int_inf^0 omega_{1,1} = (1/12) v^-1`, `Int_inf^0 omega_{2,1} =
-(1/120) v^-3` from TR. **The Gaussian structure** is thus: `G` is the c = 1 vertex; the vanishing period is its Borel sum
`Int e^-u (eps/sinh(eps u))^2 du` (the Schwinger/Gopakumar-Vafa integral of the c = 1 string), with the argument `w hbar/z` entering
non-linearly (through `artanh`) and `v` as `v/hbar` in an exponent (the size `N` of a Gaussian matrix integral).

## 3. Singular structure at v = 0 and the model lemma

* **Fixed `z`.** `omega_{g,n}(z; v)` is a polynomial in `v` over `prod (4z_i^5 - v)^K`: regular at `v = 0` for `z_i != 0`, with limit the
  closed formula of `x = z^4` (vertex `A` only). All singular behaviour comes from `z -> 0`, where `1/(4z^5 - v)` expands in powers of
  `z^5/v` (the inner region `|z| < |v|^(1/5)`, in which `x ~ v/z` and the curve is the Gaussian one).
* **(ML1) PROVED mod Thm 7.1.** For legs `Phi_j in L` the pairing is `(-1)^n Res_inf .. Res_inf prod Phi_j omega_{g,n}` (no poles at 0;
  functions of `x` drop out), a polynomial in the parameters since `1/q(z)` expands at `inf` with polynomial coefficients. This holds on
  the whole family (`L = Q[c, V][z] + Q[c, V][x]`). Note `v/z = x - z^4 in L` (the `H^2` primary) but `1/z = (x - z^4)/v` is not.
* **(ML2) PROVED mod Thm 7.1.** `F_g = c_g V^(2-2g)` exactly on the whole family (`g >= 2`); `Int_inf^0 omega_{g,1}` pure Gaussian.
* **(ML3) the singular part of a general leg.** A leg with poles at `z = 0` decomposes as `Phi = Phi_L + log z . f(x) + R`: `log z . f(x)`
  (monodromy `2 pi i f(x)` around `|z| = eps`; the `H`-descendants) pairs through weighted vanishing periods `Int_inf^0 f(x) omega`;
  `R` (zero-residue poles not removable by `Q[c,V][x]`) is rewritten with `1/z = (x - P(z))/V`, which introduces explicit `V^-1`.
  So the singular part of every correlator is a polynomial combination of Gaussian vanishing periods and of `V^-1` from this rewriting.
  (Statement: SKETCH; the weighted periods `Int f(x) omega_{g,n}` were not evaluated in closed form beyond `f = 1`.)
* **VERIFIED (`xy_legs.py`, exact `Q(zeta_5)` and both primes, identical logs).** Good legs `-z, -z^2/2, -z^3/3, v/z, -(4/5)z^5, z^9, z^10,
  z^14, z^15, z^19, z^20, x^2 + 3z^8`: all non-zero correlators at `(1,1), (2,1), (1,2), (2,2)` have exponent `>= 0` (e.g. `<omega_{1,1},
  -(4/5)z^5> = 1/8`, `<omega_{2,2}, z^10, z^10> = -165/64`), and `<omega_{2,1}, z^5>`, `<omega_{2,1}, z^10>`, 2 pairs at `(1,2)`, 7 pairs at
  `(2,2)` vanish although homogeneity would allow `v^-1 .. v^-3`. Bad legs: `<omega_{g,1}, -log z> = -(1/12)v^-1, (1/120)v^-3`;
  `<omega_{2,1}, tau_1(1)> = -(1/120) v^-2 = 2 F_2` (dilaton); `<omega_{1,1}, z^-5> = 15 v^-2`; `<omega_{1,2}, -log z, -log z> = (1/12) v^-2 =
  F_1''`, `<omega_{2,2}, -log z, -log z> = -(1/40) v^-4 = F_2''`; mixed `<omega_{1,2}, -(4/5)z^5, z^-5> = 60 v^-2`: one good leg does not suffice.
* **Reconciliation with `<tau_1(H) tau_0(1)>_h = <H>_h`.** In TR, `tau_1(H) = -(x log z - z^4/4 + v/z)` is logarithmic, so it is outside
  `L`. The TR string equation makes the unit leg act by `-d/dx` on the other legs and `-d/dx` of `tau_1(H)` is `-log z = tau_0(H)`: the
  singularity is carried by the monodromy of the leg, which the vanishing covector "at the base point" does not see (the LG period
  `Y^C = oint_{|z|=eps} e^{lambda/z}(.)` and the TR monodromy are different functionals: `Y^C` pairs non-trivially with `L` legs such as
  the unit, the monodromy does not). The lemma needs ALL legs in `L`; mixed correlators can be singular (verified above).

## 4. Which legs see only the regular part; the regular gauge

* **Answer for the class.** `L` (single valued, no `log z`; polynomial in `z` modulo polynomials in `x`, coefficients regular). Zero
  residue at the pole is NOT enough (`z^-5`, `1/z` give `v^-1`), and single valuedness is necessary (`tau_1(H)`).
* **The Euler-dilaton block.** `tau_1(1) + v tau_0(H) <-> -(4/5) z^5` (the logs cancel): the Euler-dilaton vector of UNIFORM_REDUCTION
  s.6.4 is the log-free part of `-Phi`, and TR gives `Sum_a Res_a (4/5) z^5 omega_{g,n+1}(z, J) = (2-2g-n - v d_v|_x) omega_{g,n}(J)` (`d_v` at fixed `x`, EO variational formula), which
  annihilates the `H`-sector (`<e, H^n>_g = 0`) and gives the regular genus-1 value `<e>_1 = 1/8` (the project's anomaly `-1/8`, other
  sign convention). The edge that the project computed for the Euler-dilaton dressing (s.6.3(ii) there:
  `H(x)H + (psi - psi')(H(x)1 - 1(x)H) - (psi^2 - psi psi' + psi'^2) 1(x)1`) has half-edge legs `1`, `psi 1 + t H`, `psi(psi 1 + t H)` (with
  the model point `t`), i.e. `-z`, `-(4/5) z^5`, `-(2/5) z^9 + (v/5) z^4` at `t = v`: all in `L`. The chain `psi^k e(psi)` stays in `L` for
  `k <= 4` and acquires `log z` at `k = 5` (a `z^1` term at the fourth `x`-antiderivative): an infinite power series of such edges is
  in `L` only if the `O(u^2)` corrections of the block conspire, or if the logs come with enough powers of `u`.
* **Not done (the remaining step).** Translate every edge and leaf of the regular gauge `P` (ur_gauge.py) into legs with the dictionary
  above and check: (M1) edge legs in `L`, up to terms whose `u`-valuation exceeds the singularity they produce; (M2) the leaf leg is
  `delta(u) log z + L` (pure residue shift, no `x`-dependent monodromy). *Reduction (SKETCH):* (M1)+(M2), the standard identification of a
  Givental `R`-action with `B -> B + Delta B` and `y dx -> y dx + Lambda` (arXiv:1211.4021, arXiv:1110.2949; quoted from memory, not
  checked), and (ML1)+(ML2) give Conjecture N: the leaves move `V` to `V + delta` (residue of the deformed curve; on record it becomes
  `t_c`), the vertex graph is `F_g` of the deformed curve `= c_g (V + delta)^(2-2g) +` regular (ML2 plus ML1 for the `L`-part of `Lambda`), and
  every graph with an edge is a sum of all-`L` correlators of the deformed curve, regular by ML1. CONJECTURAL until (M1), (M2) are checked.

## 5. Matrix model / Toda (brief)

* **Toda (PROVED given s.2.3).** On the residue line `F(v) = Sum_g hbar^(2g-2) F_g` with `F_g = B_2g/(2g(2g-2)) v^(2-2g)`, `F_1 = -(1/12) log v`
  and `F_0 = (1/2) v^2 log v + ...` (`c_555 = 1/a_5`) is the asymptotic expansion of `log G(1 + v/hbar)` (Barnes; classical), hence satisfies
  the second-difference (Toda) equation `F(v + hbar) + F(v - hbar) - 2 F(v) = log(v/hbar)` - the GUE partition function with `N = v/hbar`.
  By ML2 this holds with `v = a_5` at every point of the Hurwitz space for the `g >= 2` part.
* **Matrix model (SKETCH, not checked).** The dual side is trivial, so the x-y swap is a Laplace transform of `exp((1/hbar) Int x dy)`,
  `x dy = (y^4 + v/y) dy`: wave function `psi(x) = Int dy y^(-v/hbar) exp((x y - y^5/5)/hbar)`, i.e. a generalized Kontsevich-Penner model
  `Int dM det(M)^(-v/hbar) exp(tr(Lambda M - M^5/5)/hbar)`. Its `det^(-N)` factor is the Gaussian/Penner part, `M^5/5` the 4-spin (A_3) part.
  The KP side (arXiv:2211.12259, rationally 4-constrained KP, hypermaps) is quoted in LG_MODEL.md s.5.5; nothing new here.

## 6. Files (all new, in this folder)

| file | content |
|---|---|
| `xy_field.py` | coefficient fields `Q(zeta_5)` (exact) and `Z/p` (from the first attempt; reviewed, unchanged) |
| `xy_local.py`, `xy_tr.py` | own direct TR on the line (`v = 4`, residue at `a_0 = 1`, Z/5 symmetry) |
| `xy_trgen.py` | own direct TR mod `p` at a generic point of the family (no symmetry); `F_2`, `F_3` |
| `xy_closed.py`, `xy_jets.py` | Theorem 7.1 of arXiv:2212.00320: `n = 1` symbolic in `z, v`; `n >= 2` exactly at rational points (jets) |
| `xy_run1.py`, `xy_run2.py`, `xy_explicit.py`, `xy_legs.py` | drivers; logs `xy_run1_qz.log`, `xy_run1_fp89_g3.log`, `xy_run2_qz.log`, `xy_run2_fp89.log`, `xy_explicit.log`, `xy_legs_{qz,fp89,fp61}.log`, `xy_trgen_{p89,p61}_g{2,3}.log` |

## 7. Caveats
* Everything rests on Theorem 7.1 of arXiv:2212.00320 (read, hypotheses checked) for all-genus statements; direct TR checks reach
  `g <= 3`, `n <= 3`. The identification of the model CohFT with this TR, and of Givental R-actions with deformations of `B` and `y dx`,
  are quoted, not re-derived. The F_1 sign uses the EO variational formula with its sign fixed by consistency, not re-proved.
* The class `L` and ML1 are about TR legs; their relation to the project's vanishing covector `Y^C` is only partly clarified (s.3).
* Sign conventions of `F_g`, of the kernel and of `tau_k` were fixed by regression against the project's `5/24`, `-1/8`, `c_2 = -1/240`.
