# Hatted translation of (G): the vanishing covector, the contraction lemma (L), the edge statement (E), the leaf statement (T)

Date 2026-11-13; companion of `GAUGE_IDENTITY.md` (labels and conventions there). AI-derived and AI-checked only. Graph-sum conventions:
`fq_sser.py` (rows, edges, leaves), `OU_EXTENSION_genus3.md` s.1 (hatted theory), `NATURAL_FORM.md` s.5 (outer graph sum `R'.(M.omega^)`).

## 0. Statements

* **Covector (PROVED closed form).** `phi_{i,m} = (u tau^(-1/2)/5) tau^(-m) zeta^(i(m-2)) / (4^q q! (m-q)!)`, `q = floor(m/5)`, `m >= 0`, colour `i`.
* **(L) Contraction lemma (PROVED; uses only the symplecticity of the rows of R').** Contracting any covector into one end of an edge leaves
  `(1/psi') R'(psi')^-1 [R'(x) f(x)]_(>0)|_(x = psi')` on the other end, `f_i(x) = Sum_m (-1)^m phi_{i,m} x^-m`.
* **(E) Edge statement (PROVED given Theorem U, for the canonical gauge; all genera, all u-orders).** `Sum_{i,m} phi_{i,m} E'(m, n)_{ij} = 0` for all `n, j`.
* **(T) Leaf statement (PROVED given Theorem U).** `Pi(T') := Sum_{i, d >= 2} phi_{i,d} T'^(i)_d = vhat - t_c`, `vhat = u tau^(-3/2)/xhat`.
* All four VERIFIED: the closed form over Q (`m <= 80`), (L), (E), (T) on the actual rows (s.6). What they do not give: regularity of the graph sum by
  themselves (s.5).

## 1. Hatted objects

`s = u^(1/5)`, `z = s^4 zhat`, `sigma = v(u)^(1/5) = s tau^(-1/4)` (so `z/sigma^4 = tau zhat`), `mu(b, k) = (4k + c_b - c_{b+k})/5`. The hatted gauge is
`Ghat(u, zhat) = Lambda P(u, s^4 zhat) Lambda'^-1`, `Lambda = diag(s^(4-c_a))`, `Lambda' = diag(s^3 sigma^(1-c_b))`, i.e. `[zhat^k] Ghat_{b+k,b} =
tau^((1-c_b)/4) u^mu(b,k) p_{k,b}(u)` (UNIFORM_REDUCTION s.4; checked by substitution). Both `Lambda` and `Lambda'` multiply `eta` by a scalar
(`c_a + c_{3-a} = 2`), so `Ghat` is symplectic iff `P` is. Since `mu(b, k) >= (4k - 4)/5`, `[u^N] Ghat` is a POLYNOMIAL in `zhat` (degree `<= (5N+4)/4`;
`<= floor((5N+2)/4)` by Theorem D). Rows of R' (flat lower index `a` x colour `i`): `Q'(u, x) = Psi(u) G(u, x) = Ghat(u, x) Psi0`, `Psi0[b][i] = eps_b
zeta^(i(b+1))`, `eps = (1,-1,-1,1,1)`, `Psi0^T eta^-1 Psi0 = 1`. In the production graph sum `R' = R~ M_0(tau zhat)^-1`; it equals the canonical
`Ghat Psi0` at every order used through genus 9 (GAUGE_IDENTITY s.6). Edge (fq_sser): `E'(psi, psi') = (1 - Q'(-psi)^T eta^-1 Q'(-psi'))/(psi + psi')`,
`E'(m, n) = [psi^m psi'^n] E'`; leaf `T'_i(z) = z (Q'^(i)_{0,0} - Q'^(i)_0(-z))`, `T'^(i)_d = (-1)^d Q'^(i)_{0,d-1}` (`d >= 2`; the `d = 1` term vanishes).
**Hatted (G).** `Yhat^mod := Lambda' Y^C(v(u), s^4 zhat) = u tau^(-1/2) y(tau zhat)` with `y_b(t) := Y^C_b(1, t)` (homogeneity `Y^C_b(v, z) =
sigma^(c_b+1) y_b(z/sigma^4)`, Lemma V), and `Yhat^Q := Lambda Y^C_Q(u, s^4 zhat)`, `Yhat^Q_a = Sum_n u^((4 - c_a - 4n)/5) zhat^-n y^(n)_a(u)`. Applying
`Lambda` to (G) gives **`Ghat(u, zhat) Yhat^mod = Yhat^Q`** (the same terms regrouped: a slot `(a, z^m u^n)` goes to `(a, zhat^m, u^(n + (4m+4-c_a)/5))`
bijectively). By Lemma Q(a) every term of `Yhat^Q` has u-order exactly 1, so `Yhat^Q = u x (series in u with coefficients in Q[[zhat^-1]])`:
**no positive power of `zhat`**. At each u-order, `Ghat` is a polynomial and `Yhat^mod` a series in `zhat^-1`, so all products are finite sums.

## 2. The vanishing covector

**Definition.** `f(x) := Psi0^-1 Yhat^mod(x) = u tau^(-1/2) Psi0^-1 y(tau x)` (colour components), `phi_{i,m} := (-1)^m [x^-m] f_i(x)`, `m >= 0`;
`(Psi0^-1)[i][b] = (1/5) eps_b zeta^(-i(b+1))`. Then `Q'(x) f(x) = Ghat(u, x) Yhat^mod(x) = Yhat^Q(x)`.
**Closed form (PROVED).** `y_b(t) = eps_b Sum_i (-1/(4t))^i/i! (-1/t)^(c_b+4i+1)/(c_b+4i+1)!`, so `[t^-m] y_b != 0` only for `m = 5i + c_b + 1`, i.e.
`b = 1 - m (mod 5)` (one class per `m`), with `i = floor(m/5) =: q` and exponent `j = m - q`; then `(-1)^m eps_b [t^-m] y_b = (-1)^(m+q+j) /(4^q q! j!) =
1/(4^q q! (m-q)!)` (`q + j = m`) and `zeta^(-i(b+1)) = zeta^(i(m-2))`. Hence
        `phi_{i,m} = (u tau^(-1/2)/5) tau^(-m) zeta^(i(m-2)) / (4^q q! (m-q)!)`,   `q = floor(m/5)`.
(VERIFIED `m <= 80`, `gi_formulas.log` (F).) Every `phi_{i,m}` is `O(u)` exactly; there is NO decay in `m` (finiteness comes from the tables, s.4).
**Fourier/flat form.** As a function of the colour, `phi_{.,m}` is the character `zeta^(i(m-2))`: one Fourier mode `r = 2 - m (mod 5)` in the
convention `phihat(r) = (1/5) Sum_i zeta^(ri) phi_i`. On flat insertions (colour components `5 Psi(u)^-1 e_{3-c}` of `Hhat^c`; `Sum_i Psi0^-1[i][3-c]
Psi0^-1[i][b] = (1/5) eps_{3-c} eps_b delta_{bc}`, `eps_{3-c} eps_c = 1`): **`phi(tau_m(Hhat^c)) = (-1)^m Ghat(u,0)_cc [zhat^-m] Yhat^mod_c`**, the vanishing
solution itself in the natural pairing (non-zero only for `c = 1 - m mod 5`; `Ghat(u,0)_cc = tau^((1-c_c)/4) p_{0,c}(u) = 1 + O(u)`, `c_c` the hatted exponent of the class `c`). At the frozen point:
`tau_0(H) -> -u tau^(-1/2)`, `tau_1(1) -> u tau^(-3/2)`, `tau_2(H^4) -> u tau^(-5/2)/2`, `tau_3(H^3) -> u tau^(-7/2)/6`, `tau_4(H^2) -> -u tau^(-9/2)/24`,
`tau_5(H) -> -u tau^(-11/2)/96`, ..., the terms of `Y^C = (-v/z, -1, -v^4/(24 z^4), -v^3/(6 z^3), v^2/(2 z^2)) + ...` with `z^-m -> (-1)^m` (the plan's
`<tau_2(H^4)> ~ v^2/2`). On the hatted unit at `u` (`c = 0`, `Ghat(u,0)_00 = 1/xhat`): `phi(tau_0(1hat)) = 0`, `phi(tau_1(1hat)) = u tau^(-3/2)/xhat = vhat`.

## 3. The contraction lemma (L)

**Lemma L (PROVED).** Let `Q'(x)` be rows with `Sum_ab Q'^(i)_a(x) eta^ab Q'^(j)_b(-x) = delta_ij` (symplectic), and `phi = (phi_{i,m})` any covector with
`f_i(x) = Sum_m (-1)^m phi_{i,m} x^-m`. Then for all `n >= 0`, `j`:
        `Sum_{i,m} phi_{i,m} E'(m, n)_ij = [psi'^(n+1)] Sum_ab h_a(psi') eta^ab Q'^(j)_b(-psi')`,   `h(x) := [Q'(x) f(x)]_(>0)`,
equivalently the contraction is `(1/psi') Q'(psi')^-1 h(psi')` (since `Q'(-x)^T eta^-1 = Q'(x)^-1`).
*Proof.* Symplecticity at `x = psi'` gives `1 = Q'(psi')^T eta^-1 Q'(-psi')`, so `E'_ij = Sum_ab D^(i)_a eta^ab Q'^(j)_b(-psi')` with the difference quotient
`D^(i)_a = (Q'^(i)_a(psi') - Q'^(i)_a(-psi))/(psi' + psi) = Sum_N q^(i)_{a,N} Sum_{k=0}^{N-1} psi'^k (-psi)^(N-1-k)`, `Q'^(i)_a(x) = Sum_N q^(i)_{a,N} x^N`.
Contracting `psi^m` with `phi_{i,m}` picks `m = N-1-k` with the sign `(-1)^m`: `Sum_{i,m} phi_{i,m} [psi^m] D^(i)_a = Sum_{k>=0} psi'^k Sum_i Sum_m (-1)^m
phi_{i,m} q^(i)_{a,m+k+1} = Sum_{k>=0} psi'^k [x^(k+1)] (Q'(x) f(x))_a = h_a(psi')/psi'`. QED
**Where the sign comes from.** The edge is built from `R'(-psi)`; the difference quotient in `(psi', -psi)` produces `(-psi)^m`. Without the
sign the statement fails (control below). **Finiteness.** At each u-order both sides are finite sums: `E'(m, n)` vanishes at `u^N` once
`m + n + 1` exceeds the two-factor degree bound of `[u^N] Ghat` (Theorem D), and `h` is a polynomial times a series in `x^-1`.
**Interpretation (QUOTED formalism, not used in the proofs).** In the quantization formalism of arXiv:math/0108100 a covector on the descendant
insertions is the quantized linear Hamiltonian of an element `f` of the negative half of the loop space, and `R^ f^ R^-1 = (Rf)^`; the creation
part of `(Rf)^` is `(Rf)_+`. Lemma L is the edge-level form of this covariance: an edge turns a covector at one end into the insertion
`(1/psi') R'^-1 (R'f)_+` at the other.

## 4. The edge statement (E) and the leaf statement (T)

**(E) (PROVED given Theorem U).** For the canonical gauge, `R'` is symplectic (Proposition S; `Psi0` orthonormal) and, with the vanishing
covector, `Q'(x) f(x) = Yhat^Q(x)` has no positive power of `x` (hatted (G)), so `h = 0` in Lemma L: **for every `n`, `j` and every u-order,
`Sum_{i,m} phi_{i,m} E'(m, n)_ij = 0`.** The statement is about the edge factor alone and therefore holds in every graph of every genus.
**(T) (PROVED given Theorem U).** `Sum_{i, d>=2} phi_{i,d} T'^(i)_d = Sum_i Sum_{d>=2} [x^-d] f_i . [x^(d-1)] Q'^(i)_0 = [x^-1] (Q'(x) f(x))_0 - Sum_i Q'^(i)_{0,0} [x^-1] f_i`
(the `d = 0` term vanishes: `Q'` has no negative powers). Now `[x^-1] (Q' f)_0 = [zhat^-1] Yhat^Q_0 = y^(1)_0 = -t_c` (Lemma Q(d)), `Q'^(i)_{0,0} =
Ghat(u,0)_00 Psi0[0][i] = zeta^i/xhat` (the diagonal formula `p_{0,0} = tau^(-1/4)/varpi_0`, GAUGE_IDENTITY s.5), and `[x^-1] f_i = -u tau^(-3/2) zeta^-i/5`
(only `b = 0` has a `t^-1` term, `y_0 = -1/t + ...`). Hence **`Pi(T') = vhat - t_c`, `vhat = u tau^(-3/2)/xhat`**, in every conifold frame. Leading
term: `T'_2` is `u p_{1,4}(u)` in every colour (the direction `tau_2(Hhat^4)`) and `phi_{i,2} = u tau^(-5/2)/10`, so `Pi(T') = (p_{1,4}(0)/2) u^2 + O(u^3)
= (2/3 - gamma_t - w_1) u^2 + ...`; frame (0,0): `-u^2/30`. (UNIFORM_REDUCTION s.6.2 has `-p_{1,4}(0) = 2 gamma_t + 2 w_1 - 4/3`.)
**Actual versus canonical.** (E), (T) are theorems for the canonical gauge; for the production `R'` they hold at every order where the actual
gauge is canonical (VERIFIED through `zhat^24`, `u^17`; GAUGE_IDENTITY s.6).

## 5. What (E) and (T) give, and what they do not

* **Edges.** If the vertex correlator at the end of an edge depends on that half-edge only through `phi`, the edge term vanishes at all
  u-orders. **Leaves.** If a vertex depends on its leaf insertions only through `phi`, the leaves enter only through `Pi(T') = vhat - t_c`.
* **The sign is forced by the dilaton equation.** A one-point law `<xi>_h = kappa phi(xi) F_h'(vhat)` (singular part, `F_h = c_h vhat^(2-2h)`) holds for
  `xi = tau_1(1hat)` exactly, with `kappa = -1`: `<tau_1(1)>_h = (2h-2) F_h = -vhat F_h'(vhat)` and `phi(tau_1(1hat)) = vhat`. With `kappa = -1` the leaves shift
  the argument `vhat -> vhat - Pi(T') = t_c` (at the linear level), and one leaf gives `-Pi(T') F_g'(vhat) = (t_c - vhat)(2-2g) c_g vhat^(1-2g) =
  -(g-1)(c_g/15) u^(3-2g) + ...` in frame (0,0), i.e. `-1/18000`, `-1/189000` at genus 2, 3: **the single-leaf bookkeeping on record (NATURAL_FORM s.2,
  UNIFORM_REDUCTION s.6.2) follows from (T) GIVEN the one-point law for `tau_2(H^4)`** (the only insertion at this order), including its frame
  dependence `-(g-1)(2 gamma_t + 2 w_1 - 4/3) c_g`. The one-point law itself is a statement about the model theory (plan, Stage 1): CONJECTURAL here.
* **Dimension caveat (PROVED counting).** `phi_{i,d} != 0` for every `d`, but one-point descendants of a genus-h CohFT vanish for `d > 3h - 2`, so
  no one-point law can hold for all `d`. It is needed only where it reaches the polar part: a leaf `T'_d` is `O(u^mu_0(d-1))`, `mu_0(k) >= (4k-1)/5`,
  and a genus-g term is polar only if `mu_0(d-1) <= 2g - 3`, i.e. `d <= (5g-5)/2 < 3g - 2`.
* **(E) is not sufficient by itself.** The covector is LINEAR; the model's singular structure is not: by the string equation
  `<tau_1(H) tau_0(1)>_h = <tau_0(H)>_h` is singular although `phi(tau_1(H)) = phi(tau_0(1)) = 0` (support rule `c = 1 - m`). The mechanism needs the
  non-linear (genus-0 descendant) extension of the vanishing period and its compatibility with the edges (plan, Stage 1-3). Not addressed here.
* **Linear-response regression (VERIFIED, `gi_hat_table.log` (A)).** Of the 15 single-term left dressings on record
  (`natural_form/extra/nat_kregular_linear.log`, statuses re-read today), every one that annihilates the covector is invisible on record (4 of
  4); the converse fails once (`zhat^6 (H^4 (x) H^3 - H^3 (x) H^4)` is invisible but not annihilating): 14 of 15 agree. The Euler-dilaton block
  `alpha zhat (H - zhat)(x)(H + zhat)` of `G_2` annihilates the covector ((R3) in the hatted Fourier convention); the opposite relative sign, visible
  on record, does not.

## 6. Verification (all logs re-produced today)

| check | rows | result |
|---|---|---|
| hatted (G) on the referee's G table over Q: positive part (84 graded coefficients, 73 genuine sums), non-positive part = `Yhat^Q` from the periods (264), leaf over Q | `natural_form/referee/G_fourier_NZ16_K12.json`, frame (0,0), `zhat <= 16`, `u <= 11` | ALL PASS (`gi_hat_table.log`) |
| (L) with 3 random covectors (300 series identities); degree bound of the edge tables | referee rows, `zhat <= 24`, `u < 18`; p89, p61; frames (0,0), (1/3,-2/7) | ALL PASS |
| (E): 2160 coefficients (`n < 24`, all `j`, `u^N`, `N <= 17`), 965 of them sums of >= 2 non-zero entries | same | ALL PASS; control without `(-1)^m`: 880 of 2160 non-zero |
| (T) mod `u^18`; unit row `zeta^i/xhat`; leading term `(p_{1,4}(0)/2) u^2` | same | ALL PASS |
| (L), (E), (T) on the production rows (`fq3_par.setup`, symbol `M_0`) | `zhat <= 6`, `u < 6`, p89, p61, both frames: 150 edge coefficients, 65 genuine | ALL PASS |

Logs: `gi_hat_edge_ref_{89,61}_ab_{0_0,1o3_m2o7}.log` (run with `NZ = 24`, `K = 18`), `gi_hat_edge_prod_*.log`, `gi_hat_table.log`, `gi_formulas.log`.
Scripts `gi_hat_core.py`, `gi_hat_table.py`, `gi_hat_edge.py` were written by the first attempt and reviewed here line by line (edge and leaf
tables are the production class `fq_sser.EdgeLeaf`; `R'` rows `Q' = Qhat M_0(tau zhat)^-1`; the covector without its common factor `u tau^(-1/2)`).
