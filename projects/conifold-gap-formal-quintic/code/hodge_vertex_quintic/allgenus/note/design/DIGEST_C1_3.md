# DIGEST_C1 (piece 3 of 6): s.2 Statements, parts MOD (model singular structure), TR (TR and x-y swap on the model), UNI (universality)

Same format as piece 2. Local identifiers MOD-n, TR-n, UNI-n (the letters A-C of piece 2 and D1-D19 of piece 1 are taken). Label
prefixes: `mod` (model singular structure; remarks in s.5-6 of the plan), `trm` (model as TR, s.5), `mlk` (L-regularity and the Gap
Theorem, s.6), `uni` (universality, s.9).

## MOD. The singular structure of the model [s.5-s.6 remarks; evidence appendix]

**MOD-1 `prop:mod-onepoint` (prop).** For `h >= 2` and `X in K` homogeneous with `<X>_h` singular: `<X>_h = Pi_0(X) F_h'(v)` (the n = 1
case of ML_K). Off `K` this "one-point law" fails: `<tau_5(H)>_3 = 53/1209600 = (-53/192) F_3'` (e = -1) against the law's
`(1/96) F_3' = -1/604800`; at genus 4 the remainders of `tau_5(H)`, `tau_6(1)`, `tau_7(H^4)` are `(19/100800000) x (1, -1, -1/2) =
(19/100800000) Pi_5(X)`, carried by the component along `t^H_5`. Regular one-point functions need not obey the law (e.g. `<tau_3(H^3)>_2`).
At genus 1 the law holds for `H` and fails for the dilaton by the regular anomaly `1/8`. Source MS1 s.0 item 2, s.3; MRF item 2 (exact
over Q, confirmed); TRM s.3 (one-point lists); TRRF claim 14 (robust to sign and to `z -> mu z`). Status: on K PROVED as an instance of
ML_K (final, S2A s.1); VERIFIED `h = 2, 3, 4`; off K REFUTED (exact counterexample). Use: remark "the restriction to K is essential",
numbers in the evidence appendix.

**MOD-2 `rmk:mod-noV` (rmk).** No genus-independent non-linear `V` exists with singular part `F_h(V(q))` on all of `Hc_+`:
`Hess_h(X, Y) := (<XY>_h - F_h'' PiC(X) PiC(Y))/F_h'` depends on `h` for `t1(H) t2(H^3)` (`-23/12` at `h = 2`, `-13/48` at `h = 3`),
`t2(H) t1(H^3)` (`-11/6`, `-31/60`), `t3(H) t0(H^3)` (`-3/2`, `-11/60`), all four correlators singular; each pair has one insertion off K;
42 singular pair instances differ from the string-only extension `V_str`; the ancestor genus-0 point `W` is excluded by
`<tau_2(H^4)>_2 = c_2/v != 0`; the descendant genus-0 point has dilaton degree 0. All genus dependence sits in the D-directions.
Source MS1 s.0 item 3, s.4; MRF item 3; RD3 s.4.2. Status REFUTED (exact). Corrections: MS1's "rules out ... any jet expression" holds
only for genus-INDEPENDENT `V` (jets of the genus-0 solution included, also non-homogeneous `V`); genus-mixing structures (hbar-dependent
`V`, Toda corrections) are not excluded (MRF item 3; SYN01 s.2 "Model", s.7). The plan's non-linear `V(q)` is REFUTED and not needed
(SYN01 s.0). Use: remark plus evidence.

**MOD-3 `prop:mod-H2tower` (prop; superseded).** On `Q tau_0(H) + span{kappa(H^2, k)}` ML_K holds in all genera: by arXiv:2211.12259
Thm 2 (special point `f = p^(N-1) + 1/p`, calibration (1.24), `N = 5`; dictionary: `t^1 = a_1` <-> `H^2`, tower `kappa(H^2, k)`,
special point `v = 1/4`, `a_5 = 4v`), a genus-h correlator of `kappa(H^2, k_i)` has `e = 2 - 2h + 4B - n` with `B = Sum (k_i + 1)/5`,
which by Euler's formula `V - 5B + (n + B) = 2 - 2h` is the number `V >= 1` of vertices of the hypermaps counted; adding `m` insertions
`tau_0(H)` gives `v^(V-m) V!/(V-m)!`. Source MS2 s.8; MRF item 5 (Thm 2 and its proof read, pp. 7-10; precision: the homogeneous equation
allows `c H` in `S_1 H^2` with `-c log v H^2` in `S_2 H^2`; log-freeness of the `beta = 1` column of (1.24) forces `c = 0`). Status
PROVED given the quoted theorem and the dictionary. Superseded by ML_K for all of K (S2A s.1). Use: optional remark. MS2 s.8 names
the authors: violates the citation rule (MRF s.4).

**MOD-4 `rmk:mod-toda` (rmk).** `Sum_h hbar^(h-1) F_h(V)` (`F_0`, `F_1` included) is the asymptotic expansion of `log G(1 + N)`
(Barnes `G`) at `N = V (-5/hbar)^(1/2)`, since `c_h hbar^(h-1) = B_2h/(2h(2h-2)) (-hbar/5)^(h-1)`; from `G(N+2) G(N)/G(N+1)^2 = N`
the singular part on K solves a second-difference (Toda) equation in `Pi_0` with step `(-hbar/5)^(1/2)` (formal, elementary).
Source MS2 s.9(b); RD3 s.5.1 (ML-Toda, CONJECTURAL as a route); XY s.5 (Toda for `log G(1 + v/hbar)`). Correction: for the TR/CohFT
series the right-hand side is `log(v/hbar)` only up to a constant fixed by the genus-0 quadratic terms (XYRF s.5; SYN01 s.2 "xyswap").
Not used.

**MOD-5 `prop:mod-primary` (prop).** (Primary-sector lemma.) For `h >= 1` every primary correlator `<tau_0(H^(a_1)) .. tau_0(H^(a_n))>_h`
vanishes unless all `a_i = 1`, and `<H^n>_h = F_h^(n)`. Proof route: `tau_0(H^a) = kappa(H^a, 0) in K_0` for `a != 1` and
`w(0, a) = (c_a - 4)/5 < 0` (model sign: `(4 - c_a)/5 > 0`), so every primary correlator at `h >= 1` is singular by weight. Source MS2
s.9(c); MRF item 9; RD1 s.2.1(d). Second route for `h >= 2`: TR-5. Status PROVED given ML_K. Not used separately by the chain.

**MOD-6 `rmk:mod-delta` (rmk; first-order normalisation effect on totals).** For `P = P_can(1 + c z^m + ..)` (`m` odd; leaf `c z^(m+1) 1`,
edge `c Sum_j (-1)^j w^(m-1-j) z^j eta^-1`) the genus-g total moves at first order by `c Delta_g(m)`, `Delta_g(m) = <tau_(m+1)(1)>_g +
(1/2) Sum_j (-1)^j [<tau_(m-1-j)(e_a) tau_j(e^a)>_(g-1) + Sum_(h1+h2=g) <tau_(m-1-j)(e_a)>_(h1) <tau_j(e^a)>_(h2)]` (frozen values).
`Delta_4(5) = -1/8400000` (leaf) `+ 0` (loop) `+ 1/8400000` (bridges) `= 0`; `Delta_3(5) = -1/640` (exponent 0, regular);
`Delta_2(5) = 0`; `m = 1, 3, 7` vanish by the charge rule. Source MS2 s.9; MRF item 10 (formula re-derived, values confirmed).
Correction (SYN01 C1, s.7): "genus 4 is NOT a sharp test of the `z^5` normalisation" holds for the TOTAL only; per graph the response
is non-zero (vertex `-1/8400000`, bridge (1,3) `-1/1008000`, bridge (2,2) `1/900000`, loop 0 per unit `c_5`), so genus 4 IS a sharp
per-graph test (univ is right). Status VERIFIED (exact). Use: evidence appendix; supports UNI-7 (graph-locality).

## TR. The model as topological recursion; the x-y swap [s.5-s.6]

**TR-1 `thm:trm-identification` (thm; S2A "M1").** For `x = W_v = p^4/4 + v/p`, `y = p`, `B = dp1 dp2/(p1 - p2)^2`:
`omega_(h,n) = (-5)^(h-1) Sum <prod_i tau_(d_i)(H^(a_i))>^frozen_h prod_i d xif_(a_i, d_i)`, frozen = the ancestor correlators of
`M_0.omega^` at `v`; for `x = kappa W`: `omega_(g,n) = Lambda^(g-1) Sum kappa^-(3g-3+n-Sum d) <..>^frozen prod d xi`, `Lambda = -5 kappa`
(for `kappa = 1` the TR CohFT IS the frozen theory with metric `x (-1/5)`). In all genera `F_h^frozen = (-5)^(1-h) F_h^TR`. Source TRM
s.0 item 2, s.1, s.2.2-s.2.4; TRRF s.3.1, s.3.3; S2A s.1 M1. Status PROVED modulo arXiv:1509.06954 Thm 8.1 (hypotheses checked: TRM
s.1; TRRF claim 1: conditions (1)-(4) of Thm 5.1, `d_lambda p` without zeros over finite `lambda`, distinct critical values `5 p_c^4`) and
the uniqueness of the homogeneous R-matrix (TR-2); `F_h` by the dilaton equation on both sides (dilaton leg = LD at `(b,k) = (0,1)`,
S2A M1; or TR-3). Used by the chain (step 7.1). Corrections (TRRF s.6): "x = -W flips psi -> -psi" holds only up to `(-1)^(3g-3+n)`;
state the `kappa^-D` formula as the claim; "any point" -> "any point where `x` is simply ramified"; the metric-sign remark on (2-7) is
incomplete: `TR(lambda) = (-1)^(g-1+D) x` (Hurwitz CohFT in the LG convention), harmless; author names to be replaced by citations.

**TR-2 `prop:trm-Rmatrix` (prop).** With `X(z)[j][i]` := the Laplace transform (2-31) (row `j` = branch point, column `i` = xi-function):
`X(z) = M_0(z/kappa)^-1 = M_0(-z/kappa)^T`, i.e. `R_TR(z) = M_0(z/kappa)`, to all orders; (2-32) gives `M_0(z/kappa)^-1 1^`
(`1^ = Sum_k et_k/c_k`, the frozen unit), so `omega_(0,1) = y dx` produces the frozen dilaton leaf `T(z) = z(1^ - M_0^-1(z) 1^)`
(`kappa = 1`; colour-0 coefficients `z^2..z^5`: `-3/20, -9/800, 11571/80000, 2415609/6400000`). The transposed reading fails from `z^2`.
All orders: both are homogeneous R-matrices of the same `(U, V)`; `[R_(m+1), U] = (m + V) R_m` determines `R` from `R_0 = 1`; `U` agrees
(critical values) and `V = [R_1, U]` agrees because `R_1` agrees. Source TRM s.2.3; TRRF claim 2, correction 4. Status VERIFIED through
`z^8` (worker), `z^10` (referee); all orders "PROVED modulo that standard uniqueness" (TRM s.2.3, s.8). TRRF correction 4: name the inputs
(homogeneity of `R_TR`: arXiv:1509.06954 Prop. 5.2 under Thm 8.1; `M_0` satisfies `(5/4)[C_(n+1), Z] = (V - n) C_n`, LG s.3-4); a
direct route: `xit^i x' = -(1/c_i) Sum_k c_i^k p^(2-k)`, so the Laplace transform of `d xi^i` is an exact combination of the
stationary-phase series `J_k`, `k = -2..2`, of LG s.1. Used by the chain (inside M1).

**TR-3 `lem:trm-dilaton` (lem).** `P_(i,d) := Sum_j Res_(c_j) Phi d xit^i_d = -delta_(d,1)/c_i` for every `d >= 0`; hence
`F_g^TR = Lambda^(g-1) kappa^(3-3g) F_g^frozen` in all genera. Source TRRF s.3.3 (closes the worker's gap: TRM checked `d <= 3`).
Status PROVED; VERIFIED `d <= 10`, `kappa = 1, 4, -1`, both primes. Alternative in the chain: LD at `(b, k) = (0, 1)` (S2A M1, step 7.2).

**TR-4 `thm:mlk-gap` (thm; "Gap Theorem", TRM's "Theorem G").** Let `x = P(z) + v/(z - a)`, `P` a polynomial of degree `r >= 1`,
`v != 0`, `x` simply ramified, `y = z`, `B = dz1 dz2/(z1 - z2)^2`. Then for every `g >= 2`, `F_g = B_2g/(2g (2g-2)) v^(2-2g)`, in the
convention of arXiv:1509.06954 (2-28) = arXiv:2502.18115 (2.1) = the engines (`F_2 = -1/240` for `x = z + 1/z`; see s.4 CV9). For
`x = W_v`: frozen `F_g = (-5)^(1-g) F_g^TR = c_g v^(2-2g)`. Source TRM s.4.1; TRRF s.3.2 (re-derived; b-points of arXiv:2502.18115 p. 13
are exactly `a` and `infinity`); XYRF s.2 ("short proof from (1.3)"); S2A M8. Status PROVED by two routes: (1) mod arXiv:2502.18115
Cor. 1.2, eq. (1.3); (2) mod arXiv:2212.00320 Thm 7.1 plus a Borel sum (TR-8). VERIFIED (s.5). Used by the chain (step 7.6, ML_K(ii));
route 1 is the optional quote of (Q). Corrections: rename "Theorem G" -> "Gap Theorem" (SYN01 C9); TRM's label "EO/DNOPS conventions":
the kernel is MINUS the printed kernel (4-3) of arXiv:math-ph/0702045 (S2B K1); arXiv:2502.18115's printed closed forms (s.3.2, Cor. 3.3,
Prop. 3.4) carry a sign slip, `B_2g/(2g(2-2g))` instead of `B_2g/(2g(2g-2))`, from "`d^(m-1)z/dx^(m-1) = (-1)^m (m-1)! z^m`" (correct:
`(-1)^(m-1)`) (TRM s.4.1 "Sign remark"; TRRF claim 9).

**TR-5 `cor:trm-residue` (cor).** On the whole Hurwitz space `H_(0;1,r)` (`r >= 1`, at points where `x` is simply ramified; the set is
dense and connected and `F_g` extends) `F_g`, `g >= 2`, depends only on the residue of `x` at its simple pole. Hence the primary-sector
lemma for `g >= 2` (cf. MOD-5). Source TRM s.0 item 3, s.4.3; XY s.0 item 3, s.2.3 (whole family); TRRF claim 8, correction 5. Status
PROVED mod arXiv:2502.18115 Cor. 1.2 (resp. arXiv:2212.00320 Thm 7.1); the primary-sector consequence also uses arXiv:1509.06954
Thm 8.1 and "the CohFTs over the Frobenius manifold form the shifted family" (standard, not re-proved; TRM s.8; genus 1 not covered).
Not used by the chain (ML_K contains it).

**TR-6 `rmk:trm-v0` (rmk).** Homogeneity: `p -> k p`, `u -> k^5 u` gives `omega -> k^(5(2-2g-n)) omega`, hence for `|p| > |u|^(1/5)`
`omega^(u)_(g,n)(p) = Sum_(K>=0) u^(2-2g-n+K) Omega_K(p)` (PROVED; TRRF s.3.4). At fixed `p != 0` every computed `omega_(g,n)`, `n >= 1`,
is a power series in `u` with a non-zero limit: VERIFIED for 13 `(g,n)` (8 re-derived, TRRF claim 11); all `(g,n)` CONJECTURAL; limits =
4-spin correlators CONJECTURAL. Source TRM s.0 item 4, s.5. Superseded for L-legs by TR-9 (SYN01 E9: "TR s.0 item 4 understates").
Not used. TRM s.6 (Delta B deformation, contour argument) is SKETCH and must not enter the note.

**TR-7 `thm:mlk-closed` (thm; applicability of the closed formula).** arXiv:2212.00320 Thm 7.1 applies to `x = z^4 + v/z` (its
example (7.2), `r = 5`, after `z -> v^(1/5) z`) and to the whole family `x = z^4 + c3 z^3 + c2 z^2 + c1 z + c0 + V/z`: `dx` has simple
zeros (`V != 0`), `y = z` is regular with `dy != 0`, `x` rational. Consequently `omega_(g,n) = N_(g,n)(z; v) prod dz_i /
prod (4 z_i^5 - v)^(6g-4+2n)` with `N` polynomial (family: `q(z) = z^2 x'(z)` in place of `4 z^5 - v`); kernel dictionary
`omega^code_(g,n) = (-1)^n omega^ABDKS_(g,n)`. Source XY s.0 items 1-2, s.1, s.2.1; XYRF s.2 (statement read: s.1.1-1.5, s.7; dual
hypotheses hold; the kernel relation VERIFIED with own code). Status: applicability CONFIRMED; the structure PROVED from (7.9).
Correction: Remark 7.3 of the paper speaks of the discriminant; `v = 0` is a boundary point where the pole disappears (harmless, XYRF
s.2). Used by the chain (the hypothesis check of quote (Q) arXiv:2212.00320 Thm 7.1).

**TR-8 `prop:mlk-period` (prop; Borel route to the gap).** On the whole family, `Int_inf^0 omega^code_(g,1) = (B_2g/2g) V^(1-2g)`
(`g >= 1`), and `F_g = (1/(2-2g)) (-V) Int_inf^0 omega^code_(g,1) = B_2g/(2g(2g-2)) V^(2-2g)` (`g >= 2`). Steps: monomial count near 0 and
infinity; the primitives `T_r` regular at 0 and `Phi_g(inf) = 0`; the pure-pole constant term; Borel summation with
`eps s = tanh(eps u)`: `Sum_r (r-1)! [t^(r+1)] e^(E_0) = Int_0^inf e^-u eps^2/sinh^2(eps u) du`; the polynomial part of `Phi` pairs to 0 by
weighted degree; the `-V log z` part by the contour lemma `Sum_a Res_a(log z omega) = Phi(0) - Phi(inf)`. Source XY s.2.3; XYRF s.2
(confirmed; two-line proof of the contour lemma; Borel identity VERIFIED `g = 1..10`). Status PROVED mod Thm 7.1. Used: route 2 of TR-4
(makes arXiv:2502.18115 optional, S2A M8).

**TR-9 `thm:mlk-Lreg` (thm; "L-regularity", XY's "ML1", S2A "M7", SYN01 "E9").** For `2g - 2 + n > 0`, every correlator with all `n >= 1`
legs in `L` is a polynomial in the parameters (on the whole family), in particular regular at `V = 0`; on the line: legs in `Q[v][p]`
give a pairing in `Q[v]`. Source XY s.0 item 4, s.3; XYRF s.3; S2A M7. Status PROVED mod arXiv:2212.00320 Thm 7.1; VERIFIED (XYRF run B
off the line: 12 points, 36 L-correlators polynomial in `S`, controls singular). Used by the chain (step 7.5). Correction: XY's
"11 L-correlators that homogeneity would allow to be singular vanish" are forced by decay at infinity and are not evidence (XYRF s.3;
SYN01 s.7). Name clash: XY's ML1/ML2 (L-regularity / gap) vs MS1's ML1/ML2 (pure-H sector / kappa-regularity) (SYN01 C9).

**TR-10 `rmk:mlk-legs` (rmk).** In XY's leg dictionary (`tau_0(1) <-> -z`, `tau_0(H) <-> -log z`, `tau_k(a) <-> Int tau_(k-1)(a) dx`):
`tau_1(H)` is a logarithmic leg (monodromy around `|z| = eps`), outside `L`; the string equation reads `<omega_(g,2), -z, Int z^m dx> =
+<omega_(g,1), z^m>` (unit acts by `+d/dx`); `<tau_0(1) tau_1(H)>_g = <tau_0(H)>_g` = `-(1/12) v^-1` (`g = 1`), `(1/120) v^-3` (`g = 2`).
Euler-dilaton `e <-> -(4/5) z^5` (log-free part of `-Phi`) in `L`, with `Sum_a Res_a (4/5) z^5 omega_(g,n+1)(z, J) = (2 - 2g - n -
v d_v|_x) omega_(g,n)(J)`; `<e>_1 = 1/8`, `<e, H>_g = 0`; the Euler-dilaton edge (UR s.6.3(ii)) has half-edge legs `-z`, `-(4/5) z^5`,
`-(16/45) z^9 + (v/5) z^4`, all in `L`; the chain `psi^k e` is log-free for `k <= 4` and has `log z` at `k = 5`. Source XY s.0 items 5-6,
s.3-s.4; XYRF s.4. Corrections: string sign `+d/dx`, XY's prose had `-d/dx` twice, errors cancelling (XYRF s.4; SYN01 C8); coefficient
`-16/45` not `-2/5` (XYRF s.4; SYN01 C7); this dictionary is NOT the one of record: LD (S2A step 7.2) is (S2B K10). "`Y^C` and the TR
monodromy are different functionals" is REFUTED (SYN01 C5). "Mixed correlators are singular in general" WEAKENED (`(z^5, -log z)` gives 0;
XYRF s.4). Not used by the chain.

**TR-11 `rmk:mlk-F1` (rmk).** `F_1 = -(1/12) log v` on the line from the EO variational formula (`d_v(y dx)|_x = -dz/z`, `d_v F_1 =
Int_0^inf omega_(1,1) = -1/(12 v)`). Source XY s.2.3. Status WEAKENED: the orientation is fixed by consistency with the dilaton
definition at `g >= 2`, `F_1` was not computed from its definition, off the line not treated (XYRF s.2). Not used as such (the chain's
`h = 1` base of Theorem V uses homogeneity and the value `-1/12`, S2A M9; other reader).

**TR-12 `rmk:mlk-ML3` (rmk).** A leg with poles at `z = 0` decomposes as `Phi_L + log z . f(x) + R`; the singular part of every correlator is
a polynomial combination of Gaussian vanishing periods and of `V^-1` from rewriting `1/z = (x - P(z))/V`. Source XY s.3 (ML3). Status
SKETCH. Superseded by D(ii), D(iii), LD and Theorem V (S2A s.1 M2-M4, M9). Omit. Likewise XY s.4 (M1)/(M2) (edges and leaves of the gauge
as legs) was not done; superseded by B2-B3 plus D(iii)/LD (S2A step 7; SYN01 s.2 "(M1)/(M2) settled conditionally").

## UNI. Universality, Lemma C, the failing controls [s.9]

**UNI-1 `cor:uni-flat` (cor; theory's "Corollary U", S2A "Corollary (universality, all genera)").** Let `A` be a self-dual cyclic
one-pole system (D15) with the product `B` of the pulled-back model, Theorem U applied to `A`, and residues `r = -1/lambda^2`, `rt = -1`.
Then (R1) gives column 1 of `P(0, z)` equal to `P_11(0, z) e_1`, (G) reads `P Y^C = P_11(0, z) Y^C_A`, and Prop. S becomes
`P(-z)^T eta P(z) = (rt/r) P_11(0, z) P_11(0, -z) eta` (the entry `(2,1)` carries `P_22(0, z) = (rt/r) P_11(0, z)`); so the constant
normalisation `P_11(0, z) = 1/lambda` is symplectic and annihilating, `Y' = Y^C_A/lambda`, `T_P = T_A/lambda`. With quintic Theorem R
and ML_K: every edge graph of the flat sum is regular and the vertex graph is `c_g (T_A/lambda)^(2-2g)` + regular, for every `g >= 2`.
Source RD2 s.3.3; TRF s.2 (CONFIRMED: `r P_22(0,z) = rt P_11(0,z)`, `t(z) = P_22(0,-z) P_11(0,z) = (rt/r) c(z) c(-z)`, `t = 1` for
`c = 1/lambda`); S2A s.2 Corollary. Status PROVED modulo (Q) and Theorem U for the system (flat form); the hatted-engine form used in
`univ/`, `bvar/` follows from S15 with `D = 1` by construction, but matching the engine's conventions to S15's Lemma 3 is SKETCH (S2A s.2).
Corrections: `lambda` (SYN01 C10): define it by `r = -1/lambda^2`, i.e. `lambda = f(0) |g(0)|`, not `lambda = f(0)` (RD2, UN1 s.1);
they agree when `g(0) = -1`; URF s.2.1: `c0 = 1/(f(0) |g(0)|)`; general bvar systems: `|c0| = m^10/|f(0) g(0)|`, sign fixed by
`sigma > 0`, the vertex law is even in `c0` (BVRF correction 3; S2B K6). Depends on: UR (R1), GI Theorem G and Prop. S for `r != rt`
(not read), C1, D14. Vertex law in the system's own coordinate `T = Int b_0 = Int_0^u f`, not in `t_c` (SYN01 s.5; UN1 s.2.2).

**UNI-2 `lem:uni-coordinate` (lem; "Lemma C" of BVRF, B is a coordinate).** Let `X = (f, g, B)` be a self-dual cyclic system with
`ord_0 B = -1`, model pull-back `v_X`, `a = v_X'(0)`, and `phi := v_X^-1 o v_q in w Q[[w]]`, `phi'(0) = 1/a`. (a) `X~ := phi^* X`
(`b~ = phi' . b^X o phi`) has the self-dual shape and `prod b~ = B_q`; it is a univ system with `F~ = phi' f(phi)/t_c'`,
`Gm~ = phi' g(phi)/g_q`, `F~(0) = f(0)/a`, `Gm~(0) = -g(0)/a`, model pull-back `v_X o phi = v_q`, hence the quintic `tau`. (b) A Theorem-U
gauge `P(u, z)` of `X` gives the gauge `P(phi(w), z)` of `X~` with the same `P_11(0,z)`, `P_00(0,z)`, `P_0(0)`: canonical goes to canonical,
stored to stored. (c) The pure-`1/z` solutions pull back: (G) for `X` iff for `X~`, `T''_(X~) = T''_X o phi`. (d) `R_(X~) = R_X o phi`,
`M_(X~) = M_X o phi`, `R'_(X~) = R'_X o phi`, and every outer-graph amplitude `V^(X~)_Gamma(w) = V^X_Gamma(phi(w))`. (e) A graph is
regular for `X` iff for `X~`; the vertex law transports. Hatted encoding: `G~(w) = (w/phi)^(2g-2) G(phi(w))` (powers of
`lam = (w/phi(w))^(1/5)`: edge `lam^(4(m+n+1))`, leaf `lam^(4e-3)`, vertex `lam^-(2g_v - 2 + n_v)`, total `lam^(10(g-1))`).
Source BVRF s.3; S2B K7 (re-derived: `prod b~ = phi'^5 B_X(phi) = -(v_X o phi)'^5/(v_X o phi) = -v_q'^5/v_q = B_q`). Status PROVED
(elementary); VERIFIED (`rf_reparam`). Consequence: B-universality is equivalent to `(F, Gm)`-universality plus coordinate covariance, not
a separate conjecture (S2B K7, O21, O22). Rename in the note (clash with RD2's Fock-space "Lemma C", B6).

**UNI-3 `rmk:uni-annihilation` (rmk; univ's "Conjecture A").** Statement (UN1 s.5): `P''(u, z)` symplectic, Z/5-graded, flat-regular,
`(P'' Y^C(v(u)))_+ = 0`; then in the hatted outer sum of `P''.Omega^LG_(v(u))` every edge graph is regular and the vertex graph is
`c_g T''^(2-2g)` + regular, `T'' = -[z^-1 e_0] P'' Y^C`. Status in UN1: CONJECTURAL. Its flat form is quintic Theorem R (a), (b) for
`P = P''` (RD1 s.0 item 1 and RD2 s.3.4: Theorem R needs no flatness, A1-A4 "regular, as predicted"); with ML_K it is PROVED mod (Q) in
flat form. The hatted form for non-gauge `P''` needs an S15-type identification that no source writes (open; s.7). Correction: "A +
(G) => Conjecture N" also needs (H) and `c^act = 1` (URF s.0 item 1; SYN01 s.2 "Univ"), now supplied by `D_0 = 1` (MECH5 s.2).

**UNI-4 `rmk:uni-stored` (rmk).** The stored gauge `P_stored = P_can (1 + m5 z^5 + ..)` is still a Theorem-U gauge between the same
two systems (flat), symplectic only through `z^9` (first defect at `z^10`: `c(z) c(-z) = 1 - (m5^2 - 2 m10) z^10 + ..`), and violates (A):
`P_stored Y^C = (1 + m5 z^5 + ..) Y^C_Q` with positive part `m5 (z^5 Y^C_Q)_+ = -m5 (z^5 e_1 + t_c z^4 e_0 + ..)`. Symplecticity fixes only
the even part of `log c`; (G) fixes the odd part (`c_5` etc.). So flatness is not sufficient; this holds BY CONSTRUCTION (Theorem U's
freedom `c(z^5)`), the experiment's content being the size and pattern of the failure (URF s.0 item 2). Source UN1 s.2.3 "Reading";
URF s.2.1; MECH1 s.1 item 2; SYN01 E2. Status PROVED (the algebra) / VERIFIED (the window).

**UNI-5 `lem:uni-valuation` (lem; evidence-level).** If the hatted rows change by the right scalar factor `1 + m u^4 zhat^5` modulo
`(zhat^10, u^6)` and `E_out = O(u)`, then `xhat` and `E_M` are unchanged, `E_R` and the leaves change by `O(u^4)`, an outer graph with `E`
edges gets `O(u^(E+3))` from one changed edge and `O(u^(E+4))` from one changed leaf, products of changes are `O(u^8)`; so the genus-4
difference is linear in `m`, confined to the vertex graph and graphs with `E <= 2`, at `[u^4], [u^5]` (`E <= 1`) and `[u^5]` (`E = 2`).
Source UN1 s.2.3; URF s.2.4 (CONFIRMED, PROVED elementary). Use: explanation in the evidence appendix (hatted engine).

**UNI-6 `lem:uni-delta6` (lem).** Let `D` be symplectic and `W := D Y^C_Q = Sum_m delta_m z^m e_(m+1)`. If `delta_m = 0` for all `m >= 1`
except `m = 6`, then `delta_6 = 0`: the `z^6` coefficient of `W(-z)^T eta W(z) = Y(-z)^T eta Y(z)` (no positive powers) is
`2 delta_6 delta_0 eta_21`, and `delta_0 = -1 + O(u)` is a unit; by induction the same for every self-dual even `m` (16, 26, ..). Source
URF s.3.4 (PROVED non-linearly; upgrades UN2 s.3.3 "PROVED first order / VERIFIED non-linear"). Use: consistency of the annihilation
conditions; remark or appendix.

**UNI-7 `rmk:uni-graphlocal` (rmk).** In the stored normalisation six genus-4 edge graphs (1, 2 at `[u^4]`; 3, 4, 7, 8 at `[u^5]`)
become singular while the total polar part is unchanged through `u^5` (the normalisation only moves polar terms between the vertex graph
and six edge graphs). Hence no argument about total potentials (KP/Toda of the full theory, total free energies, the gap of the total)
can prove Conjecture N; the proof must be graph-local and must use the odd part of `log c`, fixed by (G). Source MECH1 s.2; SYN01 s.5;
UN1 s.0 item 1, s.2.3; S2B s.6.1 (last bullet). Status: the observation VERIFIED; the conclusion is a methodological remark. Use: s.9
remark (or introduction), numbers in the appendix.

**UNI-8 `rmk:uni-necessity` (rmk; the failing controls as necessity examples).** None of the hypotheses of quintic Theorem R can be
dropped in general (examples only, evidence appendix): (A) violated: B1, B2 (UN2 s.3.2), the stored normalisation (UNI-4), the flat
`(1 + u^m z E_21) P` (TRF s.3 X3/X4), the wrong covector xA1 (URF s.5), the `c5`-dressed gauge breaking Lemma K at every shift (MS2 s.6);
flat regularity violated: C1a-C1c (UN2 s.3.2; C1c regular at genus 2 and at first order, singular at genus 3 through products:
first-order invisibility is not stable, since `mu(b,k1) + mu(b+k1,k2) = mu(b,k1+k2)`); at first order regularity is EQUIVALENT to
annihilation at genus 2, and the two genus-3 surplus directions are singular at genus 4 (UN2 s.3.1, s.4). For (Sp) there is no example:
without it the numerator of `E_R` does not vanish at `psi = -psi'` and the edge is not a power series (RD2 s.3.4; TRF s.0 item 3(v)).
Flatness (the u-ODE) is not a hypothesis: annihilating non-gauge dressings (A1-A4, `delta_0` not constant) are regular (UN2 s.3.2-s.3.4;
TRF s.3 `exp(u^m z N) P`). Also the restriction of ML to K is essential (MOD-1, MOD-2). Status: VERIFIED (finite) only; "flatness is neither
sufficient nor necessary" WEAKENED: "not sufficient" by construction (URF s.0 item 2; SYN01 s.2).

**UNI-9 `lem:uni-linear` (lem; optional).** The `eps^1` part of the hatted graph sum is linear over `u`-series in the dressing: the response
to `p_(k,b)(u) H(k,b)` is `p_(k,b)(u) u^mu(b,k)` times the response to the constant `h(k,b)` (every ingredient is polynomial in the row
entries; `xhat` is a rational function of the `zhat^0` unit row). Around the pure model (`R' = 1 + eps X`, `E_out = O(eps)`) only the
vertex graph and the one-edge graphs move at first order, so first-order Conjecture A around the model is a finite statement per genus
about the polar parts of one- and two-point descendant correlators contracted with flat-regular annihilating bivectors. Source UN2 s.3.1,
s.3.5 (both PROVED, elementary). Not used; evidence appendix at most.
