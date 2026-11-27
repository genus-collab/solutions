# Formal quintic, genus 2: scan of the constant covariant scalar frames (2026-11-04)

Folder: hodge_vertex_quintic/. Follows FORMAL_QUINTIC_genus2.md (same engine, same conventions, same conifold frame Cf with alpha = beta = 0). New code: fq_phi.py (the scalar family and its vertex classes), fq_scalar_lib.py, fq_scalar_report.py, driver genus2_scalar.py -> logs/genus2_scalar.log (both primes 2^89-1 and 2^61-1, 615 s), post-run checks genus2_scalar_post.py -> logs/genus2_scalar_post.log. No existing module was changed. Everything is an exact finite computation (F_p with rational reconstruction; the Q(zeta)-valued entries are reconstructed from their Galois conjugates). Evidence at genus 2 only; AI-checked only.

## 1. Set-up and formulas

**The family.** FORMAL_QUINTIC s.3 showed that in every Hodge frame tried, every coloured genus-2 graph with an edge has the maximal pole t_c^{-2} at the conifold point. The question here is whether any constant, symplectic, Z/5-covariant scalar factorisation of the R-matrix does better. Write R = R' phi_i(z) at colour i with

    phi_i(z) = phi_0(zeta^{-i} z),    phi_0(z) = exp(a_1 z + a_2 z^3),

a_k constants (q-independent). Only odd powers of z keep phi symplectic (phi_0(z) phi_0(-z) = 1), and only a_1, a_2 enter at genus 2 (total z-order of a genus-2 decoration <= 3g-3 = 3). Covariance: at colour i the coefficients are a_k^{(i)} = zeta^{-i(2k-1)} a_k, which is exactly how the true weights behave (wt_i = zeta^i wt_0 gives p_{2k-1}(1/wt_i) = zeta^{-i(2k-1)} p_{2k-1}(1/wt_0)). The frame is: rows R' = R phi_i^{-1}, i.e. the z-coefficients of the rows at colour i are multiplied by zeta^{-in} c_n with c_n = [z^n] exp(-a_1 z - a_2 z^3) (c_1 = -a_1, c_2 = a_1^2/2, c_3 = -a_1^3/6 - a_2; c_n has a-weight n, with weight 1 for a_1 and 3 for a_2); vertex CohFT phi_i . 1, the action of the scalar on the trivial rank-1 theory. By the group property of the R-matrix action, R . 1 = R' . (phi . 1), so the colour-summed total is independent of (a_1, a_2) and equals F_2; only the split into coloured graph terms moves.

**Vertex classes (shortcut with the exact identification).** The action of exp(Sum_k a_k z^{2k-1}) on the trivial theory is the Hodge-type class

    phi . 1 = exp( Sum_k (2k)!/B_{2k} a_k ch_{2k-1}(E) ) = exp(12 a_1 ch_1(E) - 720 a_2 ch_3(E) + ...).

Identification with the existing Hodge machinery (fq_hodge): a weight multiset w has Hodge_w = prod_a Lambda^vee(w_a)/w_a^g = exp(Sum_{k odd} s_k ch_k) with s_k = -(k-1)! Sum_a w_a^{-k}, and its QRR constant is h_w(z) = exp(-Sum_k B_{2k}/(2k(2k-1)) z^{2k-1} p_{2k-1}(1/w)), p_m = Sum_a w_a^{-m}. So a Hodge class is the point a_k = -B_{2k}/(2k(2k-1)) p_{2k-1} of the family, and then s_{2k-1} = -(2k-2)! p_{2k-1} = (2k)!/B_{2k} a_k: 12 for k = 1, -720 for k = 2, i.e. a_1 = -p_1/12, a_2 = p_3/360. Hence the vertex integrals are

    I_phi(g; d) = Int_{Mbar_{g,n}} (phi . 1) prod psi_i^{d_i} = Sum_{odd partitions K of D} hi(g, d, K) prod_k s_k^{m_k}/m_k!,   s_1 = 12 a_1, s_3 = -720 a_2,

homogeneous of a-weight exactly D = 3g-3+n-Sum d (PhiVertexTable: fq_hodge.VertexTable with these s_k, phases zeta^{-ik} s_k at colour i). Hand values, re-verified in code: I(1; 0) = a_1/2 (at a_1 = -1/(12w) this is -1/(24w) = Int(-lambda_1/w), fixing the sign), I(1; 1) = 1/24, I(2; empty) = a_1^3/10 + a_2/48, I(2; 4) = 1/1152. Points of the family: (0, 0) = point frame F5; (-1/12, 1/360) = Hodge_(1), frame F4 (I(2; empty) = 0 there, as Int_{Mbar_2} Lambda^vee(1) = 0); (-3/20, -7/2500) = the true class Hodge_{wt_0}, frame F0 (p_1(1/wt_0) = 9/5, p_3 = -126/125; I(2; empty) = -19/48000 = F_2(0)/5 there). The family contains all Hodge classes with arbitrary weights, so the agreement of PhiVertexTable with VertexTable at random two-weight multisets (p_1, p_3 then range over a Zariski-dense set) checks the vertex table as a polynomial in (a_1, a_2), not only at the three named points.

**Polynomial structure.** The a-weight of a vertex is D_v = 3g_v-3+n_v-Sum d; a row coefficient of z-order n has a-weight <= n (c_k P_{a,n-k}); so an edge coefficient E_{mn} has a-weight <= m+n+1 and a dilaton leaf Dil(d) <= d-1. Summing over a decoration, Sum_v(3g_v-3+n_v-Sum d) + Sum_e(m_e+n_e+1) + Sum_leaves(d_k-1) = 3 Sum g_v - 3|V| + 3|E| = 3g-3 = 3. Hence every coloured term is a polynomial of a-weight <= 3:

    T_{Gamma,c}(a) = T0 + T1 a_1 + T2 a_1^2 + T3 a_1^3 + T4 a_2,    T_k in the ring N,

and the five T_k are obtained by interpolation of the (unchanged) graph sum at the nodes (0,0), (1,0), (-1,0), (2,0), (0,1). Colour-shift covariance holds term by term, T_{Gamma,c+1}(a) = T_{Gamma,c}(a), so the 115 coloured graphs fall into 23 orbits.

**Valuation bookkeeping (conifold frame, s = 1/L, u = s^5).** From FORMAL_QUINTIC s.3: val P_{a,n} = -4n + (0,-1,3,2,1)_a, and every decoration of a genus-2 graph has valuation exactly -10 in the point frame. The a-weight-k part of a row coefficient of z-order n is c_k P_{a,n-k}, of valuation 4k higher; a vertex of psi-deficit D (a-weight D) lowers the half-edge and leaf z-orders by D in total, raising the valuation by 4D. So the a-weight-W part of every decoration sits at valuation exactly -10 + 4W: val T0 = -10, val T1 >= -6, val T2 >= -2, val T3, T4 >= 2, and the exponents of Cf(T_W) lie in the class -10 + 4W = -W mod 5: T0 in class 0 (a u-series), T1 in 4, T2 in 3, T3 and T4 in 2. Exact version (post-run check A): the Cf image of a numerator monomial with L-power m is s^{-m} times a u-series, so the exponent classes of Cf(T_k) are exactly {-m mod 5} over the L-powers of T_k and val Cf(T_k) >= -max L-degree(T_k); on the interpolated NEls the classes are {0}, {4}, {3}, {2}, {2} and the maximal L-degrees are 10, 6, 2, -2, -2, to all orders. Consequently the polar part of a coloured term is

    T0_{-10} s^{-10} + T0_{-5} s^{-5} + a_1 (T1_{-6} s^{-6} + T1_{-1} s^{-1}) + a_1^2 T2_{-2} s^{-2}:

a_2 never enters a polar term, and the a-independent part T0 (= the point-frame term) is untouched by any constant scalar. Because the four classes 0, 4, 3, 2 are distinct, Cf(T(a)) is a u-series if and only if a_1 T1 = 0, a_1^2 T2 = 0 and a_1^3 T3 + a_2 T4 = 0 separately (only T3 and T4 share a class).

## 2. Checks performed (logs/genus2_scalar.log, both primes, all PASS unless stated)

- fq_phi self-tests: the four hand values of I_phi; a-weight homogeneity for g <= 2, n <= 3; the two Hodge points; a_k^{(i)} = zeta^{-i(2k-1)} a_k from wt_i; PhiVertexTable = VertexTable at w = {1}, w = wt_0 and three random two-weight multisets, all colours, all (g, d) with g <= 2, n <= 3; (0, 0) = point class; phi_inv_coeffs = h_w^{-1} to z^3 at both Hodge points, covariance zeta^{-in} c_n, c_3 = -a_1^3/6 - a_2.
- Exactly in N, per coloured term: node (0,0) = F5; the interpolant at (-1/12, 1/360) = F4 and at (-3/20, -7/2500) = F0 (frames computed independently from the Hodge weights in genus2_frames.py); the interpolant at a random F_p point = the directly computed phi frame (excludes hidden monomials beyond the five); colour-shift equality of each T_k; Sum over coloured graphs of T0 = F_2 = total(F5) = total(F0); Sum T_k = 0 for k = 1..4.
- Cf of T0..T4 for all 115 coloured graphs (precision O(s^80)): val T0 = -10 for all 115; val T1 >= -6 (minimum -6; zero series for 20), T2 >= -2 (minimum -2; zero for 20), T3 >= 2 (zero for 80), T4 >= 2 (zero for 100); exponent classes as predicted. The exact all-orders version is check A of genus2_scalar_post.log.
- Table consistency: Galois back-substitution of the Q(zeta) coordinates reproduces every entry; totals row (1/1200, -7/6000, 0, 0, 0); the T0_{-10} colour sums reproduce the t_c^{-2} coefficients of FORMAL_QUINTIC s.3 graph by graph; at a_1 = -1/12 and -3/20 the vertex-graph entries reproduce the s^{-6}, s^{-2} coefficients of F4 and F0 in genus2_frames.log (37/51200 x (-1/12) = -37/614400, -119/7680 x 1/144 = -119/1105920; -111/1024000, -357/1024000).
- Two primes: valuations, tables, solves and the optional system identical. The log of record ends with two FAIL lines, both one artefact: fq_fp.ratrec at 2^61-1 reconstructs only denominators below 10737418, and the two entries 1807/15360000 (vertex graph) and -31/15360000 (graph 4, (0,0)) exceed it, so the extra Q(zeta)-coordinate column reads None there and the derived "reconstructed" cross-prime comparison of that column fails; the values themselves agree (1393615335851750081 and 1550751631410234689 are exactly those rationals mod 2^61-1). Check B of genus2_scalar_post.log redoes the cross-prime comparison by residues: all 115 table entries agree. The 2^61-1 column of the table is therefore residue-checked, not reconstructed, for those two entries.

## 3. The table: polar coefficients of the 23 orbit representatives

Columns: T0_-10, T0_-5 (a-independent), T1_-6, T1_-1 (coefficient of a_1), T2_-2 (coefficient of a_1^2); these are all the polar coefficients. Graphs as in FORMAL_QUINTIC s.3: 1 = genus-2 vertex, 2 = genus-1 vertex with a loop, 3 = genus-0 vertex with two loops, 4 = two genus-1 vertices joined by an edge, 5 = genus-0 vertex with a loop joined to a genus-1 vertex, 6 = theta graph, 7 = dumbbell. z = zeta, w := z^2 + z^3; the mixed colourings (0, j) and (0, 5-j) of graphs 4, 6, 7 have identical entries. Entries at 2^89-1 reconstructed from the Galois conjugates; identical by residues at 2^61-1.

| graph, colouring | T0_-10 | T0_-5 | T1_-6 | T1_-1 | T2_-2 |
|---|---|---|---|---|---|
| 1 (0) | 1807/15360000 | -7/40000 | 37/51200 | 0 | -119/7680 |
| 2 (0) | -811/1280000 | 763/960000 | -343/64000 | 19/1600 | 89/1920 |
| 3 (0) | -57/320000 | 89/100000 | 87/16000 | -2/125 | -7/160 |
| 4 (0,0) | -31/15360000 | 1/240000 | -73/256000 | 1/1500 | 71/7680 |
| 4 (0,1) = 4 (0,4) | 1/40000 - w/30000 | -1/30000 + 11w/240000 | -w/4000 | -1/3000 + w/12000 | 1/100 - w/200 |
| 4 (0,2) = 4 (0,3) | 7/120000 + w/30000 | -19/240000 - 11w/240000 | 1/4000 + w/4000 | -1/2400 - w/12000 | 3/200 + w/200 |
| 5 (0,0) | -17/1280000 | 313/4800000 | 67/64000 | -97/24000 | -53/1920 |
| 5 (0,1) | -1/40000 - z/20000 - z^2/20000 - z^3/40000 | 23/300000 + 7z/37500 + 17z^2/100000 + 17z^3/200000 | 3/2000 + z/1000 + z^2/2000 | -41/15000 - z/750 + 7z^2/10000 + 3z^3/2500 | -1/50 + w/100 |
| 5 (0,2) | 1/40000 + z/40000 + z^3/20000 | -7/75000 - 17z/200000 + z^2/60000 - 17z^3/100000 | 1/1000 - z/2000 + z^2/2000 - z^3/2000 | -103/30000 + z/2000 - 61z^2/30000 - 7z^3/10000 | -3/100 - w/100 |
| 5 (0,3) | -z/40000 + z^2/40000 - z^3/40000 | -1/120000 + 17z/200000 - 17z^2/200000 + 61z^3/600000 | 3/2000 + z/2000 + z^3/1000 | -59/15000 - z/2000 - 3z^2/2500 - 19z^3/7500 | -3/100 - w/100 |
| 5 (0,4) | 1/40000 + z/20000 + z^2/40000 | -11/100000 - 7z/37500 - 61z^2/600000 - z^3/60000 | 1/2000 - z/1000 - z^2/1000 - z^3/2000 | -7/5000 + z/750 + 19z^2/7500 + 61z^3/30000 | -1/50 + w/100 |
| 6 (0,0) | 1/96000 | -3/40000 | -1/1600 | 3/1000 | 1/80 |
| 6 (0,1) = 6 (0,4) | 1/12500 - w/9375 | -1/12500 + w/6250 | 0 | 0 | 0 |
| 6 (0,2) = 6 (0,3) | 7/37500 + w/9375 | -3/12500 - w/6250 | 0 | 0 | 0 |
| 7 (0,0) | 1/64000 | -9/80000 | -3/3200 | 9/2000 | 3/160 |
| 7 (0,1) = 7 (0,4) | 1/40000 - w/80000 | -3/25000 + 29w/400000 | -1/1000 + w/2000 | 3/1250 - 17w/10000 | 1/100 - w/200 |
| 7 (0,2) = 7 (0,3) | 3/80000 + w/80000 | -77/400000 - 29w/400000 | -3/2000 - w/2000 | 41/10000 + 17w/10000 | 3/200 + w/200 |

Per-graph colour sums (all 5^V colourings, rational): graph 1: 1807/3072000, -7/8000, 37/10240, 0, -119/1536; graph 2: -811/256000, 763/192000, -343/12800, 19/320, 89/384; graph 3: -57/64000, 89/20000, 87/3200, -2/25, -7/32; graph 4: 843/1024000, -53/48000, 11/10240, -1/240, 455/1536; graph 5: 3/51200, -67/192000, 71/2560, -373/4800, -245/384; graph 6: 87/32000, -143/40000, -1/320, 3/200, 1/16; graph 7: 9/12800, -59/16000, -19/640, 7/80, 11/32. Totals: 1/1200, -7/6000, 0, 0, 0 (the first two are the polar part of F_2, FORMAL_QUINTIC s.4; the zeros are Sum T_k = 0).

For the vertex graph (every colour) the a-dependence begins T1 = 37/51200 s^-6 + ..., T2 = -119/7680 s^-2 + ..., T3 = -131/3456 s^2 + ..., T4 = 23/1152 s^2 + ..., with T0 = 1807/15360000 s^-10 - 7/40000 s^-5 + 67/1440000 + .... Leading coefficients (T3_2, T4_2) of the other monochrome representatives: graph 2 (-1/288, -1/48), 3 (1/8, 0), 4 (-397/3456, 1/1152), 5 (23/96, 0), 6 (-1/12, 0), 7 (-1/8, 0); all mixed colourings (0, 0). For the 20 mixed colourings of the theta graph the entire a-dependence vanishes (T1 = T2 = T3 = T4 = 0 as series): their term is the point-frame term in every scalar frame.

## 4. The solve for a constant covariant scalar

Full system (all negative powers of all 110 coloured edge graphs vanish): inconsistent already at s^-10, whose coefficient T0_-10 is a-independent and non-zero for 110 of 110 coloured edge graphs (values in the table; the colour sums are the t_c^-2 coefficients of FORMAL_QUINTIC s.3). No constant choice of (a_1, a_2) removes the a-independent leading pole of any edge graph; a constant scalar can only add polar terms.

Weak system (the a-dependent polar terms only): a_1 T1_-6 = a_1^2 T2_-2 = a_1 T1_-1 = 0 for all coloured edge graphs. Each of T1_-6, T1_-1, T2_-2 is non-zero for 90 of the 110 coloured edge graphs (18 of the 22 edge representatives; the 20 zeros are the mixed theta colourings). Solution set: a_1 = 0, a_2 arbitrary.

u-series statement. Among the 115 coloured graphs the Cf series T1, T2, T3, T4 are non-zero for 95, 95, 35 (the monochrome graphs) and 15 (monochrome graphs 1, 2, 4). Since the exponent classes 4, 3, 2 are distinct from each other and from 0, every coloured term is a u-series if and only if a_1 = 0 (from T1 or T2) and then a_2 T4 = 0, i.e. a_2 = 0. The only constant covariant scalar frame in which every coloured term is a u-series is the point frame (a_1, a_2) = (0, 0); the Hodge frames F4 and F0 are not (their s^-6 and s^-2 terms are a_1 T1 and a_1^2 T2 evaluated at a_1 = -1/12 and -3/20).

## 5. Optional: a q-dependent leading system (beyond the user's constant-scalar choice)

Since T0 is frame-independent, the only way a covariant scalar could cancel the s^-10 pole of an edge graph is with q-dependent a_k. The decomposition T0 + T1 a_1 + ... + T4 a_2 remains valid for series a_k (the group property of the R-matrix action holds pointwise in q). The grading forces the ansatz a_1 in L^4 Q[[u]], a_2 in L^12 Q[[u]], i.e. a_1 = alpha_1 s^-4 (1 + O(u)), a_2 = alpha_2 s^-12 (1 + O(u)) with constants alpha_1, alpha_2 and valuations exactly L^4, L^12: with a_1 of lower L-degree the term a_1 T1 cannot reach s^-10 and T0_-10 stands alone; with higher L-degree the terms a_1 T1, a_1^2 T2, a_1^3 T3 land at distinct new pole orders that the two-vertex graphs cannot cancel simultaneously (this reasoning is the referee's; not computed). The leading s^-10 coefficient of every coloured term then gives

    T0_-10 + alpha_1 T1_-6 + alpha_1^2 T2_-2 + alpha_1^3 T3_2 + alpha_2 T4_2 = 0,

one equation per representative (23; cubic in alpha_1, linear in alpha_2; coefficients in s.3). Eliminating alpha_2 with the vertex-graph equation (T4_2 = 23/1152) leaves 22 polynomials in alpha_1 of degrees {0, 2, 3}: the four mixed theta colourings give non-zero constants (their equation is T0_-10 = 0 alone), the other mixed colourings quadratics (T3_2 = T4_2 = 0 there), the monochrome graphs cubics. The gcd over F_p is 1 at both primes: no common root, so no covariant scalar frame of this form, constant or q-dependent, kills the leading pole of all coloured edge graphs at once. The conclusion rests neither on the theta colourings nor on the pivot (post-run check C): the gcd is also 1 without the four theta mixed colourings (18 equations), for the six monochrome cubics alone, and for the twelve non-theta mixed quadratics alone. Only pairs of equations share roots (68 pairs; (6,(0,0)) and (7,(0,0)) are proportional cubics; the Galois pairs (4,(0,1)), (4,(0,4)) and (4,(0,2)), (4,(0,3)) coincide; the monochrome theta and dumbbell cubics share a quadratic factor with each mixed dumbbell colouring).

## 6. Caveats and referee items

- Genus 2 only; a_3 and higher never enter. Everything is exact finite data at two primes; no human has read this. Nothing here bears on genus >= 3, where scalars of higher z-order enter and the phase argument (total z-order < 5) no longer makes the colour-summed total trivially frame-independent.
- Referee pass (AI, fresh code, from the logged data only): the vertex shortcut checked by hand against the PPZ action at (g, n) = (1, 1) (Mumford's 12 lambda_1 = kappa_1 - psi + delta); the identification a_1 = -p_1/12, a_2 = p_3/360 and the Hodge points recomputed; I(2; empty) from Faber's lambda_1^3 = 1/2880, lambda_1 lambda_2 = 1/5760; the degree bound and the valuation shift re-derived; the leading coefficients (T0_-10, T1_-6, T2_-2, T3_2, T4_2) of graphs 2 and 4 in four colourings recomputed from the logged row data by independent code, all equal to the table; the elimination and the gcd statements redone exactly over Q(zeta); the two FAIL lines traced to the reconstruction bound. Verdict: sound with fixes; the fixes (residue comparison, exact exponent classes, sub-system gcds, the ansatz restriction, the distinct-classes remark) are incorporated above. No major item.
- Not independently re-derived: the vertex graph's entries (1807/15360000, 37/51200, -119/7680, -131/3456, 23/1152) beyond their consistency with the F4 and F0 values of genus2_frames.log and the point-frame value checked earlier (genus-2 Hodge integrals with psi insertions); the subleading columns T0_-5 and T1_-1 (these need the u-corrections of the periods); graphs 3, 5, 6, 7 (their entries involve Delta^-1 and genus-0 edge structure); the in-N identities and the interpolation, taken from the log (both primes PASS).
- The Q(zeta)-coordinate column assumes sigma_m(T(Gamma, c)) = T(Gamma, m c) (lambda_i -> lambda_{mi}); justified since every ingredient is rational in the lambda_i, and verified by exact back-substitution for all entries at both primes and by the referee's leading-order values for graph 4 at (0,1), (0,2).
- The log of record genus2_scalar.log ends with "SCALAR FAILURES" naming the two reconstruction checks; a trailing note points to genus2_scalar_post.log, where the same data pass by residues. genus2_scalar.py was left unchanged so that its log remains reproducible by its driver.
