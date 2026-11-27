# Local P2: preliminary scoping of edge-graph regularity (workflow wf_e39f8d50-a99, 2026-11-01)

Scratch scripts: session scratchpad lp2_explore (frozen.pl, leafrow.pl, frob.pl). Full per-lens reports: workflow journal.

PRELIMINARY ASSESSMENT: edge-graph regularity for K_P2 x C^2 (fibrewise torus T_F, frame B, mode-B triple Hodge vertices, conifold polarisation)

Cross-checks I did myself:
- lp2_lib.m lines 16-32 (LLt0, sqDsym). By hand from sqDsym, Delta_0 = 3 eps1 eps2 E^3 t^2/(1+3t)^2 with t = y^(1/3). Only the real branch t -> -1/3 is singular.
- I reran frozen.pl on cluster_genus4_hodge/cPP_z9.m at eps = (2,1). The frame-B row (h_w^-1) is E/(E+z) through z^9, up to float error 1e-5. The h^+1 twist gives a non-geometric row.
- round1_status.md:319 and conifold_regularity_note.md:175 (third-order Picard-Fuchs, 609 odd-order decorated graphs) say what the reports quote.

(a) INGREDIENTS
- Vertex class (same triple Hodge class), leaf scale c = E = eps1+eps2: transfers verbatim. ESTABLISHED: memory refgw-hodge-vertex.md; frozen.pl.
- Lemma A = M_irr(h), Conjecture M, string/dilaton reduction, Theorem C: transfer verbatim. These are Hodge-side statements: lemmaA/round2/round2_status.md, note s.3a. They apply once the local-P2 leaves are 1/(c - j psi). That is PLAUSIBLE-strong, not yet shown beyond the identity row.
- Graph combinatorics (M_irr only needed for h <= g-1): transfers verbatim. ESTABLISHED.
- F1, Delta_0 = 3 eps1 eps2 E^3/t_c^2 times a unit, two regular colours: transfers with changes (rank 3). ESTABLISHED by a closed form from lp2_lib.m. The ratio -3 to the conifold matches the stored (-3)^(g-1).
- F2, frozen row E/(E+z): transfers verbatim. ESTABLISHED to z^9 (exact BigRat in the report, float rerun by me), and at t_c = 0 to all orders by BS Thm. cfasym. It does not depend on the polarisation.
- L1, dilaton/identity-row leaf support: [u^m] lies in span{1/(E+jz), j <= m+1}. ESTABLISHED to u^5, z^9 (leafrow.pl, exact). All orders: UNKNOWN.
- H-row (P1 = C1i(L P0 + ...)): with changes. PLAUSIBLE, since it is a derivative of the identity row up to the scalar C1i.
- H^2-row and edges (F3, L2, L3, the E_00 loops, a rank-3 analogue of the E_11 flow): new work. UNKNOWN. This is where the polarisation XX enters.
- Thimble mechanism: with changes. The Gauss 2F1 is replaced by a third-order, 3F2-type equation with exponents {0, 1, 1+c/z} at the conifold. A Frobenius check (frob.pl, numerical, sign assumed) gives simple poles at z = -c/j only: PLAUSIBLE. A closed form like the KEY FORMULA: UNKNOWN.
- Theorem R (counting): with changes, a formal transfer once LC_P2 holds. PLAUSIBLE.
- Frame-B formula (Prop. 13.2 analogue): with changes. Steps (a), (c), (e) carry over. Steps (b) and (d) have no Beta-integral Lemma 13.3, and on T_F the h_nw cancellation is lost; they rely on BS's normalisation (five-weight QRR constants, delta-curve limit). Checked exactly for g = 2-4 (hodge_vertex/README.txt). No written proof.
- Data: ESTABLISHED. Edge graphs regular at genus 2-4 (6/7, 41/42, 378/379 regular; the exception each time is the single-vertex graph). At genus 4 the vertex graph carries exactly the gap. The pole bound Sum_v (2g_v-2)_+ holds.

(b) MAIN DIFFERENCES AND RISKS
Differences:
1. Rank 3. There is one singular idempotent and two regular ones, which are complex conjugate Galois branches of t = y^(1/3). Pitfall: the singular colour is colour 0 only on the real branch.
2. Polarisation. Graph terms lie in Q(eps)(y)[XX], and regularity holds only on shell. XX(0) = -2/9 is forced by the indicial equation, so the conifold polarisation is canonical. Off-shell (1+27y)^-1 poles occur in 155 of the regular genus-4 graphs. The conifold has none of this.
3. No t -> -t symmetry. The 609 odd-order decorated graphs are harmless for Lemma A, which covers all D, but subleading ratios are not universal.
4. Gap input. Only the leading order is known TR-free in all genera (BS Thm. cfasym); the full refined gap is a conjecture in BS (Rmk. cfgap).

Risks:
1. Edges and legs may bring in leaves of a non-weight scale or polynomial remainders. Single generic-scale leaves do not vanish (note s.7), so this would break the argument.
2. Per-factor F3 is untested. The evidence is graph-level and summed over colourings.
3. The sign or orientation of the leaf scale: scale -c is the mode-A failure. The F2 check supports +c on the identity row only.
4. The rigour of the frame-B formula on T_F rests on BS's normalisation.

(c) MINIMAL PLAUSIBLY PROVABLE STATEMENT
Conditional theorem. Assume LC_P2 (leaf structure at the singular idempotent through order t_c^(2h-3), in the conifold polarisation) and F3_P2. Then every genus-g stable graph with an edge in the frame-B Hodge-vertex sum for local P2 is regular at t_c = 0. The single-vertex graph's polar part, as a Hodge-side quantity, equals (-3)^(g-1) gap_g^conifold t_c^(2-2g). This follows from Theorem C, Theorem R, Lemma A and R2 Cor. 7.

Unconditionally, within reach: F1, F2 and L1 at low order, already essentially in hand. Caution: combined with a frame-B formula for local P2, this would prove BS's conjectural refined gap. It should be scoped and refereed carefully.

(d) NEXT STEPS, RANKED
1. Perl, cheap: extend leafrow.pl to the H-row P1 (free of c1 up to C1i). Confirm the shifted leaf support and the 1/u pole. Also rerun frob.pl with lambda -> -lambda against the cPP denominators, to fix the sign once.
2. Mathematica, seconds, when a kernel is free: expand Dil, Leaf and Edge at colour pairs (0,0), (0,1), (0,2), with XX -> XXcf from lp2_tc_section.m, to O(t_c^3). Test F3 (rank-one limits) and the L2/L3 partial fractions in 1/(E - j psi), copying hvc_leaf_proof_check*.wls.
3. Check per-colouring regularity at the two regular idempotents from lp2_lib.m.
4. Port hvc_vertex_conjecture_g2.wls: in the exact singular genus-2 vertex with legs and loops, the t_c^-2 and t_c^-1 coefficients must vanish. Then try one genus-3 case.
5. Write down the T_F equivariant QDE from the notebook (grep only) and confirm the 3F2 form and the exponent 1 + E/z. Then attempt an L1 proof via Frobenius on the order-3 equation.
6. Sketch Prop. 13.2-P2, with the polarisation paragraph: algebraic independence of XX and finite generation, BS prop:algindABC and prop:fingenTLE.

CAVEATS:
- The F2 rerun by me was a float check at one weight choice (eps=(2,1)); the exact BigRat checks and leafrow.pl L1 checks are from the scoping report, not rerun.
- frob.pl (3F2-type exponents, simple poles) is numerical and assumes a sign for the O(-3) weight; not verified against project notebooks.
- BS index slip (beta_{2g} vs beta_{2g-2}) and BS line citations were not re-read by me.
- No large-radius GV check of genus-4 local-P2 F40 found on file.
- No Mathematica run; all edge/leg (F3, L2, L3) statements remain untested factor by factor.
