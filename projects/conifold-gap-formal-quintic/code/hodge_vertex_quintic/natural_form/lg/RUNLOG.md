# Run log, natural_form/lg (2026-11-10)

All runs: single process, `py -I <script>` from this folder, pure Python, exact arithmetic. No file outside this folder was written; the background production runs (genus3_natural.py) and their outputs were not touched. Times are wall-clock on the shared machine.

| # | command | output | result |
|---|---|---|---|
| 1 | `py -I lg_dump_symbol.py 89` and `... 61` | lg_symbol_89.json, lg_symbol_61.json, lg_dump_*.console | fq3_par.setup((p, "0", "0", 34, 6)) in 1.5 s per prime; checks val_ok, hat_val_ok, symbol_check True, 0 symplectic failures, covariance a+1-n True |
| 2 | `py -I lg_sp.py` | console | Stirling (8 coefficients), Airy (8), exact-Gaussian and reparametrisation checks pass; J_k, k = -2..5, through eps^7 printed |
| 3 | `py -I lg_match.py` (first version) | lg_match_run1.log | rows identified: a = 0..4 <-> k = 0, -1, 3 (= -2), 2, 1, c = -4, with N = (1, -1, -1, 1, 1). The script's EXPECTATION N_a = 1 failed for a = 1, 2 (reported FAIL), and so did the 175-residue block (70 mismatches = the rows a = 1, 2); the control (P,Q) = (1,4) unexpectedly matched at k = -2, c = -1 (the same model under x -> 1/x). Colour matrices (D) passed (signs cancel). |
| 4 | `py -I lg_match.py` (expectations corrected: N = (1,-1,-1,1,1); (1,4) recognised as the inversion; negative test of the guess H^a <-> x^a dx added) | lg_match.log | ALL PASS at both primes: (A) agnostic fit, (B) 175 residues, guess row a = J_a only for a = 0, unit only k = 0, (C) controls, (D) 175 colour-matrix entries |
| 5 | `py -I lg_gm.py` | lg_gm.log, lg_symbol_rows_Q.json | ALL PASS: (R1), (R2), (R3) on the stationary-phase series k = -4..8 through zhat^10; P_4(m) = -4m; (MC) recursion = engine for 5 rows through zhat^10; symplecticity of the recursion output through zhat^40; 175 residues per prime; C_1[0][0] = 1/20, dF_1 = -(1/12) du/u |
| 6 | `py -I lg_wkb_check.py 10`, then `... 40` | lg_wkb_check.log (Nz = 40) | ALL PASS: deg_L phi_n = 4n and [L^{4n}] phi_n = zeta^{-in} g^0_n for n <= 40, colours 0 and 1, both primes (0.03 s per WKB; the Laurent polynomials are sparse). Inline control: g^4 in place of g^0 is rejected; f_7 reconstructs to -7749624763881/512000000000 = the recursion value. |
| 7 | inline probe Nz = 8, then `py -I lg_rows_highorder.py 12`, `... 18`, `... 8 1/3 -2/7` | lg_rows_highorder.log | ALL PASS: 325 residues (Nz = 12, 3-4 s per prime), 475 residues (Nz = 18, 7-8 s per prime), and 225 residues in the conifold frame (alpha, beta) = (1/3, -2/7); symplectic failures 0, covariance True in every case |
| 8 | `py -I lg_frob.py` | lg_frob.log | ALL PASS: flat metric at 4 random rational points (and the control with the naive coordinate a_4), c_555 = 1/a_5, unit d/dt4, cyclic multiplication by d_5 on the line, degrees/spectrum bookkeeping, canonical-frame equation (5/4)[C_{n+1}, Z] = (V - n) C_n for the dump (n <= 5) and the recursion (n <= 30) at both primes |
| 9 | inline: leading terms of the conifold-frame images | lg_cf_leading.log | Cf(C0) = 1 + O(u), Cf(C1) = -1 + O(u) for (alpha, beta) = (0, 0) and (1/3, -2/7) |
| 10 | web: arXiv:2211.12259 (abstract, introduction, s.1 read), arXiv:1502.01774, hep-th/9505122, arXiv:1401.5725, math/0607012 (abstracts / search summaries) | LG_MODEL.md s.5.5 | the manifold is H_{0;1,N-1} with N = 5 of arXiv:2211.12259: same superpotential, flat coordinates, eta, d = 1/2, mu, special point p^4 + 1/p |

Negative / corrected items, for the record:
- the conjectured labelling H^a <-> x^a dx is wrong for a != 0; the data give x^{-a} dx (mod 5), i.e. k = c_a = (0, -1, 3, 2, 1);
- rows a = 1, 2 carry the sign -1 (run 3), explained by H = -d/du in LG_MODEL.md s.4;
- the potentials x^3 + u/x^2 and x^2 + u/x^3 do not reproduce the identity row for any k in -6..6; x + u/x^4 does, being the same model;
- the identity row is J_k only for k = 0 (not dx/x, not x dx);
- the rational reductions of 2D Toda (logarithmic primitive form, non-conformal) are a relative of, not the same as, the structure found (KP type, primitive form dx).
