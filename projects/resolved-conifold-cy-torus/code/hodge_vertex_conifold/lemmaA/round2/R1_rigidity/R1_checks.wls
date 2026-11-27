(* R1_checks.wls -- lemmaA round 2, route R1: independent Mathematica checks (exact; no TR of any kind).
   (C1) residue identity  Res_mu u_t(mu)^n dmu = (-1)^(n-1) t^(n-1) Exp[-n tau^2 t] Pi_n 2F1(-n tau,-n tau;n;1-Exp[-t])   [symbolic tau]
        u_t(mu)^n = u(mu)^n Exp[-n K_t(mu)],  u(mu)^n = mu^-n (1-1/mu)^(n tau),  K_t(mu) = Sum_k kappa_k t^k ((1+tau) mu^k - tau (mu-1)^k)
   (C2) branch-point identity  d/dlambda [u^-j]_+ (lambda = 1+tau) = j Pi_j   [symbolic tau]
   (C3) curve identity (FC) on the leaf curve: phi_t(lambda) = W lambda - Sum_j P_j W^j ([u^-j]_+(lambda) - [u^-j]_+(0)) satisfies
        K_t(phi_t(lambda)) - tau Log[(phi_t-1)/(lambda-1)] + (1+tau) Log[phi_t/lambda] = Log[W Exp[tau^2 t]/H_1]  (independent of lambda)
   (C4) T1 from an INDEPENDENT implementation of the genus-G linearised KP system along the leaf curve (z-expansion of Theorem D,
        genus-0 curve Y_P(z) = z - sigma Sum_n Pi_n W^(n+1) z^-n/(n+1), W = 1 + Sum_j Pi_j P_j W^j; string V_n = n U_n + Sum_j j P_j G_jn)
   (C5) forcing theorem: rank test on its whole hypothesis range, G <= 24
   (C6) zero count Z1(D) + Z2(G,D) >= G - 1 for 8 <= G <= 60, odd 3 <= D <= 2G-3 (the Perl script count_zeros2.pl does G <= 1000)
   Usage: wolframscript -file R1_checks.wls      (expected: every line ends in True / 0 failures) *)

gbin[x_, k_Integer] := If[k < 0, 0, Product[(x - i)/(i + 1), {i, 0, k - 1}]];
Pin[j_Integer, tau_] := If[j <= 1, 1, Product[j tau + a, {a, 1, j - 1}]/(j - 1)!];
kap[k_Integer, tau_] := If[k == 1, tau + 1/2, BernoulliB[k]/(k k!)];
H2F1[j_, tau_, TT_, M_] := Sum[Pochhammer[-j tau, k]^2/(Pochhammer[j, k] k!) TT^k, {k, 0, M}];

(* ---------- (C1) ---------- *)
c1[n_, M_] := Module[{Kt, e, lhs, rhs, TT},
  Kt = Sum[kap[k, tau] t^k ((1 + tau) mu^k - tau (mu - 1)^k), {k, 1, M}];
  e = Normal[Series[Exp[-n Kt], {t, 0, M}]];
  lhs = Sum[t^m Module[{pl = CoefficientList[Expand[Coefficient[e, t, m]], mu]},
        Sum[pl[[i + 1]] (-1)^(i - n + 1) gbin[n tau, i - n + 1], {i, 0, Length[pl] - 1}]], {m, 0, M}];
  TT = 1 - Exp[-t];
  rhs = Normal[Series[(-1)^(n - 1) t^(n - 1) Exp[-n tau^2 t] Pin[n, tau] H2F1[n, tau, TT, M], {t, 0, M}]];
  Simplify[Expand[lhs - rhs]] === 0];
Print["(C1) residue identity, n = 1..5, to t^7: ", And @@ Table[c1[n, 7], {n, 1, 5}]];

(* ---------- (C2) ---------- *)
c2[j_] := Module[{pp}, pp = Sum[gbin[j tau + k - 1, k] lam^(j - k), {k, 0, j}];
  Simplify[(D[pp, lam] /. lam -> 1 + tau) - j Pin[j, tau]] === 0];
Print["(C2) branch-point identity, j = 1..14: ", And @@ Table[c2[j], {j, 1, 14}]];

(* ---------- (C3) ---------- *)
c3[tauv_, M_] := Module[{tau = tauv, TT, Hj, P, W, phi, q1, q2, Dl, Kp, diff, target},
  TT = 1 - Exp[-t];
  Hj[j_] := H2F1[j, tau, TT, M + 1];
  P[j_] := P[j] = Normal[Series[(-1)^j Pin[j, tau]/j t^(j - 1) Hj[j]/Hj[1]^j, {t, 0, M}]];
  W = 1; Do[W = Normal[Series[1 + Sum[Pin[j, tau] P[j] W^j, {j, 2, M + 1}], {t, 0, M}]], {M + 2}];
  phi = Expand[Normal[Series[W lam - Sum[P[j] W^j Sum[gbin[j tau + k - 1, k] lam^(j - k), {k, 0, j - 1}], {j, 2, M + 1}], {t, 0, M}]]];
  If[Simplify[Normal[Series[(phi /. lam -> 1) - 1, {t, 0, M}]]] =!= 0, Return["phi(1) != 1"]];
  q1 = Cancel[phi/lam]; q2 = Cancel[(phi - 1)/(lam - 1)];
  Dl = Normal[Series[tau Log[q2] - (1 + tau) Log[q1], {t, 0, M}]];
  Kp = Normal[Series[Sum[kap[k, tau] t^k ((1 + tau) phi^k - tau (phi - 1)^k), {k, 1, M}], {t, 0, M}]];
  diff = Expand[Kp - Dl];
  target = Normal[Series[Log[W] + tau^2 t - Log[Hj[1]], {t, 0, M}]];
  {Simplify[D[diff, lam]] === 0, Simplify[diff - target] === 0}];
Print["(C3) curve identity, symbolic tau, to t^4 (lambda-independent, constant = log C_a): ", c3[tau, 4]];
Print["(C3) curve identity, tau = 2/7, to t^7: ", c3[2/7, 7]];
Print["(C3) curve identity, tau = -3/5, to t^7: ", c3[-3/5, 7]];

(* ---------- (C4) ---------- *)
c4[tauv_, K_, seed_] := Module[{tau = tauv, sig, NM, TT, Hj, P, W, tr, Emat, S, y, U, V, G, Gab, Phi, pred, s},
  sig = tau (1 + tau); NM = K + 2; TT = 1 - Exp[-t];
  Hj[j_] := H2F1[j, tau, TT, K + 1];
  P[j_] := P[j] = Normal[Series[(-1)^j Pin[j, tau]/j t^(j - 1) Hj[j]/Hj[1]^j, {t, 0, K}]];
  W = 1; Do[W = Normal[Series[1 + Sum[Pin[j, tau] P[j] W^j, {j, 2, NM}], {t, 0, K}]], {K + 2}];
  tr[e_] := Normal[Series[e, {t, 0, K}, {a, 0, NM}, {b, 0, NM}]];
  (* Y(w) - Y(z) = (a - b)[1/(ab) + Sum_n c_n h_{n-1}(a,b)],  a = 1/z, b = 1/w,  c_n = sigma Pi_n W^(n+1)/(n+1) *)
  Emat = Expand[a b Sum[sig Pin[n, tau] W^(n + 1)/(n + 1) Sum[a^i b^(n - 1 - i), {i, 0, n - 1}], {n, 1, NM}]];
  S = 1; Do[S = tr[1 - Emat S], {NM + 1}];   (* S = 1/(1 + Emat), built by truncated fixed-point iteration *)
  SeedRandom[seed]; y = Table[RandomInteger[{-9, 9}], {L, 1, 2 NM + 2}];
  U = Table[If[n >= 2, y[[n - 1]], 0], {n, 1, NM}]; G = ConstantArray[0, {NM, NM}];
  Do[
    V = Table[Expand[n U[[n]] + Sum[j P[j] G[[j, n]], {j, 2, NM}]], {n, 1, NM}];
    Gab = tr[a b Sum[V[[n]] Sum[a^i b^(n - 1 - i), {i, 0, n - 1}], {n, 1, NM}] S];
    G = Table[Coefficient[Coefficient[Gab, a, i], b, n], {i, 1, NM}, {n, 1, NM}];
    U = Table[Expand[If[n >= 2, y[[n - 1]], 0] + Integrate[(Sum[D[P[i], t] G[[i, n]], {i, 2, NM}]) /. t -> s, {s, 0, t}]], {n, 1, NM}];
    U = Map[Normal[Series[#, {t, 0, K}]] &, U],
    {K + 2}];
  Phi = Normal[Series[Integrate[(Sum[D[P[i], t] U[[i]], {i, 2, NM}]) /. t -> s, {s, 0, t}], {t, 0, K}]];
  pred = Table[kap[k, tau] Sum[y[[L]] (-1)^(k - L) ((1 + tau) gbin[(L + 1) tau, k - L] - tau gbin[(L + 1) (1 + tau) + k - L - 1, k - L]), {L, 1, k}], {k, 1, K}];
  Table[Simplify[Coefficient[Phi, t, k] - pred[[k]]] === 0, {k, 1, K}]];
Print["(C4) T1 from the z-expansion KP system, tau = 2/7, orders 1..7: ", c4[2/7, 7, 1]];
Print["(C4) T1 from the z-expansion KP system, tau = -3/5, orders 1..6: ", c4[-3/5, 6, 2]];
Print["(C4) T1 from the z-expansion KP system, tau = 3, orders 1..6: ", c4[3, 6, 3]];

(* ---------- (C5) ---------- *)
forced[G_, dd_, q_, k0_] := Module[{cols, rows, r1, r2},
  cols = Range[dd + 1, q - 1];
  rows = Table[Table[If[NN < n, 0, (q - k0) gbin[k0 + NN - n - 1, NN - n] + k0 gbin[k0 - 1 - n, NN - n]], {n, cols}], {NN, 1, 2 G - 3, 2}];
  r1 = MatrixRank[rows]; r2 = If[Length[cols] > 1, MatrixRank[rows[[All, 2 ;;]]], 0];
  r1 - r2 == 1];
c5[Gmax_] := Module[{n = 0, bad = 0},
  Do[If[GCD[k0, q] == 1 && 3 (q - k0) < 2 G, n++; If[! forced[G, dd, q, k0], bad++]],
    {G, 3, Gmax}, {dd, 3, 2 G - 3, 2}, {q, dd + 2, 2 G - 2}, {k0, 1, Min[dd + 1, Floor[(q - 1)/2]]}];
  {n, bad}];
Print["(C5) forcing theorem, G <= 24: {tuples, failures} = ", c5[24]];

(* ---------- (C6) ---------- *)
Z1[dd_] := Floor[dd/2] + Sum[(EulerPhi[a]/2) Max[0, Floor[dd/(a - 1)] - If[Mod[dd + 1, a] == 0, 1, 0]], {a, 3, dd + 1}];
Z2[G_, dd_] := Sum[If[GCD[k0, q] == 1 && 3 (q - k0) < 2 G, 1, 0], {k0, 1, dd + 1}, {q, Max[dd + 2, 2 k0 + 1], k0 + Ceiling[2 G/3]}];
c6[Gmin_, Gmax_] := Module[{bad = 0}, Do[If[Z1[dd] + Z2[G, dd] < G - 1, bad++], {G, Gmin, Gmax}, {dd, 3, 2 G - 3, 2}]; bad];
Print["(C6) zero count, 8 <= G <= 60: failures = ", c6[8, 60], "   (G = 4..7 at D = 3 expected to fail: ", c6[4, 7], ")"];
