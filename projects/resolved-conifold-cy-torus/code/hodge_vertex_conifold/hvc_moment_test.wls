(* hvc_moment_test.wls -- prediction of the shear / escaping-branch-point mechanism: the first non-vanishing y-moment
   Mom_g := Res_{v=tau} y^(2g-1) omega_{g,1} of the base curve C_tau is independent of tau.
   Via P2 and Conjecture M: Mom_g = (2g-1) [x^(2g-1)]F_{g,1} = (2g-1) Xi_g(2g-1) binom((2g-1)(tau+1)-1, 2g-2)/sigma^(g-1),
   Xi_g(j) from the closed one-leaf formula (6.1) (exact for g <= 4, TR-checked g = 5, 6 at tau = 1, -1/2). Also the next moment. *)
SS[y_] := (y/2)/Sinh[y/2];
xi61[h_, j_, tau_] := Module[{e1 = -tau, e2 = 1 + tau, sig = tau (1 + tau), lam},
   Product[j - l, {l, 1, 2 h - 2}] SeriesCoefficient[Series[Exp[(1 + j e2) Log[SS[e1 lam]] + (1 + j e1) Log[SS[e2 lam]]], {lam, 0, 2 h}], 2 h]/sig^h];
Fmu[g_, mu_, tau_] := xi61[g, mu, tau] Binomial[mu (tau + 1) - 1, mu - 1]/(tau (1 + tau))^(g - 1);
mom[g_, N_, tau_] := (-1)^(N + 1) N Sum[Fmu[g, mu, tau] Binomial[-(tau + 1) mu, N - mu], {mu, 1, N}];
Do[Print["g = ", g, ":  Mom_(2g-1) = ", Factor[FunctionExpand[mom[g, 2 g - 1, T]]], "    Mom_(2g) = ", Factor[FunctionExpand[mom[g, 2 g, T]]]], {g, 2, 5}];
