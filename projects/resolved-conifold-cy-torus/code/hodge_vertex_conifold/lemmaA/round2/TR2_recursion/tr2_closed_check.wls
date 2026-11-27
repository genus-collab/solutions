(* tr2_closed_check.wls -- requested independent check (NOT run by route TR2: no wolframscript for subagents).
   (a) Closed one-point formula of C_tau from ABDKS Prop. 3.1 (n = 1), symbolic in tau and j, for g <= gmax:
         [x^j] F_{g,1} = (1/j^2) Sum_{p=0}^{g} A_{g-p} B_p Binomial[j(tau+1) - 2p, j + 1 - 2g],
         A_m = [t^(2m)] Exp[Sum_k t^(2k) (2k-1)! j tau c_k(j, 1/tau)],  B_p = [t^(2p)] Exp[Sum_k t^(2k) (2k-1)! (-j(tau+1)) c_k(j, -1/(tau+1))],
         c_k(u, a) = [t^(2k)] S(u t)/S(a t),  S(t) = 2 Sinh[t/2]/t.
   (b) Norlund forms: A_m = Pochhammer-type (-N)_(2m) [y^(2m)] S1(y)^(1-N) tau^(-2m), B_p = (M)_(2p) [y^(2p)] S1(y)^(1+M) (tau+1)^(-2p),
       N = j tau, M = j(tau+1), S1 = 1/S, (a)_n = a(a-1)...(a-n+1)  -- as polynomial identities in j and tau.
   (c) (6.1): sigma^g Xi_g(j) = Product[j - l, {l, 1, 2g-2}] [lam^(2g)] S1(tau lam)^(1+M) S1((1+tau) lam)^(1-N), with
       Xi_g(j) = sigma^(g-1) [x^j]F_{g,1} / Binomial[j(tau+1)-1, j-1]  (P2 normalisation)  -- as an identity of functions of j (use
       FunctionExpand / Pochhammer form of the binomials) and tau.
   (d) Optional: compare (a) at tau = 2, 3 with hvc_tr_tau.wls output (Xi^TR_g(mu), g <= 4).
   Usage: wolframscript -file tr2_closed_check.wls <gmax>      Expected: all differences 0. *)
gmax = ToExpression[$ScriptCommandLine[[2]]];
S[t_] := 2 Sinh[t/2]/t; S1[t_] := 1/S[t];
ck[u_, a_, k_] := SeriesCoefficient[S[u t]/S[a t], {t, 0, 2 k}];
falling[a_, n_] := Product[a - i, {i, 0, n - 1}];
Aser[j_, tau_, g_] := Series[Exp[Sum[t^(2 k) (2 k - 1)! j tau ck[j, 1/tau, k], {k, 1, g}]], {t, 0, 2 g}];
Bser[j_, tau_, g_] := Series[Exp[Sum[t^(2 k) (2 k - 1)! (-j (tau + 1)) ck[j, -1/(tau + 1), k], {k, 1, g}]], {t, 0, 2 g}];
Do[
  Module[{A, B, An, Bn, N1 = j tau, M1 = j (tau + 1), diffA, diffB, closed, rhs61, xi, d61},
   A = Aser[j, tau, g]; B = Bser[j, tau, g];
   An = Table[falling[-N1, 2 m] SeriesCoefficient[S1[y]^(1 - N1), {y, 0, 2 m}] tau^(-2 m), {m, 0, g}];
   Bn = Table[falling[M1, 2 m] SeriesCoefficient[S1[y]^(1 + M1), {y, 0, 2 m}] (tau + 1)^(-2 m), {m, 0, g}];
   diffA = Simplify[Table[SeriesCoefficient[A, {t, 0, 2 m}] - An[[m + 1]], {m, 0, g}]];
   diffB = Simplify[Table[SeriesCoefficient[B, {t, 0, 2 m}] - Bn[[m + 1]], {m, 0, g}]];
   (* (6.1) termwise: compare sigma^(2g-1) [x^j]F / Binomial[M-1,j-1] with Product (j-l) [lam^(2g)] S1(tau lam)^(1+M) S1((1+tau)lam)^(1-N) *)
   closed = Sum[An[[g - p + 1]] Bn[[p + 1]] Pochhammer[M1 - 2 p - (j + 1 - 2 g) + 1, j + 1 - 2 g]/Gamma[j + 2 - 2 g], {p, 0, g}]/j^2;
   xi = (tau (1 + tau))^(2 g - 1) closed Gamma[j] Gamma[M1 - j + 1]/Gamma[M1];
   rhs61 = Pochhammer[j - 2 g + 2, 2 g - 2] SeriesCoefficient[S1[tau lam]^(1 + M1) S1[(1 + tau) lam]^(1 - N1), {lam, 0, 2 g}];
   d61 = FullSimplify[FunctionExpand[xi - rhs61], Assumptions -> j > 2 g && tau > 0];
   Print["g=", g, ": Norlund A diff ", diffA, ", B diff ", diffB, ", (6.1) diff ", d61]],
  {g, 1, gmax}];
