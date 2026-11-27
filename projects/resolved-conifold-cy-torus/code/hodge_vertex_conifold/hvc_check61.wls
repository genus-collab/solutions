(* hvc_check61.wls -- compare the conjectural closed one-leaf formula (6.1) of proof_M_residue.md with the exact
   one-leaf functional Xi_h(j) (symbolic in j and x) for h = 2, 3, 4.
   (6.1), c = 1: sigma^h Xi_h(j) = prod_{l=1}^{2h-2}(j-l) [lam^(2h)] S(e1 lam)^(1+j e2) S(e2 lam)^(1+j e1),
   S(y) = (y/2)/sinh(y/2), e1 = -tau, e2 = 1+tau, tau = -x/(1+x), sigma = tau(1+tau); Xi^(c) = c^-(3h-1) Xi^(c=1) for r = 1. *)
SetDirectory[DirectoryName[$InputFileName]];
Get["hvc_mp_header.m"];
SS[y_] := (y/2)/Sinh[y/2];
form61[h_, j_] := Module[{tau = -x/(1 + x), e1, e2, sig, lam, ser},
   e1 = -tau; e2 = 1 + tau; sig = tau (1 + tau);
   ser = Series[Exp[(1 + j e2) Log[SS[e1 lam]] + (1 + j e1) Log[SS[e2 lam]]], {lam, 0, 2 h}];
   Product[j - l, {l, 1, 2 h - 2}] SeriesCoefficient[ser, 2 h]/sig^h/(1 + x)^(3 h - 1)];
Do[Module[{t1 = AbsoluteTime[], ex = Xi[h, 1] /. jj[1] -> jv, fm, d},
   fm = form61[h, jv]; d = Together[ex - fm];
   Print["genus ", h, ": Xi_h(j) - (6.1) = ", d, "   (", Round[AbsoluteTime[] - t1, 0.1], " s)"]], {h, 2, 4}];
