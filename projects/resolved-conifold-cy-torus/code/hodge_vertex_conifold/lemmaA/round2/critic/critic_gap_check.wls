(* critic_gap_check.wls -- lemmaA round 2, completeness critic (2026-11-01). Self-contained; nothing here runs TR.
   Claim checked: R2's Corollary 7 (proved modulo ABDKS Prop. 3.1 and Eynard 2014, TR for the Hodge CohFT only),
       sigma^h Int_{Mbar_h} K Omega = -(2h-3)! [l^(2h)] S(e1 l) S(e2 l),   c = 1, e1 = -tau, e2 = 1 + tau, S(y) = (y/2)/Sinh[y/2],
   converted to the frame-B normalisation (c = 1 + x, tau = -x/c, Z[h,0] = Fpot_h(q_0) = c^-(3h-3) Int K Omega), is the vertex
   identity F4 at n = 0 in every genus:
       (-x c^3)^(h-1) Z[h,0] = gap_h := t^(2-2h) coefficient of P_h Li_{3-2h}(e^t)/(e1 e2) at e1 = x, e2 = 1,
   with Pg, polarOfRational, cyPolar copied verbatim from rc_lib.m (lines 241-245). Also checks that the subleading polar
   coefficients of P_h Li_{3-2h}(e^t)/(e1 e2) vanish (the single-pole form of (H)).
   Usage: wolframscript -file critic_gap_check.wls        Expected: every line "True ... True"; run time well under a minute. *)
Pg[g_] := Module[{tv}, SeriesCoefficient[Series[e1 e2 tv^2/(4 Sinh[e1 tv/2] Sinh[e2 tv/2]), {tv, 0, 2 g}], 2 g]];
polarOfRational[rq_, gg_] := Module[{ser = Series[rq /. Q -> Exp[t], {t, 0, -1}]}, Table[SeriesCoefficient[ser, -n], {n, 2 gg - 2, 1, -1}]];
cyPolar[gg_] := polarOfRational[Pg[gg] PolyLog[3 - 2 gg, Q]/(e1 e2), gg];
SS[y_] := (y/2)/Sinh[y/2];
cor7[h_] := Module[{tau = -x/(1 + x), sig},
   sig = tau (1 + tau);
   -(2 h - 3)! SeriesCoefficient[SS[-tau l] SS[(1 + tau) l], {l, 0, 2 h}]/sig^h];
ok = Table[Module[{c = 1 + x, lhs, pol, gap, rest, t1, t2},
    lhs = (-x c^3)^(h - 1) c^(-(3 h - 3)) cor7[h];
    pol = cyPolar[h] /. {e1 -> x, e2 -> 1};
    gap = First[pol]; rest = Rest[pol];
    t1 = Together[lhs - gap] === 0;
    t2 = AllTrue[Together /@ rest, # === 0 &];
    Print["h = ", h, ": (-x c^3)^(h-1) Z[h,0] == gap_h: ", t1, ";  subleading polar coefficients vanish: ", t2];
    t1 && t2], {h, 2, 10}];
Print["ALL TRUE: ", And @@ ok];
