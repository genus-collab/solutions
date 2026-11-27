(* hvc_pf_gauss.wls -- PF equation of the conifold thimble versus Gauss: does (z d/dt)^2 - (c/tau) z d/dt + x/tau annihilate
   I = (1 - tau)^(-x/z) tau^(c/z + 1) F(tau) when F solves Gauss with alpha = beta = 1/z + 1, gamma = c/z + 2? (d/dt = (1 - tau) d/dtau) *)
c = 1 + x; al = 1/z + 1; be = 1/z + 1; ga = c/z + 2;
Dtt[f_] := z (1 - tau) D[f, tau];
I0 = (1 - tau)^(-x/z) tau^(c/z + 1) F[tau];
pf = Dtt[Dtt[I0]] - (c/tau) Dtt[I0] + (x/tau) I0;
gauss = tau (1 - tau) F''[tau] + (ga - (al + be + 1) tau) F'[tau] - al be F[tau];
ratio = Simplify[pf/(I0/F[tau])];
(* express F'' through Gauss and test *)
red = Simplify[ratio /. F''[tau] -> Solve[gauss == 0, F''[tau]][[1, 1, 2]]];
Print["PF applied to I, with F'' eliminated by Gauss: ", red];
Print["PF/(Gauss operator) = ", Simplify[ratio/gauss]];
(* the same for the thimble of the regular critical point: other solution tau^(1 - gamma) F2, i.e. I2 = (1 - tau)^(-x/z) F~ with F~ the exponent-(1-gamma) solution *)
(* Delta_1 in terms of s*: x^2 W2 = x^3 (L1 - L2)/(x - c L2) with L2 = x s, L1 = c/tau - x s, s* root of x tau s^2 - c s + 1 = 0 *)
d1 = x^3 ((c/tau - x s) - x s)/(x - c x s);
w2 = (1 - tau s)^2 (-1/s^2 - x/(1 - s)^2 + tau^2/(1 - tau s)^2)/tau^2;
Print["x^2 W2 - Delta_1 modulo the critical equation: ", Simplify[PolynomialRemainder[Numerator[Together[x^2 w2 - d1]], x tau s^2 - c s + 1, s]]];
