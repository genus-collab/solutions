(* One branch point, x = -w^2/2, B = dw1 dw2/(w1-w2)^2 + beta dw1 dw2, y = a w + b w^3.
   DNOPS (2-31): R^{-1}(1/th) = -(1/th) L[dxi], dxi = (w^-2 + beta) dw  -> R^{-1}(z) = 1 - beta z.
   DNOPS (2-32): Delta^{-1/2} R^{-1}(1/th) = L[dy] = a + 3 b /th  -> a = Delta^{-1/2}, b = -beta a/3.
   Kernel of DNOPS (2-28) = manuscript eq:mod-kernel. Compute omega_{1,1}(z0) to first order in beta. *)
a = del^(-1/2); b = -beta a/3;
y[w_] := a w + b w^3;
Bint[z0_, w_] := (1/(z0 - w) - 1/(z0 + w)) + 2 beta w;  (* int_{-w}^{w} B(z0,.) / dz0 *)
K = (1/2) Bint[z0, w]/((y[w] - y[-w]) (-w));             (* / dx(w) = -w dw ; result per dz0 per dw^{-1} *)
w02 = -(1/(4 w^2) + beta);                              (* B(w,-w) per dw^2 *)
om11 = SeriesCoefficient[Series[K w02, {w, 0, 2}], -1];   (* residue at w=0 *)
om11 = Normal[Series[Simplify[om11], {beta, 0, 1}]];
(* Givental side with dilaton sign s: Delta^{1/2} [ c0 dxi + (1/24) d(dxi/dx) ], c0 = beta(-1/24 + s/24 + 1/2) *)
xi = -1/z0 + beta z0; dxi = D[xi, z0]; dxidx = D[xi, z0]/(-z0); ddxidx = D[dxidx, z0];
giv[s_] := Normal[Series[del^(1/2) (beta (-1/24 + s/24 + 1/2) dxi + (1/24) ddxidx), {beta, 0, 1}]];
Print["TR omega11 = ", Expand[om11]];
Print["PPZ sign (s=+1) minus TR: ", Simplify[giv[1] - om11]];
Print["printed DNOPS sign (s=-1) minus TR: ", Simplify[giv[-1] - om11]];
