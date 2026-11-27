Get["m0.mx"];
Print[Chop[N[Normal[Series[M0[[2,1]] /. t -> -psi, {psi,0,2}]],8]]];
Print[Chop[N[Normal[Series[M0[[1,2]] /. t -> -psi, {psi,0,2}]],8]]];
Print["M0 t^2 antisym part (1,2): ", Chop[N[Coefficient[M0[[1,2]]-M0[[2,1]], t, 2], 8]]];
