Get["rtr.wl"];
Print[N[Rinv[[1, 1]], 8]];
Print[N[Rinv[[2, 1]], 8]];
Print[N[Normal[Series[Mt[[1, 1]], {t, 0, 3}]], 8]];
Print[N[Normal[Series[Mt[[1, 2]], {t, 0, 3}]], 8]];
Print[N[Normal[Series[Mt[[2, 1]], {t, 0, 3}]], 8]];
