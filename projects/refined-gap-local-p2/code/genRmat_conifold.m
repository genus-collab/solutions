(* Plain-text copy of genRmat_conifold.mx (2026-11-27), written with Save; Get this file instead of the .mx. Symbols: C1cf, cf, invcf, q1sol, q2sol, tc *)

C1cf = 27*y*(1 + (11*(1 + 27*y))/9 + (109*(1 + 27*y)^2)/81 + 
      (9389*(1 + 27*y)^3)/6561 + (88351*(1 + 27*y)^4)/59049 + 
      (823187*(1 + 27*y)^5)/531441 + (68584051*(1 + 27*y)^6)/43046721 + 
      (632036939*(1 + 27*y)^7)/387420489 + (5805190621*(1 + 27*y)^8)/
       3486784401 + (4307941773209*(1 + 27*y)^9)/2541865828329)
 
cf[2] = SeriesData[tcv, 0, {(-7*\[Epsilon][1]^4 + 12*Qrrc[0, 0]*
         \[Epsilon][1]^5 - 132*Qrrc[0, 0]*\[Epsilon][1]^4*\[Epsilon][2] + 
        144*Qrrc[0, 0]^2*\[Epsilon][1]^5*\[Epsilon][2] - 
        10*\[Epsilon][1]^2*\[Epsilon][2]^2 - 120*Qrrc[0, 0]*\[Epsilon][1]^3*
         \[Epsilon][2]^2 - 576*Qrrc[0, 0]^2*\[Epsilon][1]^4*\[Epsilon][2]^2 + 
        576*Qrrc[0, 0]^3*\[Epsilon][1]^5*\[Epsilon][2]^2 + 
        120*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^2 - 
        120*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2]^3 - 
        1440*Qrrc[0, 0]^2*\[Epsilon][1]^3*\[Epsilon][2]^3 + 
        1728*Qrrc[0, 0]^3*\[Epsilon][1]^4*\[Epsilon][2]^3 + 
        360*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^3 - 7*\[Epsilon][2]^4 - 
        132*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^4 - 
        576*Qrrc[0, 0]^2*\[Epsilon][1]^2*\[Epsilon][2]^4 + 
        1728*Qrrc[0, 0]^3*\[Epsilon][1]^3*\[Epsilon][2]^4 + 
        360*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^4 + 
        12*Qrrc[0, 0]*\[Epsilon][2]^5 + 144*Qrrc[0, 0]^2*\[Epsilon][1]*
         \[Epsilon][2]^5 + 576*Qrrc[0, 0]^3*\[Epsilon][1]^2*\[Epsilon][2]^5 + 
        120*Qrrc[0, 1]*\[Epsilon][1]^2*\[Epsilon][2]^5)/
       (1920*\[Epsilon][1]*\[Epsilon][2]), (-(Qrrc[0, 0]*\[Epsilon][1]^5) - 
        3*Qrrc[0, 0]*\[Epsilon][1]^4*\[Epsilon][2] - 12*Qrrc[0, 0]^2*
         \[Epsilon][1]^5*\[Epsilon][2] - 4*Qrrc[0, 0]*\[Epsilon][1]^3*
         \[Epsilon][2]^2 - 36*Qrrc[0, 0]^2*\[Epsilon][1]^4*\[Epsilon][2]^2 - 
        48*Qrrc[0, 0]^3*\[Epsilon][1]^5*\[Epsilon][2]^2 - 
        10*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^2 - 
        4*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2]^3 - 
        48*Qrrc[0, 0]^2*\[Epsilon][1]^3*\[Epsilon][2]^3 - 
        144*Qrrc[0, 0]^3*\[Epsilon][1]^4*\[Epsilon][2]^3 - 
        30*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^3 - 
        3*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^4 - 36*Qrrc[0, 0]^2*
         \[Epsilon][1]^2*\[Epsilon][2]^4 - 144*Qrrc[0, 0]^3*\[Epsilon][1]^3*
         \[Epsilon][2]^4 - 30*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^4 - 
        Qrrc[0, 0]*\[Epsilon][2]^5 - 12*Qrrc[0, 0]^2*\[Epsilon][1]*
         \[Epsilon][2]^5 - 48*Qrrc[0, 0]^3*\[Epsilon][1]^2*\[Epsilon][2]^5 - 
        10*Qrrc[0, 1]*\[Epsilon][1]^2*\[Epsilon][2]^5)/
       (1440*\[Epsilon][1]*\[Epsilon][2]), 
      (-953*\[Epsilon][1]^4 + 492*Qrrc[0, 0]*\[Epsilon][1]^5 - 
        2844*\[Epsilon][1]^3*\[Epsilon][2] - 21396*Qrrc[0, 0]*\[Epsilon][1]^4*
         \[Epsilon][2] + 5904*Qrrc[0, 0]^2*\[Epsilon][1]^5*\[Epsilon][2] - 
        4718*\[Epsilon][1]^2*\[Epsilon][2]^2 - 55032*Qrrc[0, 0]*
         \[Epsilon][1]^3*\[Epsilon][2]^2 - 119520*Qrrc[0, 0]^2*
         \[Epsilon][1]^4*\[Epsilon][2]^2 + 23616*Qrrc[0, 0]^3*\[Epsilon][1]^5*
         \[Epsilon][2]^2 + 4920*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^2 - 
        2844*\[Epsilon][1]*\[Epsilon][2]^3 - 55032*Qrrc[0, 0]*\[Epsilon][1]^2*
         \[Epsilon][2]^3 - 250848*Qrrc[0, 0]^2*\[Epsilon][1]^3*
         \[Epsilon][2]^3 + 70848*Qrrc[0, 0]^3*\[Epsilon][1]^4*
         \[Epsilon][2]^3 + 14760*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^3 - 
        953*\[Epsilon][2]^4 - 21396*Qrrc[0, 0]*\[Epsilon][1]*
         \[Epsilon][2]^4 - 119520*Qrrc[0, 0]^2*\[Epsilon][1]^2*
         \[Epsilon][2]^4 + 70848*Qrrc[0, 0]^3*\[Epsilon][1]^3*
         \[Epsilon][2]^4 + 14760*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^4 + 
        492*Qrrc[0, 0]*\[Epsilon][2]^5 + 5904*Qrrc[0, 0]^2*\[Epsilon][1]*
         \[Epsilon][2]^5 + 23616*Qrrc[0, 0]^3*\[Epsilon][1]^2*
         \[Epsilon][2]^5 + 4920*Qrrc[0, 1]*\[Epsilon][1]^2*\[Epsilon][2]^5)/
       (1866240*\[Epsilon][1]*\[Epsilon][2]), 
      (1169*\[Epsilon][1]^4 - 1941*Qrrc[0, 0]*\[Epsilon][1]^5 + 
        7200*\[Epsilon][1]^3*\[Epsilon][2] + 22233*Qrrc[0, 0]*\[Epsilon][1]^4*
         \[Epsilon][2] - 23292*Qrrc[0, 0]^2*\[Epsilon][1]^5*\[Epsilon][2] + 
        12710*\[Epsilon][1]^2*\[Epsilon][2]^2 + 106692*Qrrc[0, 0]*
         \[Epsilon][1]^3*\[Epsilon][2]^2 + 98460*Qrrc[0, 0]^2*\[Epsilon][1]^4*
         \[Epsilon][2]^2 - 93168*Qrrc[0, 0]^3*\[Epsilon][1]^5*
         \[Epsilon][2]^2 - 19410*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^2 + 
        7200*\[Epsilon][1]*\[Epsilon][2]^3 + 106692*Qrrc[0, 0]*
         \[Epsilon][1]^2*\[Epsilon][2]^3 + 243504*Qrrc[0, 0]^2*
         \[Epsilon][1]^3*\[Epsilon][2]^3 - 279504*Qrrc[0, 0]^3*
         \[Epsilon][1]^4*\[Epsilon][2]^3 - 58230*Qrrc[0, 1]*\[Epsilon][1]^4*
         \[Epsilon][2]^3 + 1169*\[Epsilon][2]^4 + 22233*Qrrc[0, 0]*
         \[Epsilon][1]*\[Epsilon][2]^4 + 98460*Qrrc[0, 0]^2*\[Epsilon][1]^2*
         \[Epsilon][2]^4 - 279504*Qrrc[0, 0]^3*\[Epsilon][1]^3*
         \[Epsilon][2]^4 - 58230*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^4 - 
        1941*Qrrc[0, 0]*\[Epsilon][2]^5 - 23292*Qrrc[0, 0]^2*\[Epsilon][1]*
         \[Epsilon][2]^5 - 93168*Qrrc[0, 0]^3*\[Epsilon][1]^2*
         \[Epsilon][2]^5 - 19410*Qrrc[0, 1]*\[Epsilon][1]^2*\[Epsilon][2]^5)/
       (12597120*\[Epsilon][1]*\[Epsilon][2])}, -2, 2, 1]
 
cf[3] = SeriesData[tcv, 0, {(-465*\[Epsilon][1]^6 + 2260*Qrrc[0, 0]*
         \[Epsilon][1]^7 - 2880*Qrrc[0, 0]^2*\[Epsilon][1]^8 + 
        960*Qrrc[0, 0]^3*\[Epsilon][1]^9 + 200*Qrrc[0, 1]*\[Epsilon][1]^9 - 
        10660*Qrrc[0, 0]*\[Epsilon][1]^6*\[Epsilon][2] + 
        40080*Qrrc[0, 0]^2*\[Epsilon][1]^7*\[Epsilon][2] - 
        40320*Qrrc[0, 0]^3*\[Epsilon][1]^8*\[Epsilon][2] - 
        1200*Qrrc[0, 1]*\[Epsilon][1]^8*\[Epsilon][2] + 
        11520*Qrrc[0, 0]^4*\[Epsilon][1]^9*\[Epsilon][2] + 
        2400*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^9*\[Epsilon][2] - 
        735*\[Epsilon][1]^4*\[Epsilon][2]^2 - 7980*Qrrc[0, 0]*\[Epsilon][1]^5*
         \[Epsilon][2]^2 - 52800*Qrrc[0, 0]^2*\[Epsilon][1]^6*
         \[Epsilon][2]^2 + 221760*Qrrc[0, 0]^3*\[Epsilon][1]^7*
         \[Epsilon][2]^2 + 13800*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^2 - 
        207360*Qrrc[0, 0]^4*\[Epsilon][1]^8*\[Epsilon][2]^2 - 
        43200*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^8*\[Epsilon][2]^2 + 
        55296*Qrrc[0, 0]^5*\[Epsilon][1]^9*\[Epsilon][2]^2 + 
        28800*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^9*\[Epsilon][2]^2 + 
        420*Qrrc[0, 2]*\[Epsilon][1]^9*\[Epsilon][2]^2 - 
        13860*Qrrc[0, 0]*\[Epsilon][1]^4*\[Epsilon][2]^3 - 
        156240*Qrrc[0, 0]^2*\[Epsilon][1]^5*\[Epsilon][2]^3 + 
        298560*Qrrc[0, 0]^3*\[Epsilon][1]^6*\[Epsilon][2]^3 + 
        41800*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^3 + 
        115200*Qrrc[0, 0]^4*\[Epsilon][1]^7*\[Epsilon][2]^3 + 
        24000*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^3 - 
        331776*Qrrc[0, 0]^5*\[Epsilon][1]^8*\[Epsilon][2]^3 - 
        172800*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^8*\[Epsilon][2]^3 - 
        2520*Qrrc[0, 2]*\[Epsilon][1]^8*\[Epsilon][2]^3 + 
        110592*Qrrc[0, 0]^6*\[Epsilon][1]^9*\[Epsilon][2]^3 + 
        115200*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^9*\[Epsilon][2]^3 - 
        24000*Qrrc[0, 1]^2*\[Epsilon][1]^9*\[Epsilon][2]^3 + 
        5040*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^9*\[Epsilon][2]^3 - 
        735*\[Epsilon][1]^2*\[Epsilon][2]^4 - 13860*Qrrc[0, 0]*
         \[Epsilon][1]^3*\[Epsilon][2]^4 - 120960*Qrrc[0, 0]^2*
         \[Epsilon][1]^4*\[Epsilon][2]^4 - 504000*Qrrc[0, 0]^3*
         \[Epsilon][1]^5*\[Epsilon][2]^4 + 46200*Qrrc[0, 1]*\[Epsilon][1]^5*
         \[Epsilon][2]^4 + 2511360*Qrrc[0, 0]^4*\[Epsilon][1]^6*
         \[Epsilon][2]^4 + 523200*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^6*
         \[Epsilon][2]^4 - 2433024*Qrrc[0, 0]^5*\[Epsilon][1]^7*
         \[Epsilon][2]^4 - 1267200*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^7*
         \[Epsilon][2]^4 - 18480*Qrrc[0, 2]*\[Epsilon][1]^7*\[Epsilon][2]^4 + 
        663552*Qrrc[0, 0]^6*\[Epsilon][1]^8*\[Epsilon][2]^4 + 
        691200*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^8*\[Epsilon][2]^4 - 
        144000*Qrrc[0, 1]^2*\[Epsilon][1]^8*\[Epsilon][2]^4 + 
        30240*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^8*\[Epsilon][2]^4 - 
        7980*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2]^5 - 
        156240*Qrrc[0, 0]^2*\[Epsilon][1]^3*\[Epsilon][2]^5 - 
        504000*Qrrc[0, 0]^3*\[Epsilon][1]^4*\[Epsilon][2]^5 + 
        46200*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^5 + 
        4354560*Qrrc[0, 0]^4*\[Epsilon][1]^5*\[Epsilon][2]^5 + 
        907200*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^5 - 
        5253120*Qrrc[0, 0]^5*\[Epsilon][1]^6*\[Epsilon][2]^5 - 
        2736000*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^5 - 
        39900*Qrrc[0, 2]*\[Epsilon][1]^6*\[Epsilon][2]^5 + 
        1658880*Qrrc[0, 0]^6*\[Epsilon][1]^7*\[Epsilon][2]^5 + 
        1728000*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^5 - 
        360000*Qrrc[0, 1]^2*\[Epsilon][1]^7*\[Epsilon][2]^5 + 
        75600*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^7*\[Epsilon][2]^5 - 
        465*\[Epsilon][2]^6 - 10660*Qrrc[0, 0]*\[Epsilon][1]*
         \[Epsilon][2]^6 - 52800*Qrrc[0, 0]^2*\[Epsilon][1]^2*
         \[Epsilon][2]^6 + 298560*Qrrc[0, 0]^3*\[Epsilon][1]^3*
         \[Epsilon][2]^6 + 41800*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^6 + 
        2511360*Qrrc[0, 0]^4*\[Epsilon][1]^4*\[Epsilon][2]^6 + 
        523200*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^6 - 
        5253120*Qrrc[0, 0]^5*\[Epsilon][1]^5*\[Epsilon][2]^6 - 
        2736000*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^6 - 
        39900*Qrrc[0, 2]*\[Epsilon][1]^5*\[Epsilon][2]^6 + 
        2211840*Qrrc[0, 0]^6*\[Epsilon][1]^6*\[Epsilon][2]^6 + 
        2304000*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^6 - 
        480000*Qrrc[0, 1]^2*\[Epsilon][1]^6*\[Epsilon][2]^6 + 
        100800*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^6*\[Epsilon][2]^6 + 
        2260*Qrrc[0, 0]*\[Epsilon][2]^7 + 40080*Qrrc[0, 0]^2*\[Epsilon][1]*
         \[Epsilon][2]^7 + 221760*Qrrc[0, 0]^3*\[Epsilon][1]^2*
         \[Epsilon][2]^7 + 13800*Qrrc[0, 1]*\[Epsilon][1]^2*\[Epsilon][2]^7 + 
        115200*Qrrc[0, 0]^4*\[Epsilon][1]^3*\[Epsilon][2]^7 + 
        24000*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^7 - 
        2433024*Qrrc[0, 0]^5*\[Epsilon][1]^4*\[Epsilon][2]^7 - 
        1267200*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^7 - 
        18480*Qrrc[0, 2]*\[Epsilon][1]^4*\[Epsilon][2]^7 + 
        1658880*Qrrc[0, 0]^6*\[Epsilon][1]^5*\[Epsilon][2]^7 + 
        1728000*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^7 - 
        360000*Qrrc[0, 1]^2*\[Epsilon][1]^5*\[Epsilon][2]^7 + 
        75600*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^5*\[Epsilon][2]^7 - 
        2880*Qrrc[0, 0]^2*\[Epsilon][2]^8 - 40320*Qrrc[0, 0]^3*\[Epsilon][1]*
         \[Epsilon][2]^8 - 1200*Qrrc[0, 1]*\[Epsilon][1]*\[Epsilon][2]^8 - 
        207360*Qrrc[0, 0]^4*\[Epsilon][1]^2*\[Epsilon][2]^8 - 
        43200*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^2*\[Epsilon][2]^8 - 
        331776*Qrrc[0, 0]^5*\[Epsilon][1]^3*\[Epsilon][2]^8 - 
        172800*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^8 - 
        2520*Qrrc[0, 2]*\[Epsilon][1]^3*\[Epsilon][2]^8 + 
        663552*Qrrc[0, 0]^6*\[Epsilon][1]^4*\[Epsilon][2]^8 + 
        691200*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^8 - 
        144000*Qrrc[0, 1]^2*\[Epsilon][1]^4*\[Epsilon][2]^8 + 
        30240*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^4*\[Epsilon][2]^8 + 
        960*Qrrc[0, 0]^3*\[Epsilon][2]^9 + 200*Qrrc[0, 1]*\[Epsilon][2]^9 + 
        11520*Qrrc[0, 0]^4*\[Epsilon][1]*\[Epsilon][2]^9 + 
        2400*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]*\[Epsilon][2]^9 + 
        55296*Qrrc[0, 0]^5*\[Epsilon][1]^2*\[Epsilon][2]^9 + 
        28800*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^2*\[Epsilon][2]^9 + 
        420*Qrrc[0, 2]*\[Epsilon][1]^2*\[Epsilon][2]^9 + 
        110592*Qrrc[0, 0]^6*\[Epsilon][1]^3*\[Epsilon][2]^9 + 
        115200*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^9 - 
        24000*Qrrc[0, 1]^2*\[Epsilon][1]^3*\[Epsilon][2]^9 + 
        5040*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^3*\[Epsilon][2]^9)/
       (268800*\[Epsilon][1]*\[Epsilon][2]), 
      (-205*Qrrc[0, 0]*\[Epsilon][1]^7 + 720*Qrrc[0, 0]^2*\[Epsilon][1]^8 - 
        480*Qrrc[0, 0]^3*\[Epsilon][1]^9 - 100*Qrrc[0, 1]*\[Epsilon][1]^9 - 
        615*Qrrc[0, 0]*\[Epsilon][1]^6*\[Epsilon][2] - 
        1440*Qrrc[0, 0]^2*\[Epsilon][1]^7*\[Epsilon][2] + 
        8640*Qrrc[0, 0]^3*\[Epsilon][1]^8*\[Epsilon][2] - 
        5760*Qrrc[0, 0]^4*\[Epsilon][1]^9*\[Epsilon][2] - 
        1200*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^9*\[Epsilon][2] - 
        925*Qrrc[0, 0]*\[Epsilon][1]^5*\[Epsilon][2]^2 - 
        9180*Qrrc[0, 0]^2*\[Epsilon][1]^6*\[Epsilon][2]^2 + 
        9120*Qrrc[0, 0]^3*\[Epsilon][1]^7*\[Epsilon][2]^2 - 
        650*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^2 + 
        34560*Qrrc[0, 0]^4*\[Epsilon][1]^8*\[Epsilon][2]^2 + 
        7200*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^8*\[Epsilon][2]^2 - 
        27648*Qrrc[0, 0]^5*\[Epsilon][1]^9*\[Epsilon][2]^2 - 
        14400*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^9*\[Epsilon][2]^2 - 
        210*Qrrc[0, 2]*\[Epsilon][1]^9*\[Epsilon][2]^2 - 
        1135*Qrrc[0, 0]*\[Epsilon][1]^4*\[Epsilon][2]^3 - 
        16560*Qrrc[0, 0]^2*\[Epsilon][1]^5*\[Epsilon][2]^3 - 
        34320*Qrrc[0, 0]^3*\[Epsilon][1]^6*\[Epsilon][2]^3 - 
        2650*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^3 + 
        129600*Qrrc[0, 0]^4*\[Epsilon][1]^7*\[Epsilon][2]^3 + 
        27000*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^3 - 
        55296*Qrrc[0, 0]^6*\[Epsilon][1]^9*\[Epsilon][2]^3 - 
        57600*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^9*\[Epsilon][2]^3 + 
        12000*Qrrc[0, 1]^2*\[Epsilon][1]^9*\[Epsilon][2]^3 - 
        2520*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^9*\[Epsilon][2]^3 - 
        1135*Qrrc[0, 0]*\[Epsilon][1]^3*\[Epsilon][2]^4 - 
        19080*Qrrc[0, 0]^2*\[Epsilon][1]^4*\[Epsilon][2]^4 - 
        83760*Qrrc[0, 0]^3*\[Epsilon][1]^5*\[Epsilon][2]^4 - 
        3800*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^4 + 
        115200*Qrrc[0, 0]^4*\[Epsilon][1]^6*\[Epsilon][2]^4 + 
        24000*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^4 + 
        387072*Qrrc[0, 0]^5*\[Epsilon][1]^7*\[Epsilon][2]^4 + 
        201600*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^4 + 
        2940*Qrrc[0, 2]*\[Epsilon][1]^7*\[Epsilon][2]^4 - 
        331776*Qrrc[0, 0]^6*\[Epsilon][1]^8*\[Epsilon][2]^4 - 
        345600*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^8*\[Epsilon][2]^4 + 
        72000*Qrrc[0, 1]^2*\[Epsilon][1]^8*\[Epsilon][2]^4 - 
        15120*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^8*\[Epsilon][2]^4 - 
        925*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2]^5 - 
        16560*Qrrc[0, 0]^2*\[Epsilon][1]^3*\[Epsilon][2]^5 - 
        83760*Qrrc[0, 0]^3*\[Epsilon][1]^4*\[Epsilon][2]^5 - 
        3800*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^5 + 
        51840*Qrrc[0, 0]^4*\[Epsilon][1]^5*\[Epsilon][2]^5 + 
        10800*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^5 + 
        967680*Qrrc[0, 0]^5*\[Epsilon][1]^6*\[Epsilon][2]^5 + 
        504000*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^5 + 
        7350*Qrrc[0, 2]*\[Epsilon][1]^6*\[Epsilon][2]^5 - 
        829440*Qrrc[0, 0]^6*\[Epsilon][1]^7*\[Epsilon][2]^5 - 
        864000*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^5 + 
        180000*Qrrc[0, 1]^2*\[Epsilon][1]^7*\[Epsilon][2]^5 - 
        37800*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^7*\[Epsilon][2]^5 - 
        615*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^6 - 
        9180*Qrrc[0, 0]^2*\[Epsilon][1]^2*\[Epsilon][2]^6 - 
        34320*Qrrc[0, 0]^3*\[Epsilon][1]^3*\[Epsilon][2]^6 - 
        2650*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^6 + 
        115200*Qrrc[0, 0]^4*\[Epsilon][1]^4*\[Epsilon][2]^6 + 
        24000*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^6 + 
        967680*Qrrc[0, 0]^5*\[Epsilon][1]^5*\[Epsilon][2]^6 + 
        504000*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^6 + 
        7350*Qrrc[0, 2]*\[Epsilon][1]^5*\[Epsilon][2]^6 - 
        1105920*Qrrc[0, 0]^6*\[Epsilon][1]^6*\[Epsilon][2]^6 - 
        1152000*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^6 + 
        240000*Qrrc[0, 1]^2*\[Epsilon][1]^6*\[Epsilon][2]^6 - 
        50400*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^6*\[Epsilon][2]^6 - 
        205*Qrrc[0, 0]*\[Epsilon][2]^7 - 1440*Qrrc[0, 0]^2*\[Epsilon][1]*
         \[Epsilon][2]^7 + 9120*Qrrc[0, 0]^3*\[Epsilon][1]^2*
         \[Epsilon][2]^7 - 650*Qrrc[0, 1]*\[Epsilon][1]^2*\[Epsilon][2]^7 + 
        129600*Qrrc[0, 0]^4*\[Epsilon][1]^3*\[Epsilon][2]^7 + 
        27000*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^7 + 
        387072*Qrrc[0, 0]^5*\[Epsilon][1]^4*\[Epsilon][2]^7 + 
        201600*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^7 + 
        2940*Qrrc[0, 2]*\[Epsilon][1]^4*\[Epsilon][2]^7 - 
        829440*Qrrc[0, 0]^6*\[Epsilon][1]^5*\[Epsilon][2]^7 - 
        864000*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^7 + 
        180000*Qrrc[0, 1]^2*\[Epsilon][1]^5*\[Epsilon][2]^7 - 
        37800*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^5*\[Epsilon][2]^7 + 
        720*Qrrc[0, 0]^2*\[Epsilon][2]^8 + 8640*Qrrc[0, 0]^3*\[Epsilon][1]*
         \[Epsilon][2]^8 + 34560*Qrrc[0, 0]^4*\[Epsilon][1]^2*
         \[Epsilon][2]^8 + 7200*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^2*
         \[Epsilon][2]^8 - 331776*Qrrc[0, 0]^6*\[Epsilon][1]^4*
         \[Epsilon][2]^8 - 345600*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^4*
         \[Epsilon][2]^8 + 72000*Qrrc[0, 1]^2*\[Epsilon][1]^4*
         \[Epsilon][2]^8 - 15120*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^4*
         \[Epsilon][2]^8 - 480*Qrrc[0, 0]^3*\[Epsilon][2]^9 - 
        100*Qrrc[0, 1]*\[Epsilon][2]^9 - 5760*Qrrc[0, 0]^4*\[Epsilon][1]*
         \[Epsilon][2]^9 - 1200*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]*
         \[Epsilon][2]^9 - 27648*Qrrc[0, 0]^5*\[Epsilon][1]^2*
         \[Epsilon][2]^9 - 14400*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^2*
         \[Epsilon][2]^9 - 210*Qrrc[0, 2]*\[Epsilon][1]^2*\[Epsilon][2]^9 - 
        55296*Qrrc[0, 0]^6*\[Epsilon][1]^3*\[Epsilon][2]^9 - 
        57600*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^9 + 
        12000*Qrrc[0, 1]^2*\[Epsilon][1]^3*\[Epsilon][2]^9 - 
        2520*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^3*\[Epsilon][2]^9)/
       (604800*\[Epsilon][1]*\[Epsilon][2]), 
      (50*Qrrc[0, 0]*\[Epsilon][1]^7 - 1440*Qrrc[0, 0]^2*\[Epsilon][1]^8 + 
        2640*Qrrc[0, 0]^3*\[Epsilon][1]^9 + 550*Qrrc[0, 1]*\[Epsilon][1]^9 - 
        155*Qrrc[0, 0]*\[Epsilon][1]^6*\[Epsilon][2] - 
        5340*Qrrc[0, 0]^2*\[Epsilon][1]^7*\[Epsilon][2] - 
        7200*Qrrc[0, 0]^3*\[Epsilon][1]^8*\[Epsilon][2] + 
        2100*Qrrc[0, 1]*\[Epsilon][1]^8*\[Epsilon][2] + 
        31680*Qrrc[0, 0]^4*\[Epsilon][1]^9*\[Epsilon][2] + 
        6600*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^9*\[Epsilon][2] - 
        795*Qrrc[0, 0]*\[Epsilon][1]^5*\[Epsilon][2]^2 - 
        13440*Qrrc[0, 0]^2*\[Epsilon][1]^6*\[Epsilon][2]^2 - 
        55440*Qrrc[0, 0]^3*\[Epsilon][1]^7*\[Epsilon][2]^2 + 
        3300*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^2 + 
        51840*Qrrc[0, 0]^4*\[Epsilon][1]^8*\[Epsilon][2]^2 + 
        10800*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^8*\[Epsilon][2]^2 + 
        152064*Qrrc[0, 0]^5*\[Epsilon][1]^9*\[Epsilon][2]^2 + 
        79200*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^9*\[Epsilon][2]^2 + 
        1155*Qrrc[0, 2]*\[Epsilon][1]^9*\[Epsilon][2]^2 - 
        1260*Qrrc[0, 0]*\[Epsilon][1]^4*\[Epsilon][2]^3 - 
        24660*Qrrc[0, 0]^2*\[Epsilon][1]^5*\[Epsilon][2]^3 - 
        138000*Qrrc[0, 0]^3*\[Epsilon][1]^6*\[Epsilon][2]^3 + 
        200*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^3 - 
        123840*Qrrc[0, 0]^4*\[Epsilon][1]^7*\[Epsilon][2]^3 - 
        25800*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^3 + 
        580608*Qrrc[0, 0]^5*\[Epsilon][1]^8*\[Epsilon][2]^3 + 
        302400*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^8*\[Epsilon][2]^3 + 
        4410*Qrrc[0, 2]*\[Epsilon][1]^8*\[Epsilon][2]^3 + 
        304128*Qrrc[0, 0]^6*\[Epsilon][1]^9*\[Epsilon][2]^3 + 
        316800*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^9*\[Epsilon][2]^3 - 
        66000*Qrrc[0, 1]^2*\[Epsilon][1]^9*\[Epsilon][2]^3 + 
        13860*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^9*\[Epsilon][2]^3 - 
        1260*Qrrc[0, 0]*\[Epsilon][1]^3*\[Epsilon][2]^4 - 
        30240*Qrrc[0, 0]^2*\[Epsilon][1]^4*\[Epsilon][2]^4 - 
        210960*Qrrc[0, 0]^3*\[Epsilon][1]^5*\[Epsilon][2]^4 - 
        6150*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^4 - 
        455040*Qrrc[0, 0]^4*\[Epsilon][1]^6*\[Epsilon][2]^4 - 
        94800*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^4 + 
        774144*Qrrc[0, 0]^5*\[Epsilon][1]^7*\[Epsilon][2]^4 + 
        403200*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^4 + 
        5880*Qrrc[0, 2]*\[Epsilon][1]^7*\[Epsilon][2]^4 + 
        1824768*Qrrc[0, 0]^6*\[Epsilon][1]^8*\[Epsilon][2]^4 + 
        1900800*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^8*\[Epsilon][2]^4 - 
        396000*Qrrc[0, 1]^2*\[Epsilon][1]^8*\[Epsilon][2]^4 + 
        83160*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^8*\[Epsilon][2]^4 - 
        795*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2]^5 - 
        24660*Qrrc[0, 0]^2*\[Epsilon][1]^3*\[Epsilon][2]^5 - 
        210960*Qrrc[0, 0]^3*\[Epsilon][1]^4*\[Epsilon][2]^5 - 
        6150*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^5 - 
        622080*Qrrc[0, 0]^4*\[Epsilon][1]^5*\[Epsilon][2]^5 - 
        129600*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^5 + 
        483840*Qrrc[0, 0]^5*\[Epsilon][1]^6*\[Epsilon][2]^5 + 
        252000*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^5 + 
        3675*Qrrc[0, 2]*\[Epsilon][1]^6*\[Epsilon][2]^5 + 
        4561920*Qrrc[0, 0]^6*\[Epsilon][1]^7*\[Epsilon][2]^5 + 
        4752000*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^5 - 
        990000*Qrrc[0, 1]^2*\[Epsilon][1]^7*\[Epsilon][2]^5 + 
        207900*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^7*\[Epsilon][2]^5 - 
        155*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^6 - 
        13440*Qrrc[0, 0]^2*\[Epsilon][1]^2*\[Epsilon][2]^6 - 
        138000*Qrrc[0, 0]^3*\[Epsilon][1]^3*\[Epsilon][2]^6 + 
        200*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^6 - 
        455040*Qrrc[0, 0]^4*\[Epsilon][1]^4*\[Epsilon][2]^6 - 
        94800*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^6 + 
        483840*Qrrc[0, 0]^5*\[Epsilon][1]^5*\[Epsilon][2]^6 + 
        252000*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^6 + 
        3675*Qrrc[0, 2]*\[Epsilon][1]^5*\[Epsilon][2]^6 + 
        6082560*Qrrc[0, 0]^6*\[Epsilon][1]^6*\[Epsilon][2]^6 + 
        6336000*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^6 - 
        1320000*Qrrc[0, 1]^2*\[Epsilon][1]^6*\[Epsilon][2]^6 + 
        277200*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^6*\[Epsilon][2]^6 + 
        50*Qrrc[0, 0]*\[Epsilon][2]^7 - 5340*Qrrc[0, 0]^2*\[Epsilon][1]*
         \[Epsilon][2]^7 - 55440*Qrrc[0, 0]^3*\[Epsilon][1]^2*
         \[Epsilon][2]^7 + 3300*Qrrc[0, 1]*\[Epsilon][1]^2*\[Epsilon][2]^7 - 
        123840*Qrrc[0, 0]^4*\[Epsilon][1]^3*\[Epsilon][2]^7 - 
        25800*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^7 + 
        774144*Qrrc[0, 0]^5*\[Epsilon][1]^4*\[Epsilon][2]^7 + 
        403200*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^7 + 
        5880*Qrrc[0, 2]*\[Epsilon][1]^4*\[Epsilon][2]^7 + 
        4561920*Qrrc[0, 0]^6*\[Epsilon][1]^5*\[Epsilon][2]^7 + 
        4752000*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^7 - 
        990000*Qrrc[0, 1]^2*\[Epsilon][1]^5*\[Epsilon][2]^7 + 
        207900*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^5*\[Epsilon][2]^7 - 
        1440*Qrrc[0, 0]^2*\[Epsilon][2]^8 - 7200*Qrrc[0, 0]^3*\[Epsilon][1]*
         \[Epsilon][2]^8 + 2100*Qrrc[0, 1]*\[Epsilon][1]*\[Epsilon][2]^8 + 
        51840*Qrrc[0, 0]^4*\[Epsilon][1]^2*\[Epsilon][2]^8 + 
        10800*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^2*\[Epsilon][2]^8 + 
        580608*Qrrc[0, 0]^5*\[Epsilon][1]^3*\[Epsilon][2]^8 + 
        302400*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^8 + 
        4410*Qrrc[0, 2]*\[Epsilon][1]^3*\[Epsilon][2]^8 + 
        1824768*Qrrc[0, 0]^6*\[Epsilon][1]^4*\[Epsilon][2]^8 + 
        1900800*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^8 - 
        396000*Qrrc[0, 1]^2*\[Epsilon][1]^4*\[Epsilon][2]^8 + 
        83160*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^4*\[Epsilon][2]^8 + 
        2640*Qrrc[0, 0]^3*\[Epsilon][2]^9 + 550*Qrrc[0, 1]*\[Epsilon][2]^9 + 
        31680*Qrrc[0, 0]^4*\[Epsilon][1]*\[Epsilon][2]^9 + 
        6600*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]*\[Epsilon][2]^9 + 
        152064*Qrrc[0, 0]^5*\[Epsilon][1]^2*\[Epsilon][2]^9 + 
        79200*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^2*\[Epsilon][2]^9 + 
        1155*Qrrc[0, 2]*\[Epsilon][1]^2*\[Epsilon][2]^9 + 
        304128*Qrrc[0, 0]^6*\[Epsilon][1]^3*\[Epsilon][2]^9 + 
        316800*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^9 - 
        66000*Qrrc[0, 1]^2*\[Epsilon][1]^3*\[Epsilon][2]^9 + 
        13860*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^3*\[Epsilon][2]^9)/
       (32659200*\[Epsilon][1]*\[Epsilon][2]), 
      (-2020*Qrrc[0, 0]*\[Epsilon][1]^7 - 1440*Qrrc[0, 0]^2*\[Epsilon][1]^8 - 
        1680*Qrrc[0, 0]^3*\[Epsilon][1]^9 - 350*Qrrc[0, 1]*\[Epsilon][1]^9 - 
        17030*Qrrc[0, 0]*\[Epsilon][1]^6*\[Epsilon][2] - 
        56100*Qrrc[0, 0]^2*\[Epsilon][1]^7*\[Epsilon][2] - 
        33120*Qrrc[0, 0]^3*\[Epsilon][1]^8*\[Epsilon][2] - 
        3300*Qrrc[0, 1]*\[Epsilon][1]^8*\[Epsilon][2] - 
        20160*Qrrc[0, 0]^4*\[Epsilon][1]^9*\[Epsilon][2] - 
        4200*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^9*\[Epsilon][2] - 
        43320*Qrrc[0, 0]*\[Epsilon][1]^5*\[Epsilon][2]^2 - 
        320700*Qrrc[0, 0]^2*\[Epsilon][1]^6*\[Epsilon][2]^2 - 
        534960*Qrrc[0, 0]^3*\[Epsilon][1]^7*\[Epsilon][2]^2 - 
        31800*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^2 - 
        259200*Qrrc[0, 0]^4*\[Epsilon][1]^8*\[Epsilon][2]^2 - 
        54000*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^8*\[Epsilon][2]^2 - 
        96768*Qrrc[0, 0]^5*\[Epsilon][1]^9*\[Epsilon][2]^2 - 
        50400*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^9*\[Epsilon][2]^2 - 
        735*Qrrc[0, 2]*\[Epsilon][1]^9*\[Epsilon][2]^2 - 
        58950*Qrrc[0, 0]*\[Epsilon][1]^4*\[Epsilon][2]^3 - 
        719640*Qrrc[0, 0]^2*\[Epsilon][1]^5*\[Epsilon][2]^3 - 
        2313120*Qrrc[0, 0]^3*\[Epsilon][1]^6*\[Epsilon][2]^3 - 
        191050*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^3 - 
        2197440*Qrrc[0, 0]^4*\[Epsilon][1]^7*\[Epsilon][2]^3 - 
        457800*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^3 - 
        912384*Qrrc[0, 0]^5*\[Epsilon][1]^8*\[Epsilon][2]^3 - 
        475200*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^8*\[Epsilon][2]^3 - 
        6930*Qrrc[0, 2]*\[Epsilon][1]^8*\[Epsilon][2]^3 - 
        193536*Qrrc[0, 0]^6*\[Epsilon][1]^9*\[Epsilon][2]^3 - 
        201600*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^9*\[Epsilon][2]^3 + 
        42000*Qrrc[0, 1]^2*\[Epsilon][1]^9*\[Epsilon][2]^3 - 
        8820*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^9*\[Epsilon][2]^3 - 
        58950*Qrrc[0, 0]*\[Epsilon][1]^3*\[Epsilon][2]^4 - 
        907200*Qrrc[0, 0]^2*\[Epsilon][1]^4*\[Epsilon][2]^4 - 
        4481280*Qrrc[0, 0]^3*\[Epsilon][1]^5*\[Epsilon][2]^4 - 
        434100*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^4 - 
        6779520*Qrrc[0, 0]^4*\[Epsilon][1]^6*\[Epsilon][2]^4 - 
        1412400*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^4 - 
        3207168*Qrrc[0, 0]^5*\[Epsilon][1]^7*\[Epsilon][2]^4 - 
        1670400*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^4 - 
        24360*Qrrc[0, 2]*\[Epsilon][1]^7*\[Epsilon][2]^4 - 
        1161216*Qrrc[0, 0]^6*\[Epsilon][1]^8*\[Epsilon][2]^4 - 
        1209600*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^8*\[Epsilon][2]^4 + 
        252000*Qrrc[0, 1]^2*\[Epsilon][1]^8*\[Epsilon][2]^4 - 
        52920*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^8*\[Epsilon][2]^4 - 
        43320*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2]^5 - 
        719640*Qrrc[0, 0]^2*\[Epsilon][1]^3*\[Epsilon][2]^5 - 
        4481280*Qrrc[0, 0]^3*\[Epsilon][1]^4*\[Epsilon][2]^5 - 
        434100*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^5 - 
        9642240*Qrrc[0, 0]^4*\[Epsilon][1]^5*\[Epsilon][2]^5 - 
        2008800*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^5 - 
        5736960*Qrrc[0, 0]^5*\[Epsilon][1]^6*\[Epsilon][2]^5 - 
        2988000*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^5 - 
        43575*Qrrc[0, 2]*\[Epsilon][1]^6*\[Epsilon][2]^5 - 
        2903040*Qrrc[0, 0]^6*\[Epsilon][1]^7*\[Epsilon][2]^5 - 
        3024000*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^5 + 
        630000*Qrrc[0, 1]^2*\[Epsilon][1]^7*\[Epsilon][2]^5 - 
        132300*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^7*\[Epsilon][2]^5 - 
        17030*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^6 - 
        320700*Qrrc[0, 0]^2*\[Epsilon][1]^2*\[Epsilon][2]^6 - 
        2313120*Qrrc[0, 0]^3*\[Epsilon][1]^3*\[Epsilon][2]^6 - 
        191050*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^6 - 
        6779520*Qrrc[0, 0]^4*\[Epsilon][1]^4*\[Epsilon][2]^6 - 
        1412400*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^6 - 
        5736960*Qrrc[0, 0]^5*\[Epsilon][1]^5*\[Epsilon][2]^6 - 
        2988000*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^6 - 
        43575*Qrrc[0, 2]*\[Epsilon][1]^5*\[Epsilon][2]^6 - 
        3870720*Qrrc[0, 0]^6*\[Epsilon][1]^6*\[Epsilon][2]^6 - 
        4032000*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^6 + 
        840000*Qrrc[0, 1]^2*\[Epsilon][1]^6*\[Epsilon][2]^6 - 
        176400*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^6*\[Epsilon][2]^6 - 
        2020*Qrrc[0, 0]*\[Epsilon][2]^7 - 56100*Qrrc[0, 0]^2*\[Epsilon][1]*
         \[Epsilon][2]^7 - 534960*Qrrc[0, 0]^3*\[Epsilon][1]^2*
         \[Epsilon][2]^7 - 31800*Qrrc[0, 1]*\[Epsilon][1]^2*\[Epsilon][2]^7 - 
        2197440*Qrrc[0, 0]^4*\[Epsilon][1]^3*\[Epsilon][2]^7 - 
        457800*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^7 - 
        3207168*Qrrc[0, 0]^5*\[Epsilon][1]^4*\[Epsilon][2]^7 - 
        1670400*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^7 - 
        24360*Qrrc[0, 2]*\[Epsilon][1]^4*\[Epsilon][2]^7 - 
        2903040*Qrrc[0, 0]^6*\[Epsilon][1]^5*\[Epsilon][2]^7 - 
        3024000*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^7 + 
        630000*Qrrc[0, 1]^2*\[Epsilon][1]^5*\[Epsilon][2]^7 - 
        132300*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^5*\[Epsilon][2]^7 - 
        1440*Qrrc[0, 0]^2*\[Epsilon][2]^8 - 33120*Qrrc[0, 0]^3*\[Epsilon][1]*
         \[Epsilon][2]^8 - 3300*Qrrc[0, 1]*\[Epsilon][1]*\[Epsilon][2]^8 - 
        259200*Qrrc[0, 0]^4*\[Epsilon][1]^2*\[Epsilon][2]^8 - 
        54000*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^2*\[Epsilon][2]^8 - 
        912384*Qrrc[0, 0]^5*\[Epsilon][1]^3*\[Epsilon][2]^8 - 
        475200*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^8 - 
        6930*Qrrc[0, 2]*\[Epsilon][1]^3*\[Epsilon][2]^8 - 
        1161216*Qrrc[0, 0]^6*\[Epsilon][1]^4*\[Epsilon][2]^8 - 
        1209600*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^8 + 
        252000*Qrrc[0, 1]^2*\[Epsilon][1]^4*\[Epsilon][2]^8 - 
        52920*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^4*\[Epsilon][2]^8 - 
        1680*Qrrc[0, 0]^3*\[Epsilon][2]^9 - 350*Qrrc[0, 1]*\[Epsilon][2]^9 - 
        20160*Qrrc[0, 0]^4*\[Epsilon][1]*\[Epsilon][2]^9 - 
        4200*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]*\[Epsilon][2]^9 - 
        96768*Qrrc[0, 0]^5*\[Epsilon][1]^2*\[Epsilon][2]^9 - 
        50400*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^2*\[Epsilon][2]^9 - 
        735*Qrrc[0, 2]*\[Epsilon][1]^2*\[Epsilon][2]^9 - 
        193536*Qrrc[0, 0]^6*\[Epsilon][1]^3*\[Epsilon][2]^9 - 
        201600*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^9 + 
        42000*Qrrc[0, 1]^2*\[Epsilon][1]^3*\[Epsilon][2]^9 - 
        8820*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^3*\[Epsilon][2]^9)/
       (440899200*\[Epsilon][1]*\[Epsilon][2]), 
      (1165215*\[Epsilon][1]^6 - 8135260*Qrrc[0, 0]*\[Epsilon][1]^7 + 
        19615680*Qrrc[0, 0]^2*\[Epsilon][1]^8 - 6394560*Qrrc[0, 0]^3*
         \[Epsilon][1]^9 - 1332200*Qrrc[0, 1]*\[Epsilon][1]^9 + 
        27025200*\[Epsilon][1]^5*\[Epsilon][2] - 20616140*Qrrc[0, 0]*
         \[Epsilon][1]^6*\[Epsilon][2] - 98766480*Qrrc[0, 0]^2*
         \[Epsilon][1]^7*\[Epsilon][2] + 275483520*Qrrc[0, 0]^3*
         \[Epsilon][1]^8*\[Epsilon][2] + 8353200*Qrrc[0, 1]*\[Epsilon][1]^8*
         \[Epsilon][2] - 76734720*Qrrc[0, 0]^4*\[Epsilon][1]^9*
         \[Epsilon][2] - 15986400*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^9*
         \[Epsilon][2] + 113692785*\[Epsilon][1]^4*\[Epsilon][2]^2 + 
        538404780*Qrrc[0, 0]*\[Epsilon][1]^5*\[Epsilon][2]^2 - 
        410491200*Qrrc[0, 0]^2*\[Epsilon][1]^6*\[Epsilon][2]^2 - 
        114212160*Qrrc[0, 0]^3*\[Epsilon][1]^7*\[Epsilon][2]^2 - 
        20935800*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^2 + 
        1422696960*Qrrc[0, 0]^4*\[Epsilon][1]^8*\[Epsilon][2]^2 + 
        296395200*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^8*\[Epsilon][2]^2 - 
        368326656*Qrrc[0, 0]^5*\[Epsilon][1]^9*\[Epsilon][2]^2 - 
        191836800*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^9*\[Epsilon][2]^2 - 
        2797620*Qrrc[0, 2]*\[Epsilon][1]^9*\[Epsilon][2]^2 + 
        178279200*\[Epsilon][1]^3*\[Epsilon][2]^3 + 1854295260*Qrrc[0, 0]*
         \[Epsilon][1]^4*\[Epsilon][2]^3 + 2726335440*Qrrc[0, 0]^2*
         \[Epsilon][1]^5*\[Epsilon][2]^3 - 2321981760*Qrrc[0, 0]^3*
         \[Epsilon][1]^6*\[Epsilon][2]^3 - 495479800*Qrrc[0, 1]*
         \[Epsilon][1]^6*\[Epsilon][2]^3 + 3425126400*Qrrc[0, 0]^4*
         \[Epsilon][1]^7*\[Epsilon][2]^3 + 713568000*Qrrc[0, 0]*Qrrc[0, 1]*
         \[Epsilon][1]^7*\[Epsilon][2]^3 + 2309492736*Qrrc[0, 0]^5*
         \[Epsilon][1]^8*\[Epsilon][2]^3 + 1202860800*Qrrc[0, 0]^2*Qrrc[0, 1]*
         \[Epsilon][1]^8*\[Epsilon][2]^3 + 17541720*Qrrc[0, 2]*
         \[Epsilon][1]^8*\[Epsilon][2]^3 - 736653312*Qrrc[0, 0]^6*
         \[Epsilon][1]^9*\[Epsilon][2]^3 - 767347200*Qrrc[0, 0]^3*Qrrc[0, 1]*
         \[Epsilon][1]^9*\[Epsilon][2]^3 + 159864000*Qrrc[0, 1]^2*
         \[Epsilon][1]^9*\[Epsilon][2]^3 - 33571440*Qrrc[0, 0]*Qrrc[0, 2]*
         \[Epsilon][1]^9*\[Epsilon][2]^3 + 113692785*\[Epsilon][1]^2*
         \[Epsilon][2]^4 + 1854295260*Qrrc[0, 0]*\[Epsilon][1]^3*
         \[Epsilon][2]^4 + 6036888960*Qrrc[0, 0]^2*\[Epsilon][1]^4*
         \[Epsilon][2]^4 - 4408948800*Qrrc[0, 0]^3*\[Epsilon][1]^5*
         \[Epsilon][2]^4 - 1311298200*Qrrc[0, 1]*\[Epsilon][1]^5*
         \[Epsilon][2]^4 - 331522560*Qrrc[0, 0]^4*\[Epsilon][1]^6*
         \[Epsilon][2]^4 - 69067200*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^6*
         \[Epsilon][2]^4 + 16704036864*Qrrc[0, 0]^5*\[Epsilon][1]^7*
         \[Epsilon][2]^4 + 8700019200*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^7*
         \[Epsilon][2]^4 + 126875280*Qrrc[0, 2]*\[Epsilon][1]^7*
         \[Epsilon][2]^4 - 4419919872*Qrrc[0, 0]^6*\[Epsilon][1]^8*
         \[Epsilon][2]^4 - 4604083200*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^8*
         \[Epsilon][2]^4 + 959184000*Qrrc[0, 1]^2*\[Epsilon][1]^8*
         \[Epsilon][2]^4 - 201428640*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^8*
         \[Epsilon][2]^4 + 27025200*\[Epsilon][1]*\[Epsilon][2]^5 + 
        538404780*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2]^5 + 
        2726335440*Qrrc[0, 0]^2*\[Epsilon][1]^3*\[Epsilon][2]^5 - 
        4408948800*Qrrc[0, 0]^3*\[Epsilon][1]^4*\[Epsilon][2]^5 - 
        1311298200*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^5 - 
        4514434560*Qrrc[0, 0]^4*\[Epsilon][1]^5*\[Epsilon][2]^5 - 
        940507200*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^5 + 
        35986360320*Qrrc[0, 0]^5*\[Epsilon][1]^6*\[Epsilon][2]^5 + 
        18742896000*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^5 + 
        273333900*Qrrc[0, 2]*\[Epsilon][1]^6*\[Epsilon][2]^5 - 
        11049799680*Qrrc[0, 0]^6*\[Epsilon][1]^7*\[Epsilon][2]^5 - 
        11510208000*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^5 + 
        2397960000*Qrrc[0, 1]^2*\[Epsilon][1]^7*\[Epsilon][2]^5 - 
        503571600*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^7*\[Epsilon][2]^5 + 
        1165215*\[Epsilon][2]^6 - 20616140*Qrrc[0, 0]*\[Epsilon][1]*
         \[Epsilon][2]^6 - 410491200*Qrrc[0, 0]^2*\[Epsilon][1]^2*
         \[Epsilon][2]^6 - 2321981760*Qrrc[0, 0]^3*\[Epsilon][1]^3*
         \[Epsilon][2]^6 - 495479800*Qrrc[0, 1]*\[Epsilon][1]^3*
         \[Epsilon][2]^6 - 331522560*Qrrc[0, 0]^4*\[Epsilon][1]^4*
         \[Epsilon][2]^6 - 69067200*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^4*
         \[Epsilon][2]^6 + 35986360320*Qrrc[0, 0]^5*\[Epsilon][1]^5*
         \[Epsilon][2]^6 + 18742896000*Qrrc[0, 0]^2*Qrrc[0, 1]*
         \[Epsilon][1]^5*\[Epsilon][2]^6 + 273333900*Qrrc[0, 2]*
         \[Epsilon][1]^5*\[Epsilon][2]^6 - 14733066240*Qrrc[0, 0]^6*
         \[Epsilon][1]^6*\[Epsilon][2]^6 - 15346944000*Qrrc[0, 0]^3*
         Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^6 + 3197280000*Qrrc[0, 1]^2*
         \[Epsilon][1]^6*\[Epsilon][2]^6 - 671428800*Qrrc[0, 0]*Qrrc[0, 2]*
         \[Epsilon][1]^6*\[Epsilon][2]^6 - 8135260*Qrrc[0, 0]*
         \[Epsilon][2]^7 - 98766480*Qrrc[0, 0]^2*\[Epsilon][1]*
         \[Epsilon][2]^7 - 114212160*Qrrc[0, 0]^3*\[Epsilon][1]^2*
         \[Epsilon][2]^7 - 20935800*Qrrc[0, 1]*\[Epsilon][1]^2*
         \[Epsilon][2]^7 + 3425126400*Qrrc[0, 0]^4*\[Epsilon][1]^3*
         \[Epsilon][2]^7 + 713568000*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^3*
         \[Epsilon][2]^7 + 16704036864*Qrrc[0, 0]^5*\[Epsilon][1]^4*
         \[Epsilon][2]^7 + 8700019200*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^4*
         \[Epsilon][2]^7 + 126875280*Qrrc[0, 2]*\[Epsilon][1]^4*
         \[Epsilon][2]^7 - 11049799680*Qrrc[0, 0]^6*\[Epsilon][1]^5*
         \[Epsilon][2]^7 - 11510208000*Qrrc[0, 0]^3*Qrrc[0, 1]*
         \[Epsilon][1]^5*\[Epsilon][2]^7 + 2397960000*Qrrc[0, 1]^2*
         \[Epsilon][1]^5*\[Epsilon][2]^7 - 503571600*Qrrc[0, 0]*Qrrc[0, 2]*
         \[Epsilon][1]^5*\[Epsilon][2]^7 + 19615680*Qrrc[0, 0]^2*
         \[Epsilon][2]^8 + 275483520*Qrrc[0, 0]^3*\[Epsilon][1]*
         \[Epsilon][2]^8 + 8353200*Qrrc[0, 1]*\[Epsilon][1]*\[Epsilon][2]^8 + 
        1422696960*Qrrc[0, 0]^4*\[Epsilon][1]^2*\[Epsilon][2]^8 + 
        296395200*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^2*\[Epsilon][2]^8 + 
        2309492736*Qrrc[0, 0]^5*\[Epsilon][1]^3*\[Epsilon][2]^8 + 
        1202860800*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^8 + 
        17541720*Qrrc[0, 2]*\[Epsilon][1]^3*\[Epsilon][2]^8 - 
        4419919872*Qrrc[0, 0]^6*\[Epsilon][1]^4*\[Epsilon][2]^8 - 
        4604083200*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^8 + 
        959184000*Qrrc[0, 1]^2*\[Epsilon][1]^4*\[Epsilon][2]^8 - 
        201428640*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^4*\[Epsilon][2]^8 - 
        6394560*Qrrc[0, 0]^3*\[Epsilon][2]^9 - 1332200*Qrrc[0, 1]*
         \[Epsilon][2]^9 - 76734720*Qrrc[0, 0]^4*\[Epsilon][1]*
         \[Epsilon][2]^9 - 15986400*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]*
         \[Epsilon][2]^9 - 368326656*Qrrc[0, 0]^5*\[Epsilon][1]^2*
         \[Epsilon][2]^9 - 191836800*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^2*
         \[Epsilon][2]^9 - 2797620*Qrrc[0, 2]*\[Epsilon][1]^2*
         \[Epsilon][2]^9 - 736653312*Qrrc[0, 0]^6*\[Epsilon][1]^3*
         \[Epsilon][2]^9 - 767347200*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^3*
         \[Epsilon][2]^9 + 159864000*Qrrc[0, 1]^2*\[Epsilon][1]^3*
         \[Epsilon][2]^9 - 33571440*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^3*
         \[Epsilon][2]^9)/(1269789696000*\[Epsilon][1]*\[Epsilon][2]), 
      (-1131150*\[Epsilon][1]^6 + 8764735*Qrrc[0, 0]*\[Epsilon][1]^7 - 
        14217840*Qrrc[0, 0]^2*\[Epsilon][1]^8 + 123360*Qrrc[0, 0]^3*
         \[Epsilon][1]^9 + 25700*Qrrc[0, 1]*\[Epsilon][1]^9 - 
        16380000*\[Epsilon][1]^5*\[Epsilon][2] + 24916205*Qrrc[0, 0]*
         \[Epsilon][1]^6*\[Epsilon][2] + 139295280*Qrrc[0, 0]^2*
         \[Epsilon][1]^7*\[Epsilon][2] - 226745280*Qrrc[0, 0]^3*
         \[Epsilon][1]^8*\[Epsilon][2] - 11694000*Qrrc[0, 1]*\[Epsilon][1]^8*
         \[Epsilon][2] + 1480320*Qrrc[0, 0]^4*\[Epsilon][1]^9*\[Epsilon][2] + 
        308400*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^9*\[Epsilon][2] - 
        67448850*\[Epsilon][1]^4*\[Epsilon][2]^2 - 274352025*Qrrc[0, 0]*
         \[Epsilon][1]^5*\[Epsilon][2]^2 + 549241860*Qrrc[0, 0]^2*
         \[Epsilon][1]^6*\[Epsilon][2]^2 + 547745760*Qrrc[0, 0]^3*
         \[Epsilon][1]^7*\[Epsilon][2]^2 + 28817550*Qrrc[0, 1]*
         \[Epsilon][1]^7*\[Epsilon][2]^2 - 1356030720*Qrrc[0, 0]^4*
         \[Epsilon][1]^8*\[Epsilon][2]^2 - 282506400*Qrrc[0, 0]*Qrrc[0, 1]*
         \[Epsilon][1]^8*\[Epsilon][2]^2 + 7105536*Qrrc[0, 0]^5*
         \[Epsilon][1]^9*\[Epsilon][2]^2 + 3700800*Qrrc[0, 0]^2*Qrrc[0, 1]*
         \[Epsilon][1]^9*\[Epsilon][2]^2 + 53970*Qrrc[0, 2]*\[Epsilon][1]^9*
         \[Epsilon][2]^2 - 106344000*\[Epsilon][1]^3*\[Epsilon][2]^3 - 
        1031594355*Qrrc[0, 0]*\[Epsilon][1]^4*\[Epsilon][2]^3 - 
        941536080*Qrrc[0, 0]^2*\[Epsilon][1]^5*\[Epsilon][2]^3 + 
        3633317040*Qrrc[0, 0]^3*\[Epsilon][1]^6*\[Epsilon][2]^3 + 
        538536550*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^3 - 
        1750910400*Qrrc[0, 0]^4*\[Epsilon][1]^7*\[Epsilon][2]^3 - 
        364773000*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^7*\[Epsilon][2]^3 - 
        3233157120*Qrrc[0, 0]^5*\[Epsilon][1]^8*\[Epsilon][2]^3 - 
        1683936000*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^8*\[Epsilon][2]^3 - 
        24557400*Qrrc[0, 2]*\[Epsilon][1]^8*\[Epsilon][2]^3 + 
        14211072*Qrrc[0, 0]^6*\[Epsilon][1]^9*\[Epsilon][2]^3 + 
        14803200*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^9*\[Epsilon][2]^3 - 
        3084000*Qrrc[0, 1]^2*\[Epsilon][1]^9*\[Epsilon][2]^3 + 
        647640*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^9*\[Epsilon][2]^3 - 
        67448850*\[Epsilon][1]^2*\[Epsilon][2]^4 - 1031594355*Qrrc[0, 0]*
         \[Epsilon][1]^3*\[Epsilon][2]^4 - 2674529640*Qrrc[0, 0]^2*
         \[Epsilon][1]^4*\[Epsilon][2]^4 + 6572895120*Qrrc[0, 0]^3*
         \[Epsilon][1]^5*\[Epsilon][2]^4 + 1389432600*Qrrc[0, 1]*
         \[Epsilon][1]^5*\[Epsilon][2]^4 + 5224320000*Qrrc[0, 0]^4*
         \[Epsilon][1]^6*\[Epsilon][2]^4 + 1088400000*Qrrc[0, 0]*Qrrc[0, 1]*
         \[Epsilon][1]^6*\[Epsilon][2]^4 - 16265263104*Qrrc[0, 0]^5*
         \[Epsilon][1]^7*\[Epsilon][2]^4 - 8471491200*Qrrc[0, 0]^2*Qrrc[0, 1]*
         \[Epsilon][1]^7*\[Epsilon][2]^4 - 123542580*Qrrc[0, 2]*
         \[Epsilon][1]^7*\[Epsilon][2]^4 + 85266432*Qrrc[0, 0]^6*
         \[Epsilon][1]^8*\[Epsilon][2]^4 + 88819200*Qrrc[0, 0]^3*Qrrc[0, 1]*
         \[Epsilon][1]^8*\[Epsilon][2]^4 - 18504000*Qrrc[0, 1]^2*
         \[Epsilon][1]^8*\[Epsilon][2]^4 + 3885840*Qrrc[0, 0]*Qrrc[0, 2]*
         \[Epsilon][1]^8*\[Epsilon][2]^4 - 16380000*\[Epsilon][1]*
         \[Epsilon][2]^5 - 274352025*Qrrc[0, 0]*\[Epsilon][1]^2*
         \[Epsilon][2]^5 - 941536080*Qrrc[0, 0]^2*\[Epsilon][1]^3*
         \[Epsilon][2]^5 + 6572895120*Qrrc[0, 0]^3*\[Epsilon][1]^4*
         \[Epsilon][2]^5 + 1389432600*Qrrc[0, 1]*\[Epsilon][1]^4*
         \[Epsilon][2]^5 + 11235438720*Qrrc[0, 0]^4*\[Epsilon][1]^5*
         \[Epsilon][2]^5 + 2340716400*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^5*
         \[Epsilon][2]^5 - 32580264960*Qrrc[0, 0]^5*\[Epsilon][1]^6*
         \[Epsilon][2]^5 - 16968888000*Qrrc[0, 0]^2*Qrrc[0, 1]*
         \[Epsilon][1]^6*\[Epsilon][2]^5 - 247462950*Qrrc[0, 2]*
         \[Epsilon][1]^6*\[Epsilon][2]^5 + 213166080*Qrrc[0, 0]^6*
         \[Epsilon][1]^7*\[Epsilon][2]^5 + 222048000*Qrrc[0, 0]^3*Qrrc[0, 1]*
         \[Epsilon][1]^7*\[Epsilon][2]^5 - 46260000*Qrrc[0, 1]^2*
         \[Epsilon][1]^7*\[Epsilon][2]^5 + 9714600*Qrrc[0, 0]*Qrrc[0, 2]*
         \[Epsilon][1]^7*\[Epsilon][2]^5 - 1131150*\[Epsilon][2]^6 + 
        24916205*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^6 + 
        549241860*Qrrc[0, 0]^2*\[Epsilon][1]^2*\[Epsilon][2]^6 + 
        3633317040*Qrrc[0, 0]^3*\[Epsilon][1]^3*\[Epsilon][2]^6 + 
        538536550*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^6 + 
        5224320000*Qrrc[0, 0]^4*\[Epsilon][1]^4*\[Epsilon][2]^6 + 
        1088400000*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^4*\[Epsilon][2]^6 - 
        32580264960*Qrrc[0, 0]^5*\[Epsilon][1]^5*\[Epsilon][2]^6 - 
        16968888000*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^5*\[Epsilon][2]^6 - 
        247462950*Qrrc[0, 2]*\[Epsilon][1]^5*\[Epsilon][2]^6 + 
        284221440*Qrrc[0, 0]^6*\[Epsilon][1]^6*\[Epsilon][2]^6 + 
        296064000*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^6*\[Epsilon][2]^6 - 
        61680000*Qrrc[0, 1]^2*\[Epsilon][1]^6*\[Epsilon][2]^6 + 
        12952800*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^6*\[Epsilon][2]^6 + 
        8764735*Qrrc[0, 0]*\[Epsilon][2]^7 + 139295280*Qrrc[0, 0]^2*
         \[Epsilon][1]*\[Epsilon][2]^7 + 547745760*Qrrc[0, 0]^3*
         \[Epsilon][1]^2*\[Epsilon][2]^7 + 28817550*Qrrc[0, 1]*
         \[Epsilon][1]^2*\[Epsilon][2]^7 - 1750910400*Qrrc[0, 0]^4*
         \[Epsilon][1]^3*\[Epsilon][2]^7 - 364773000*Qrrc[0, 0]*Qrrc[0, 1]*
         \[Epsilon][1]^3*\[Epsilon][2]^7 - 16265263104*Qrrc[0, 0]^5*
         \[Epsilon][1]^4*\[Epsilon][2]^7 - 8471491200*Qrrc[0, 0]^2*Qrrc[0, 1]*
         \[Epsilon][1]^4*\[Epsilon][2]^7 - 123542580*Qrrc[0, 2]*
         \[Epsilon][1]^4*\[Epsilon][2]^7 + 213166080*Qrrc[0, 0]^6*
         \[Epsilon][1]^5*\[Epsilon][2]^7 + 222048000*Qrrc[0, 0]^3*Qrrc[0, 1]*
         \[Epsilon][1]^5*\[Epsilon][2]^7 - 46260000*Qrrc[0, 1]^2*
         \[Epsilon][1]^5*\[Epsilon][2]^7 + 9714600*Qrrc[0, 0]*Qrrc[0, 2]*
         \[Epsilon][1]^5*\[Epsilon][2]^7 - 14217840*Qrrc[0, 0]^2*
         \[Epsilon][2]^8 - 226745280*Qrrc[0, 0]^3*\[Epsilon][1]*
         \[Epsilon][2]^8 - 11694000*Qrrc[0, 1]*\[Epsilon][1]*
         \[Epsilon][2]^8 - 1356030720*Qrrc[0, 0]^4*\[Epsilon][1]^2*
         \[Epsilon][2]^8 - 282506400*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]^2*
         \[Epsilon][2]^8 - 3233157120*Qrrc[0, 0]^5*\[Epsilon][1]^3*
         \[Epsilon][2]^8 - 1683936000*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^3*
         \[Epsilon][2]^8 - 24557400*Qrrc[0, 2]*\[Epsilon][1]^3*
         \[Epsilon][2]^8 + 85266432*Qrrc[0, 0]^6*\[Epsilon][1]^4*
         \[Epsilon][2]^8 + 88819200*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^4*
         \[Epsilon][2]^8 - 18504000*Qrrc[0, 1]^2*\[Epsilon][1]^4*
         \[Epsilon][2]^8 + 3885840*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^4*
         \[Epsilon][2]^8 + 123360*Qrrc[0, 0]^3*\[Epsilon][2]^9 + 
        25700*Qrrc[0, 1]*\[Epsilon][2]^9 + 1480320*Qrrc[0, 0]^4*\[Epsilon][1]*
         \[Epsilon][2]^9 + 308400*Qrrc[0, 0]*Qrrc[0, 1]*\[Epsilon][1]*
         \[Epsilon][2]^9 + 7105536*Qrrc[0, 0]^5*\[Epsilon][1]^2*
         \[Epsilon][2]^9 + 3700800*Qrrc[0, 0]^2*Qrrc[0, 1]*\[Epsilon][1]^2*
         \[Epsilon][2]^9 + 53970*Qrrc[0, 2]*\[Epsilon][1]^2*\[Epsilon][2]^9 + 
        14211072*Qrrc[0, 0]^6*\[Epsilon][1]^3*\[Epsilon][2]^9 + 
        14803200*Qrrc[0, 0]^3*Qrrc[0, 1]*\[Epsilon][1]^3*\[Epsilon][2]^9 - 
        3084000*Qrrc[0, 1]^2*\[Epsilon][1]^3*\[Epsilon][2]^9 + 
        647640*Qrrc[0, 0]*Qrrc[0, 2]*\[Epsilon][1]^3*\[Epsilon][2]^9)/
       (2857026816000*\[Epsilon][1]*\[Epsilon][2])}, -4, 2, 1]
 
invcf = SeriesData[tcv, 0, {1, -11/18, 145/486, -6733/52488, 120127/2361960, 
      -2431777/127545840, 11024687/1607077584, -4151904751/1735643790720, 
      114405084629/140587147048320, -30899177349251/113875589109139200}, 1, 
     11, 1]
 
q1sol = {{Qrrc[0, 1] -> (-(Qrrc[0, 0]*\[Epsilon][1]^3) - 
         Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] - 12*Qrrc[0, 0]^2*
          \[Epsilon][1]^3*\[Epsilon][2] - Qrrc[0, 0]*\[Epsilon][1]*
          \[Epsilon][2]^2 - 12*Qrrc[0, 0]^2*\[Epsilon][1]^2*\[Epsilon][2]^2 - 
         48*Qrrc[0, 0]^3*\[Epsilon][1]^3*\[Epsilon][2]^2 - 
         Qrrc[0, 0]*\[Epsilon][2]^3 - 12*Qrrc[0, 0]^2*\[Epsilon][1]*
          \[Epsilon][2]^3 - 48*Qrrc[0, 0]^3*\[Epsilon][1]^2*\[Epsilon][2]^3)/
        (10*\[Epsilon][1]^2*\[Epsilon][2]^2*(\[Epsilon][1] + \[Epsilon][2]))}}
 
q2sol[-3] = 
    {{Qrrc[0, 2] -> (-1/30240*(Qrrc[0, 0]*(\[Epsilon][1] + \[Epsilon][2])^3*
           (\[Epsilon][1]^2 - 5*\[Epsilon][1]*\[Epsilon][2] + 
            12*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + \[Epsilon][2]^2 + 
            12*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2)) - 
         (Qrrc[0, 0]^2*\[Epsilon][1]*(\[Epsilon][1] + \[Epsilon][2])^3*
           (\[Epsilon][1]^2 - 5*\[Epsilon][1]*\[Epsilon][2] + 
            12*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + \[Epsilon][2]^2 + 
            12*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/1008 - 
         (Qrrc[0, 0]^3*\[Epsilon][1]^2*(\[Epsilon][1] + \[Epsilon][2])^3*
           (\[Epsilon][1]^2 - 5*\[Epsilon][1]*\[Epsilon][2] + 
            12*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + \[Epsilon][2]^2 + 
            12*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/105 - 
         (Qrrc[0, 0]^4*\[Epsilon][1]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (\[Epsilon][1]^2 - 5*\[Epsilon][1]*\[Epsilon][2] + 
            12*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + \[Epsilon][2]^2 + 
            12*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/35 - 
         (Qrrc[0, 0]*\[Epsilon][1]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (\[Epsilon][1]^2 - 5*\[Epsilon][1]*\[Epsilon][2] + 
            12*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + \[Epsilon][2]^2 + 
            12*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/
          (60480*\[Epsilon][2]^3) - (Qrrc[0, 0]*\[Epsilon][1]^2*
           (\[Epsilon][1] + \[Epsilon][2])^3*(\[Epsilon][1]^2 - 
            5*\[Epsilon][1]*\[Epsilon][2] + 12*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + \[Epsilon][2]^2 + 12*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(30240*\[Epsilon][2]^2) - 
         (Qrrc[0, 0]^2*\[Epsilon][1]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (\[Epsilon][1]^2 - 5*\[Epsilon][1]*\[Epsilon][2] + 
            12*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + \[Epsilon][2]^2 + 
            12*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/
          (2520*\[Epsilon][2]^2) - (Qrrc[0, 0]*\[Epsilon][1]*
           (\[Epsilon][1] + \[Epsilon][2])^3*(\[Epsilon][1]^2 - 
            5*\[Epsilon][1]*\[Epsilon][2] + 12*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + \[Epsilon][2]^2 + 12*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(30240*\[Epsilon][2]) - 
         (Qrrc[0, 0]^2*\[Epsilon][1]^2*(\[Epsilon][1] + \[Epsilon][2])^3*
           (\[Epsilon][1]^2 - 5*\[Epsilon][1]*\[Epsilon][2] + 
            12*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + \[Epsilon][2]^2 + 
            12*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/
          (1260*\[Epsilon][2]) - (Qrrc[0, 0]^3*\[Epsilon][1]^3*
           (\[Epsilon][1] + \[Epsilon][2])^3*(\[Epsilon][1]^2 - 
            5*\[Epsilon][1]*\[Epsilon][2] + 12*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + \[Epsilon][2]^2 + 12*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(210*\[Epsilon][2]) - 
         (Qrrc[0, 0]^2*\[Epsilon][2]*(\[Epsilon][1] + \[Epsilon][2])^3*
           (\[Epsilon][1]^2 - 5*\[Epsilon][1]*\[Epsilon][2] + 
            12*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + \[Epsilon][2]^2 + 
            12*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/1008 - 
         (Qrrc[0, 0]*\[Epsilon][2]*(\[Epsilon][1] + \[Epsilon][2])^3*
           (\[Epsilon][1]^2 - 5*\[Epsilon][1]*\[Epsilon][2] + 
            12*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + \[Epsilon][2]^2 + 
            12*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/
          (30240*\[Epsilon][1]) - (Qrrc[0, 0]^3*\[Epsilon][1]*\[Epsilon][2]*
           (\[Epsilon][1] + \[Epsilon][2])^3*(\[Epsilon][1]^2 - 
            5*\[Epsilon][1]*\[Epsilon][2] + 12*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + \[Epsilon][2]^2 + 12*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/84 - (2*Qrrc[0, 0]^4*\[Epsilon][1]^2*
           \[Epsilon][2]*(\[Epsilon][1] + \[Epsilon][2])^3*
           (\[Epsilon][1]^2 - 5*\[Epsilon][1]*\[Epsilon][2] + 
            12*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + \[Epsilon][2]^2 + 
            12*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/35 - 
         (12*Qrrc[0, 0]^5*\[Epsilon][1]^3*\[Epsilon][2]*
           (\[Epsilon][1] + \[Epsilon][2])^3*(\[Epsilon][1]^2 - 
            5*\[Epsilon][1]*\[Epsilon][2] + 12*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + \[Epsilon][2]^2 + 12*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/175 - (Qrrc[0, 0]^3*\[Epsilon][2]^2*
           (\[Epsilon][1] + \[Epsilon][2])^3*(\[Epsilon][1]^2 - 
            5*\[Epsilon][1]*\[Epsilon][2] + 12*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + \[Epsilon][2]^2 + 12*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/105 - (Qrrc[0, 0]*\[Epsilon][2]^2*
           (\[Epsilon][1] + \[Epsilon][2])^3*(\[Epsilon][1]^2 - 
            5*\[Epsilon][1]*\[Epsilon][2] + 12*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + \[Epsilon][2]^2 + 12*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(30240*\[Epsilon][1]^2) - 
         (Qrrc[0, 0]^2*\[Epsilon][2]^2*(\[Epsilon][1] + \[Epsilon][2])^3*
           (\[Epsilon][1]^2 - 5*\[Epsilon][1]*\[Epsilon][2] + 
            12*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + \[Epsilon][2]^2 + 
            12*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/
          (1260*\[Epsilon][1]) - (2*Qrrc[0, 0]^4*\[Epsilon][1]*
           \[Epsilon][2]^2*(\[Epsilon][1] + \[Epsilon][2])^3*
           (\[Epsilon][1]^2 - 5*\[Epsilon][1]*\[Epsilon][2] + 
            12*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + \[Epsilon][2]^2 + 
            12*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/35 - 
         (24*Qrrc[0, 0]^5*\[Epsilon][1]^2*\[Epsilon][2]^2*
           (\[Epsilon][1] + \[Epsilon][2])^3*(\[Epsilon][1]^2 - 
            5*\[Epsilon][1]*\[Epsilon][2] + 12*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + \[Epsilon][2]^2 + 12*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/175 - (Qrrc[0, 0]^4*\[Epsilon][2]^3*
           (\[Epsilon][1] + \[Epsilon][2])^3*(\[Epsilon][1]^2 - 
            5*\[Epsilon][1]*\[Epsilon][2] + 12*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + \[Epsilon][2]^2 + 12*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/35 - (Qrrc[0, 0]*\[Epsilon][2]^3*
           (\[Epsilon][1] + \[Epsilon][2])^3*(\[Epsilon][1]^2 - 
            5*\[Epsilon][1]*\[Epsilon][2] + 12*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + \[Epsilon][2]^2 + 12*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(60480*\[Epsilon][1]^3) - 
         (Qrrc[0, 0]^2*\[Epsilon][2]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (\[Epsilon][1]^2 - 5*\[Epsilon][1]*\[Epsilon][2] + 
            12*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + \[Epsilon][2]^2 + 
            12*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/
          (2520*\[Epsilon][1]^2) - (Qrrc[0, 0]^3*\[Epsilon][2]^3*
           (\[Epsilon][1] + \[Epsilon][2])^3*(\[Epsilon][1]^2 - 
            5*\[Epsilon][1]*\[Epsilon][2] + 12*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + \[Epsilon][2]^2 + 12*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(210*\[Epsilon][1]) - 
         (12*Qrrc[0, 0]^5*\[Epsilon][1]*\[Epsilon][2]^3*
           (\[Epsilon][1] + \[Epsilon][2])^3*(\[Epsilon][1]^2 - 
            5*\[Epsilon][1]*\[Epsilon][2] + 12*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + \[Epsilon][2]^2 + 12*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/175)/(-1/2880*(\[Epsilon][1]^3*\[Epsilon][2]*
           (\[Epsilon][1] + \[Epsilon][2])^3*(\[Epsilon][1]^2 - 
            5*\[Epsilon][1]*\[Epsilon][2] + 12*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + \[Epsilon][2]^2 + 12*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2)) - (\[Epsilon][1]^2*\[Epsilon][2]^2*
           (\[Epsilon][1] + \[Epsilon][2])^3*(\[Epsilon][1]^2 - 
            5*\[Epsilon][1]*\[Epsilon][2] + 12*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + \[Epsilon][2]^2 + 12*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/1440 - (\[Epsilon][1]*\[Epsilon][2]^3*
           (\[Epsilon][1] + \[Epsilon][2])^3*(\[Epsilon][1]^2 - 
            5*\[Epsilon][1]*\[Epsilon][2] + 12*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + \[Epsilon][2]^2 + 12*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/2880)}}
 
q2sol[-2] = {{Qrrc[0, 2] -> ((Qrrc[0, 0]*(\[Epsilon][1] + \[Epsilon][2])^3*
           (11*\[Epsilon][1]^2 - 13*\[Epsilon][1]*\[Epsilon][2] + 
            132*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 
            11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/3265920 + (Qrrc[0, 0]^2*\[Epsilon][1]*
           (\[Epsilon][1] + \[Epsilon][2])^3*(11*\[Epsilon][1]^2 - 
            13*\[Epsilon][1]*\[Epsilon][2] + 132*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*
             \[Epsilon][1]*\[Epsilon][2]^2))/108864 + 
         (Qrrc[0, 0]^3*\[Epsilon][1]^2*(\[Epsilon][1] + \[Epsilon][2])^3*
           (11*\[Epsilon][1]^2 - 13*\[Epsilon][1]*\[Epsilon][2] + 
            132*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 
            11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/11340 + (Qrrc[0, 0]^4*\[Epsilon][1]^3*
           (\[Epsilon][1] + \[Epsilon][2])^3*(11*\[Epsilon][1]^2 - 
            13*\[Epsilon][1]*\[Epsilon][2] + 132*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*
             \[Epsilon][1]*\[Epsilon][2]^2))/3780 + 
         (Qrrc[0, 0]*\[Epsilon][1]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (11*\[Epsilon][1]^2 - 13*\[Epsilon][1]*\[Epsilon][2] + 
            132*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 
            11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(6531840*\[Epsilon][2]^3) + 
         (Qrrc[0, 0]*\[Epsilon][1]^2*(\[Epsilon][1] + \[Epsilon][2])^3*
           (11*\[Epsilon][1]^2 - 13*\[Epsilon][1]*\[Epsilon][2] + 
            132*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 
            11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(3265920*\[Epsilon][2]^2) + 
         (Qrrc[0, 0]^2*\[Epsilon][1]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (11*\[Epsilon][1]^2 - 13*\[Epsilon][1]*\[Epsilon][2] + 
            132*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 
            11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(272160*\[Epsilon][2]^2) + 
         (Qrrc[0, 0]*\[Epsilon][1]*(\[Epsilon][1] + \[Epsilon][2])^3*
           (11*\[Epsilon][1]^2 - 13*\[Epsilon][1]*\[Epsilon][2] + 
            132*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 
            11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(3265920*\[Epsilon][2]) + 
         (Qrrc[0, 0]^2*\[Epsilon][1]^2*(\[Epsilon][1] + \[Epsilon][2])^3*
           (11*\[Epsilon][1]^2 - 13*\[Epsilon][1]*\[Epsilon][2] + 
            132*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 
            11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(136080*\[Epsilon][2]) + 
         (Qrrc[0, 0]^3*\[Epsilon][1]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (11*\[Epsilon][1]^2 - 13*\[Epsilon][1]*\[Epsilon][2] + 
            132*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 
            11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(22680*\[Epsilon][2]) + 
         (Qrrc[0, 0]^2*\[Epsilon][2]*(\[Epsilon][1] + \[Epsilon][2])^3*
           (11*\[Epsilon][1]^2 - 13*\[Epsilon][1]*\[Epsilon][2] + 
            132*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 
            11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/108864 + (Qrrc[0, 0]*\[Epsilon][2]*
           (\[Epsilon][1] + \[Epsilon][2])^3*(11*\[Epsilon][1]^2 - 
            13*\[Epsilon][1]*\[Epsilon][2] + 132*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*
             \[Epsilon][1]*\[Epsilon][2]^2))/(3265920*\[Epsilon][1]) + 
         (Qrrc[0, 0]^3*\[Epsilon][1]*\[Epsilon][2]*
           (\[Epsilon][1] + \[Epsilon][2])^3*(11*\[Epsilon][1]^2 - 
            13*\[Epsilon][1]*\[Epsilon][2] + 132*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*
             \[Epsilon][1]*\[Epsilon][2]^2))/9072 + 
         (Qrrc[0, 0]^4*\[Epsilon][1]^2*\[Epsilon][2]*
           (\[Epsilon][1] + \[Epsilon][2])^3*(11*\[Epsilon][1]^2 - 
            13*\[Epsilon][1]*\[Epsilon][2] + 132*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*
             \[Epsilon][1]*\[Epsilon][2]^2))/1890 + 
         (Qrrc[0, 0]^5*\[Epsilon][1]^3*\[Epsilon][2]*
           (\[Epsilon][1] + \[Epsilon][2])^3*(11*\[Epsilon][1]^2 - 
            13*\[Epsilon][1]*\[Epsilon][2] + 132*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*
             \[Epsilon][1]*\[Epsilon][2]^2))/1575 + 
         (Qrrc[0, 0]^3*\[Epsilon][2]^2*(\[Epsilon][1] + \[Epsilon][2])^3*
           (11*\[Epsilon][1]^2 - 13*\[Epsilon][1]*\[Epsilon][2] + 
            132*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 
            11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/11340 + (Qrrc[0, 0]*\[Epsilon][2]^2*
           (\[Epsilon][1] + \[Epsilon][2])^3*(11*\[Epsilon][1]^2 - 
            13*\[Epsilon][1]*\[Epsilon][2] + 132*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*
             \[Epsilon][1]*\[Epsilon][2]^2))/(3265920*\[Epsilon][1]^2) + 
         (Qrrc[0, 0]^2*\[Epsilon][2]^2*(\[Epsilon][1] + \[Epsilon][2])^3*
           (11*\[Epsilon][1]^2 - 13*\[Epsilon][1]*\[Epsilon][2] + 
            132*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 
            11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(136080*\[Epsilon][1]) + 
         (Qrrc[0, 0]^4*\[Epsilon][1]*\[Epsilon][2]^2*
           (\[Epsilon][1] + \[Epsilon][2])^3*(11*\[Epsilon][1]^2 - 
            13*\[Epsilon][1]*\[Epsilon][2] + 132*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*
             \[Epsilon][1]*\[Epsilon][2]^2))/1890 + 
         (2*Qrrc[0, 0]^5*\[Epsilon][1]^2*\[Epsilon][2]^2*
           (\[Epsilon][1] + \[Epsilon][2])^3*(11*\[Epsilon][1]^2 - 
            13*\[Epsilon][1]*\[Epsilon][2] + 132*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*
             \[Epsilon][1]*\[Epsilon][2]^2))/1575 + 
         (Qrrc[0, 0]^4*\[Epsilon][2]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (11*\[Epsilon][1]^2 - 13*\[Epsilon][1]*\[Epsilon][2] + 
            132*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 
            11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/3780 + (Qrrc[0, 0]*\[Epsilon][2]^3*
           (\[Epsilon][1] + \[Epsilon][2])^3*(11*\[Epsilon][1]^2 - 
            13*\[Epsilon][1]*\[Epsilon][2] + 132*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*
             \[Epsilon][1]*\[Epsilon][2]^2))/(6531840*\[Epsilon][1]^3) + 
         (Qrrc[0, 0]^2*\[Epsilon][2]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (11*\[Epsilon][1]^2 - 13*\[Epsilon][1]*\[Epsilon][2] + 
            132*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 
            11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(272160*\[Epsilon][1]^2) + 
         (Qrrc[0, 0]^3*\[Epsilon][2]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (11*\[Epsilon][1]^2 - 13*\[Epsilon][1]*\[Epsilon][2] + 
            132*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 
            11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(22680*\[Epsilon][1]) + 
         (Qrrc[0, 0]^5*\[Epsilon][1]*\[Epsilon][2]^3*
           (\[Epsilon][1] + \[Epsilon][2])^3*(11*\[Epsilon][1]^2 - 
            13*\[Epsilon][1]*\[Epsilon][2] + 132*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*
             \[Epsilon][1]*\[Epsilon][2]^2))/1575)/
        ((\[Epsilon][1]^3*\[Epsilon][2]*(\[Epsilon][1] + \[Epsilon][2])^3*
           (11*\[Epsilon][1]^2 - 13*\[Epsilon][1]*\[Epsilon][2] + 
            132*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 
            11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/311040 + (\[Epsilon][1]^2*\[Epsilon][2]^2*
           (\[Epsilon][1] + \[Epsilon][2])^3*(11*\[Epsilon][1]^2 - 
            13*\[Epsilon][1]*\[Epsilon][2] + 132*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*
             \[Epsilon][1]*\[Epsilon][2]^2))/155520 + 
         (\[Epsilon][1]*\[Epsilon][2]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (11*\[Epsilon][1]^2 - 13*\[Epsilon][1]*\[Epsilon][2] + 
            132*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 
            11*\[Epsilon][2]^2 + 132*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/311040)}}
 
q2sol[-1] = 
    {{Qrrc[0, 2] -> (-1/44089920*(Qrrc[0, 0]*(\[Epsilon][1] + \[Epsilon][2])^
            3*(7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2)) - 
         (Qrrc[0, 0]^2*\[Epsilon][1]*(\[Epsilon][1] + \[Epsilon][2])^3*
           (7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/1469664 - 
         (Qrrc[0, 0]^3*\[Epsilon][1]^2*(\[Epsilon][1] + \[Epsilon][2])^3*
           (7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/153090 - 
         (Qrrc[0, 0]^4*\[Epsilon][1]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/51030 - 
         (Qrrc[0, 0]*\[Epsilon][1]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/
          (88179840*\[Epsilon][2]^3) - (Qrrc[0, 0]*\[Epsilon][1]^2*
           (\[Epsilon][1] + \[Epsilon][2])^3*(7*\[Epsilon][1]^2 + 
            31*\[Epsilon][1]*\[Epsilon][2] + 84*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 7*\[Epsilon][2]^2 + 84*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(44089920*\[Epsilon][2]^2) - 
         (Qrrc[0, 0]^2*\[Epsilon][1]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/
          (3674160*\[Epsilon][2]^2) - (Qrrc[0, 0]*\[Epsilon][1]*
           (\[Epsilon][1] + \[Epsilon][2])^3*(7*\[Epsilon][1]^2 + 
            31*\[Epsilon][1]*\[Epsilon][2] + 84*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 7*\[Epsilon][2]^2 + 84*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(44089920*\[Epsilon][2]) - 
         (Qrrc[0, 0]^2*\[Epsilon][1]^2*(\[Epsilon][1] + \[Epsilon][2])^3*
           (7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/
          (1837080*\[Epsilon][2]) - (Qrrc[0, 0]^3*\[Epsilon][1]^3*
           (\[Epsilon][1] + \[Epsilon][2])^3*(7*\[Epsilon][1]^2 + 
            31*\[Epsilon][1]*\[Epsilon][2] + 84*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 7*\[Epsilon][2]^2 + 84*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(306180*\[Epsilon][2]) - 
         (Qrrc[0, 0]^2*\[Epsilon][2]*(\[Epsilon][1] + \[Epsilon][2])^3*
           (7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/1469664 - 
         (Qrrc[0, 0]*\[Epsilon][2]*(\[Epsilon][1] + \[Epsilon][2])^3*
           (7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/
          (44089920*\[Epsilon][1]) - (Qrrc[0, 0]^3*\[Epsilon][1]*
           \[Epsilon][2]*(\[Epsilon][1] + \[Epsilon][2])^3*
           (7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/122472 - 
         (Qrrc[0, 0]^4*\[Epsilon][1]^2*\[Epsilon][2]*
           (\[Epsilon][1] + \[Epsilon][2])^3*(7*\[Epsilon][1]^2 + 
            31*\[Epsilon][1]*\[Epsilon][2] + 84*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 7*\[Epsilon][2]^2 + 84*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/25515 - (2*Qrrc[0, 0]^5*\[Epsilon][1]^3*
           \[Epsilon][2]*(\[Epsilon][1] + \[Epsilon][2])^3*
           (7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/42525 - 
         (Qrrc[0, 0]^3*\[Epsilon][2]^2*(\[Epsilon][1] + \[Epsilon][2])^3*
           (7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/153090 - 
         (Qrrc[0, 0]*\[Epsilon][2]^2*(\[Epsilon][1] + \[Epsilon][2])^3*
           (7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/
          (44089920*\[Epsilon][1]^2) - (Qrrc[0, 0]^2*\[Epsilon][2]^2*
           (\[Epsilon][1] + \[Epsilon][2])^3*(7*\[Epsilon][1]^2 + 
            31*\[Epsilon][1]*\[Epsilon][2] + 84*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 7*\[Epsilon][2]^2 + 84*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(1837080*\[Epsilon][1]) - 
         (Qrrc[0, 0]^4*\[Epsilon][1]*\[Epsilon][2]^2*
           (\[Epsilon][1] + \[Epsilon][2])^3*(7*\[Epsilon][1]^2 + 
            31*\[Epsilon][1]*\[Epsilon][2] + 84*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 7*\[Epsilon][2]^2 + 84*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/25515 - (4*Qrrc[0, 0]^5*\[Epsilon][1]^2*
           \[Epsilon][2]^2*(\[Epsilon][1] + \[Epsilon][2])^3*
           (7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/42525 - 
         (Qrrc[0, 0]^4*\[Epsilon][2]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/51030 - 
         (Qrrc[0, 0]*\[Epsilon][2]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/
          (88179840*\[Epsilon][1]^3) - (Qrrc[0, 0]^2*\[Epsilon][2]^3*
           (\[Epsilon][1] + \[Epsilon][2])^3*(7*\[Epsilon][1]^2 + 
            31*\[Epsilon][1]*\[Epsilon][2] + 84*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 7*\[Epsilon][2]^2 + 84*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/(3674160*\[Epsilon][1]^2) - 
         (Qrrc[0, 0]^3*\[Epsilon][2]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/
          (306180*\[Epsilon][1]) - (2*Qrrc[0, 0]^5*\[Epsilon][1]*
           \[Epsilon][2]^3*(\[Epsilon][1] + \[Epsilon][2])^3*
           (7*\[Epsilon][1]^2 + 31*\[Epsilon][1]*\[Epsilon][2] + 
            84*Qrrc[0, 0]*\[Epsilon][1]^2*\[Epsilon][2] + 7*\[Epsilon][2]^2 + 
            84*Qrrc[0, 0]*\[Epsilon][1]*\[Epsilon][2]^2))/42525)/
        (-1/4199040*(\[Epsilon][1]^3*\[Epsilon][2]*
           (\[Epsilon][1] + \[Epsilon][2])^3*(7*\[Epsilon][1]^2 + 
            31*\[Epsilon][1]*\[Epsilon][2] + 84*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 7*\[Epsilon][2]^2 + 84*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2)) - (\[Epsilon][1]^2*\[Epsilon][2]^2*
           (\[Epsilon][1] + \[Epsilon][2])^3*(7*\[Epsilon][1]^2 + 
            31*\[Epsilon][1]*\[Epsilon][2] + 84*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 7*\[Epsilon][2]^2 + 84*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/2099520 - (\[Epsilon][1]*\[Epsilon][2]^3*
           (\[Epsilon][1] + \[Epsilon][2])^3*(7*\[Epsilon][1]^2 + 
            31*\[Epsilon][1]*\[Epsilon][2] + 84*Qrrc[0, 0]*\[Epsilon][1]^2*
             \[Epsilon][2] + 7*\[Epsilon][2]^2 + 84*Qrrc[0, 0]*\[Epsilon][1]*
             \[Epsilon][2]^2))/4199040)}}
 
tc = Del + (11*Del^2)/18 + (109*Del^3)/243 + (9389*Del^4)/26244 + 
     (88351*Del^5)/295245 + (823187*Del^6)/3188646 + 
     (68584051*Del^7)/301327047 + (632036939*Del^8)/3099363912 + 
     (5805190621*Del^9)/31381059609 + (4307941773209*Del^10)/25418658283290
