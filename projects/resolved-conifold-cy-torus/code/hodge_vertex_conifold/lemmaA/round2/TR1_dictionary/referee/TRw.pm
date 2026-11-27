package TRw;
# Referee's INDEPENDENT exact Eynard-Orantin recursion on the framed-vertex curve C_tau (round-2 TR1 referee).
# Written from scratch; shares no code with TRxi.pm / TRtau.pm.  Differences in method:
#   - global coordinate w, basis e_m(p) = dw_p/(w_p - a)^m (m >= 2) instead of the Laplace basis d xi_k or a zeta-basis;
#   - the local involution is computed directly as a power series s(u) in u = w - a (q = a+u, sigma(q) = a+s(u)) from the
#     ODE  s s' P(u) = u P(s),  P(u) = (a+u)(1+a+u)  (equivalent to X(a+s) = X(a+u), s = -u + O(u^2)); no zeta coordinate;
#   - residues are taken directly in u.
# Curve: X = log w - (tau+1) log(1+w), dX = (1 - tau w) dw/(w(1+w)), branch point a = 1/tau, B = dw dw'/(w-w')^2.
# Kernel (orientation of Zhou (23) and of the author's code):
#   K(w0,q) = (1/2) dw0 (1/(w0-q) - 1/(w0-sigma q)) / ((y(q) - y(sigma q)) dX(q))
#           = Sum_{m>=1} e_{m+1}(w0) V_m(u) / du,   V_m = (u^m - s^m) / (2 (y(a+u) - y(a+s)) X'(a+u)).
# omega_{g,n+1}(w0,J) = Res_{u=0} K(w0,q) E(q;J), E = omega_{g-1,n+2}(q,sq,J) + Sum' omega omega (omega_{0,2} = B included).
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
our ($tau, $a, $L, $ycase, $eparm, $sigma);
our (@s, @sp, @Ri, @Fh, %Wm, %Tq, %RG, %W, @SPOW, @RPOW, %Lc);

# ---- power series helpers (array refs, index = exponent, truncated at $L) ----
sub pmul { my ($A, $B, $N) = @_; $N //= $L; my @C = map { R(0) } 0..$N;
  for my $i (0..$#$A) { last if $i > $N; my $x = $A->[$i]; next if $x->is_zero;
    for my $j (0..$#$B) { last if $i + $j > $N; my $y = $B->[$j]; next if $y->is_zero; $C[$i+$j] += $x * $y; } }
  return \@C; }
sub pinv { my ($A, $N) = @_; $N //= $L; die "pinv: zero constant term" if $A->[0]->is_zero;
  my $i0 = R(1) / $A->[0]; my @B = ($i0);
  for my $n (1..$N) { my $t = R(0); for my $k (1..$n) { last if $k > $#$A; next if $A->[$k]->is_zero; $t += $A->[$k] * $B[$n-$k]; } push @B, -$t * $i0; }
  return \@B; }
sub pderiv { my $A = shift; return [ map { $A->[$_] * $_ } 1..$#$A ]; }
sub pscal { my ($c, $A) = @_; return [ map { $_ * $c } @$A ]; }
sub padd { my ($A, $B) = @_; my $n = $#$A > $#$B ? $#$A : $#$B; return [ map { ($A->[$_] // R(0)) + ($B->[$_] // R(0)) } 0..$n ]; }
sub plog1p { my ($T) = @_;   # log(1 + T(u)), T(0) = 0: integral of T'/(1+T)
  my @one = @$T; $one[0] = $one[0] + 1; my $d = pmul(pderiv($T), pinv(\@one, $L), $L - 1);
  return [ R(0), map { $d->[$_ - 1] / $_ } 1..$L ]; }

sub setup { my %o = @_;
  $tau = R($o{tau}); $L = $o{L}; $ycase = $o{y} // 'base'; $eparm = R($o{e} // 0);
  die "tau = 0, -1 excluded" if $tau->is_zero || $tau == -1;
  $a = R(1) / $tau; $sigma = $tau * ($tau + 1);
  %Wm = (); %Tq = (); %RG = (); %W = (); @SPOW = (); @RPOW = (); %Lc = ();
  # s(u) = Sum_{k>=1} s_k u^k from  s s' P(u) = u P(s);  p0 = a(1+a), p1 = 1+2a.
  my ($p0, $p1) = ($a * (1 + $a), 1 + 2 * $a);
  @s = (R(0), R(-1)); my @c = (R(0), R(0), R(1));     # c_n = [u^n] s^2  (s computed to order L+2)
  for my $m (2..$L + 2) {
    my $cm1 = $c[$m-1] // R(0); my $cm = $c[$m] // R(0);
    my $cnext = (2 * ($p1 * $s[$m-1] + $cm1) - $p1 * $m * $cm - ($m - 1) * $cm1) / ($p0 * ($m + 1));
    my $sig = R(0); for my $i (2..$m-1) { $sig += $s[$i] * $s[$m+1-$i]; }
    $s[$m] = ($sig - $cnext) / 2;
    # update c_{m+1} with the new s_m (consistency) and c_{m+2} partial terms are recomputed next round
    my $t = R(0); for my $i (1..$m) { my $j = $m + 1 - $i; next if $j < 1 || $j > $m; $t += $s[$i] * $s[$j]; } $c[$m+1] = $t;
    die "s recursion inconsistent at $m" unless $c[$m+1] == $cnext; }
  # check X(a+s(u)) - X(a+u) = O(u^(Lc+1)) with the Taylor coefficients X_j of X(a+u) - X(a)
  { my $Lc = $L < 16 ? $L : 16; my $t1 = $tau / ($tau + 1); my @Xj;
    for my $j (2..$Lc) { $Xj[$j] = R((-1) ** ($j - 1)) / $j * ($tau ** $j - ($tau + 1) * $t1 ** $j); }
    my @sv = map { $s[$_] } 0..$Lc; my @pw = @sv; my @acc = map { R(0) } 0..$Lc;
    for my $j (2..$Lc) { @pw = @{ pmul(\@pw, \@sv, $Lc) }; $acc[$_] += $Xj[$j] * $pw[$_] for 0..$Lc; $acc[$j] -= $Xj[$j]; }
    for my $i (0..$Lc) { die "involution check failed at order $i" unless $acc[$i]->is_zero; } }
  @sp = @{ pderiv(\@s) }; push @sp, R(0) while @sp <= $L;
  # s = -u (1 + r(u)):  onepr = -s/u ;  Ri = (1+r)^(-1)
  my @onepr = map { -($s[$_ + 1] // R(0)) } 0..$L;
  @Ri = @{ pinv(\@onepr, $L) };
  # y(a+u) - y(a+s(u)) = u * D1(u)
  my @U = (R(0), R(1), map { R(0) } 2..$L);
  my $Delta;
  my $deltaY0 = sub {   # -1/(a+u) + 1/(a+s) = (u - s)/((a+u)(a+s))
    my @den1 = map { R(0) } 0..$L; $den1[0] = $a; $den1[1] = R(1);
    my @den2 = map { ($s[$_] // R(0)) } 0..$L; $den2[0] = $a;
    my @num = map { ($U[$_] // R(0)) - ($s[$_] // R(0)) } 0..$L;
    return pmul(pmul(\@num, pinv(\@den1, $L), $L), pinv(\@den2, $L), $L); };
  my $deltaYB = sub {   # log(1 + u/(1+a)) - log(1 + s/(1+a))
    my $t1 = pscal(R(1) / (1 + $a), \@U); my $t2 = pscal(R(1) / (1 + $a), [ map { $s[$_] // R(0) } 0..$L ]);
    return padd(plog1p($t1), pscal(R(-1), plog1p($t2))); };
  if    ($ycase eq 'base') { $Delta = $deltaY0->(); }
  elsif ($ycase eq 'bkmp') { $Delta = $deltaYB->(); }
  elsif ($ycase eq 'zhou') { $Delta = pscal(R(-1), $deltaYB->()); }
  elsif ($ycase eq 'mix')  { $Delta = padd($deltaYB->(), pscal($eparm, $deltaY0->())); }
  else { die "unknown y case $ycase" }
  die "Delta(0) != 0" unless $Delta->[0]->is_zero;
  my @D1 = map { $Delta->[$_ + 1] // R(0) } 0..$L - 1; push @D1, R(0);
  die "dy(a) = 0 (irregular point)" if $D1[0]->is_zero;
  # F(u) = 1/(2 Delta X'(a+u)),  X'(a+u) = -tau u/P(u)  =>  F = u^(-2) Fh,  Fh = -P(u)/(2 tau D1(u))
  my @P = map { R(0) } 0..$L; $P[0] = $p0; $P[1] = $p1; $P[2] = R(1);
  @Fh = @{ pscal(R(-1) / (2 * $tau), pmul(\@P, pinv(\@D1, $L), $L)) };
  return 1; }

# cached powers s^m and Ri^q (incremental)
sub spow { my $m = shift; if (!@SPOW) { @SPOW = ([ R(1), map { R(0) } 1..$L ]); }
  my $sm = [ map { $s[$_] // R(0) } 0..$L ]; for my $k (scalar(@SPOW)..$m) { $SPOW[$k] = pmul($SPOW[$k-1], $sm, $L); } return $SPOW[$m]; }
sub ripow { my $q = shift; if (!@RPOW) { @RPOW = ([ R(1), map { R(0) } 1..$L ]); }
  for my $k (scalar(@RPOW)..$q) { $RPOW[$k] = pmul($RPOW[$k-1], \@Ri, $L); } return $RPOW[$q]; }
# Wm(m) = (u^m - s^m) Fh  (power series);  V_m(u) = u^(-2) Wm(u)
sub Wm { my $m = shift; return $Wm{$m} if exists $Wm{$m};
  my $pw = spow($m);
  my @d = map { -$pw->[$_] } 0..$L; $d[$m] += 1 if $m <= $L;
  $Wm{$m} = pmul(\@d, \@Fh, $L); return $Wm{$m}; }
# S_q(u) = s'(u) s(u)^(-q):  q >= 1: S_q[j] = (-1)^q [u^(j+q)] (s' Ri^q);  q <= 0: S_q[j] = [u^j] (s' s^(-q))
sub Tq { my $q = shift; return $Tq{$q} if exists $Tq{$q}; my $pw = [ R(1), map { R(0) } 1..$L ];
  if ($q >= 1) { $pw = ripow($q); my $t = pmul(\@sp, $pw, $L); $t = pscal(R(-1), $t) if $q % 2; $Tq{$q} = $t; }
  else { $pw = spow(-$q); $Tq{$q} = pmul(\@sp, $pw, $L); }
  return $Tq{$q}; }
sub Sq_coef { my ($q, $j) = @_; my $T = Tq($q); my $idx = $q >= 1 ? $j + $q : $j; return R(0) if $idx < 0; die "Sq precision ($q,$j)" if $idx > $L - 1; return $T->[$idx]; }
# RG(m,p,q) = Res_{u=0} V_m(u) u^(-p) S_q(u) du = Sum_i V_m[i] S_q[p-1-i],  V_m[i] = Wm[i+2]
sub RG { my ($m, $p, $q) = @_; my $key = "$m,$p,$q"; return $RG{$key} if exists $RG{$key};
  my $Wv = Wm($m); my $tot = R(0);
  my $qlow = $q >= 1 ? -$q : 0;           # S_q[j] = 0 for j < qlow
  for my $i (-2..$L - 2) { my $j = $p - 1 - $i; last if $j < $qlow; my $v = $Wv->[$i + 2]; next if $v->is_zero;
    die "V precision ($m,$p,$q)" if $i + 2 > $L - 1; $tot += $v * Sq_coef($q, $j); }
  $RG{$key} = $tot; return $tot; }

sub add_to { my ($h, $k, $v) = @_; return if $v->is_zero; if (exists $h->{$k}) { $h->{$k} += $v; delete $h->{$k} if $h->{$k}->is_zero; } else { $h->{$k} = $v->copy; } }
sub stable { my ($g, $n) = @_; return 2*$g - 2 + $n > 0; }
sub initial { %W = ();
  # (0,3): E = B(q,p1)B(sq,p2) + B(q,p2)B(sq,p1); coefficient (l1+1)(l2+1)[RG(m,-l1,-l2) + RG(m,-l2,-l1)] on e_{m+1} e_{l1+2} e_{l2+2}
  my %w03;
  for my $l1 (0..3) { for my $l2 (0..3) { for my $m (1..4) { my $r = ($l1 + 1) * ($l2 + 1) * (RG($m, -$l1, -$l2) + RG($m, -$l2, -$l1));
        add_to(\%w03, join(',', $m + 1, $l1 + 2, $l2 + 2), $r); } } }
  $W{"0,3"} = \%w03;
  # (1,1): E = B(q, sq) = s'(u)/(u - s)^2 du^2 = u^(-2) s' (2 + r)^(-2) du^2;  Res V_m E = [u^3] (Wm s' (2+r)^(-2))
  my @twopr = map { -($s[$_ + 1] // R(0)) } 0..$L - 1; push @twopr, R(0); $twopr[0] += 1;   # 1 + (1+r) = 2 + r
  my $i2 = pinv(pmul(\@twopr, \@twopr, $L), $L); my $f = pmul(\@sp, $i2, $L);
  my %w11; for my $m (1..6) { my $t = pmul(Wm($m), $f, $L); add_to(\%w11, $m + 1, $t->[3]); }
  $W{"1,1"} = \%w11; }
sub compute { my ($g, $np1) = @_; my $n = $np1 - 1; my %res;
  # (A) omega_{g-1,n+2}(q, sigma q, J)
  if ($g >= 1 && stable($g - 1, $n + 2)) { my $Wa = $W{($g-1) . ',' . ($n+2)} or die "missing";
    for my $key (keys %$Wa) { my @k = split /,/, $key; my $ma = shift @k; my $mb = shift @k; my $rest = join(',', @k); my $cf = $Wa->{$key};
      for my $m (1..$ma + $mb + 1) { my $r = RG($m, $ma, $mb); next if $r->is_zero; add_to(\%res, ($rest eq '' ? ($m+1) : ($m+1) . ",$rest"), $cf * $r); } } }
  # (B) stable products omega_h(q,I) omega_{g-h}(sigma q, J\I)
  for my $h (0..$g) { for my $mask (0..(2**$n) - 1) {
      my @I = grep { ($mask >> $_) & 1 } 0..$n-1; my @Ic = grep { !(($mask >> $_) & 1) } 0..$n-1;
      next unless stable($h, 1 + @I) && stable($g - $h, 1 + @Ic);
      my $W1 = $W{"$h," . (1 + @I)} or die "missing $h"; my $W2 = $W{($g - $h) . ',' . (1 + @Ic)} or die "missing";
      for my $k1 (keys %$W1) { my @p = split /,/, $k1; my $ma = shift @p; my $c1 = $W1->{$k1};
        for my $k2 (keys %$W2) { my @qq = split /,/, $k2; my $mb = shift @qq; my $cc = $c1 * $W2->{$k2};
          my @kJ; my ($ia, $ib) = (0, 0); for my $j (0..$n-1) { if (($mask >> $j) & 1) { push @kJ, $p[$ia++] } else { push @kJ, $qq[$ib++] } }
          my $rest = join(',', @kJ);
          for my $m (1..$ma + $mb + 1) { my $r = RG($m, $ma, $mb); next if $r->is_zero; add_to(\%res, ($rest eq '' ? ($m+1) : ($m+1) . ",$rest"), $cc * $r); } } } } }
  # (C) B(q,p_j) omega_{g,n}(sigma q, J\j) + omega_{g,n}(q, J\j) B(sigma q, p_j)
  if ($n >= 1 && stable($g, $n)) { my $Wc = $W{"$g,$n"} or die "missing";
    for my $jpos (0..$n-1) { for my $key (keys %$Wc) { my @k = split /,/, $key; my $mb = shift @k; my $cf = $Wc->{$key};
        for my $l (0..$mb + 1) { for my $m (1..$mb + 1 - $l) {
            my $r = ($l + 1) * (RG($m, -$l, $mb) + RG($m, $mb, -$l)); next if $r->is_zero;
            my @kJ = @k; splice(@kJ, $jpos, 0, $l + 2); add_to(\%res, join(',', $m + 1, @kJ), $cf * $r); } } } } }
  $W{"$g,$np1"} = \%res; }
sub run_all { my ($S, $gmax) = @_; initial(); my @order;
  for my $chi (2..3*$S) { for my $g (0..$gmax) { my $n = $chi - 2*$g + 2; next if $n < 1; next if $g + $n > $S;
      next if ($g == 0 && $n <= 3) || ($g == 1 && $n == 1); push @order, [$g, $n]; } }
  compute(@$_) for @order; return @order; }

# ---- x-expansion at w = 0 along x = w (1+w)^(-tau-1):  [x^mu] of the primitive of e_m vanishing at w = 0
#      = (1/mu) [w^(mu-1)] (w-a)^(-m) (1+w)^(mu(tau+1))   (Lagrange inversion)
sub binomR { my ($N, $k) = @_; my $r = R(1); for my $t (1..$k) { $r *= ($N - $t + 1) / $t; } return $r; }
sub binomI { my ($n, $k) = @_; return R(0) if $k < 0 || $k > $n; my $r = R(1); for my $t (1..$k) { $r *= R($n - $k + $t) / $t; } return $r; }
sub Lcoef { my ($m, $mu) = @_; my $key = "$m,$mu"; return $Lc{$key} if exists $Lc{$key};
  my $N = $mu * ($tau + 1); my $tot = R(0);
  for my $i (0..$mu-1) { my $k = $mu - 1 - $i; $tot += binomI($m + $i - 1, $i) / ($a ** $i) * binomR($N, $k); }
  $tot *= ((-$a) ** (-$m)) / $mu; $Lc{$key} = $tot; return $tot; }
sub xcoef { my ($Wh, @mu) = @_; my $tot = R(0);
  for my $key (keys %$Wh) { my @m = split /,/, $key; my $p = $Wh->{$key}->copy; for my $i (0..$#mu) { $p *= Lcoef($m[$i], $mu[$i]); last if $p->is_zero; } $tot += $p; }
  return $tot; }
sub Pi_mu { my $mu = shift; my $r = R(1); for my $t (1..$mu-1) { $r *= ($mu * $tau + $t) / $t; } return $r; }
1;
