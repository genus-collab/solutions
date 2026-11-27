package TRxi;
# Exact Eynard-Orantin recursion on the framed-vertex curve C_tau, run directly in the Laplace basis d xi_k
# (xi_0 = w/(1 - tau w), xi_k = D^k xi_0, D = d/dX;  [x^mu] xi_k = mu^k Pi_mu(tau)).  Same conventions as TRtau.pm:
#   X = log w - (tau+1) log(1+w), branch point a = 1/tau, B = dw dw'/(w-w')^2,
#   K(p0,q) = (1/2) Int_{sigma q}^{q} B(p0,.) / ((y(q) - y(sigma q)) dX(q)).
#   Local coordinate: X(a+u) - X(a) = -c2 zeta^2, c2 = tau^3/(2(tau+1)), so D = -(1/(2 c2 zeta)) d/dzeta.
#   Output of one recursion step: omega(p0) = Sum_l e_l(p0) (1/(2l+1)) Res zeta^(2l+1) Phi, e_l = (2l+1) zeta^(-2l-2) dzeta + hol.,
#   converted to the xi-basis with d xi_k = Sum_l T_{kl} e_l, T_{kl} = [zeta^(-2l-2)] xi_k'(zeta) / (2l+1).
#   B(q,p) = Sum_l e_l(p) zeta(q)^(2l) dzeta + (odd terms, which drop out by the linear loop equation).
# Correlators: W{"g,n"}{"k1,...,kn"} = coefficient of prod d xi_{k_i}(p_i).
use strict; use warnings;
use Math::BigRat; use Time::HiRes;
our $ZERO = Math::BigRat->new(0);
our $VERB = $ENV{TRVERB} // 0; our $T0 = 0;
sub R { Math::BigRat->new(@_) }

sub ps_mul { my ($a, $b, $L) = @_; my @c = map { R(0) } 0..$L;
  for my $i (0..$#$a) { last if $i > $L; my $ai = $a->[$i]; next if $ai->is_zero;
    for my $j (0..$#$b) { last if $i + $j > $L; my $bj = $b->[$j]; next if $bj->is_zero; $c[$i+$j] += $ai * $bj; } }
  return \@c; }
sub ps_inv { my ($a, $L) = @_; die "ps_inv: a0 = 0" if $a->[0]->is_zero; my $i0 = R(1) / $a->[0]; my @b = ($i0);
  for my $n (1..$L) { my $s = R(0); for my $k (1..$n) { last if $k > $#$a; next if $a->[$k]->is_zero; $s += $a->[$k] * $b[$n-$k]; } $b[$n] = -$s * $i0; }
  return \@b; }
sub ps_deriv { my $a = shift; return [ map { $a->[$_] * $_ } 1..$#$a ]; }
sub ps_negarg { my $a = shift; return [ map { ($_ % 2) ? -$a->[$_] : $a->[$_]->copy } 0..$#$a ]; }
# Laurent series { v => lowest exponent, c => [...] }
sub ls { my ($v, $c) = @_; return { v => $v, c => $c }; }
sub ls_top { my $A = shift; return $A->{v} + $#{$A->{c}}; }
sub ls_coef { my ($A, $e) = @_; my $i = $e - $A->{v}; return ($i < 0 || $i > $#{$A->{c}}) ? $ZERO : $A->{c}[$i]; }
sub ls_negarg { my $A = shift; my $v = $A->{v}; return ls($v, [ map { (($v + $_) % 2) ? -$A->{c}[$_] : $A->{c}[$_]->copy } 0..$#{$A->{c}} ]); }

our ($tau, $a, $c2, $kmax, $L, $ycase, $eparm, $sigma);
our (@u, @Xs, @Xm, @invYh, $F11, @xi, @T, @Tinv);
our (%rAB, %rC, %rhoAB, %rhoC);

sub setup { my %o = @_; $T0 = Time::HiRes::time();
  $tau = R($o{tau}); $kmax = $o{kmax}; $ycase = $o{y} // 'base'; $eparm = R($o{e} // 0);
  die "tau = 0 or -1 not allowed" if $tau->is_zero || $tau == -1;
  $sigma = $tau * ($tau + 1); $a = R(1) / $tau; $c2 = $tau ** 3 / (2 * ($tau + 1));
  %rAB = (); %rC = (); %rhoAB = (); %rhoC = ();
  # needed precision: X_a up to exponent 2 kmax + 4, so xi_0 up to exponent 4 kmax + 6; invYh up to 4 kmax + 6
  $L = 4 * $kmax + 10;
  # u(zeta) from the ODE  u u' = kappa zeta (a+u)(1+a+u),  kappa = tau^2/(tau+1)
  my $kap = $tau * $tau / ($tau + 1);
  @u = (R(0), R(1));
  for my $m (2..$L+2) { my $s = (1 + 2*$a) * $u[$m-1]; for my $i (1..$m-2) { $s += $u[$i] * $u[$m-1-$i]; }
    $s *= $kap; for my $i (2..$m-1) { $s -= ($m + 1 - $i) * $u[$i] * $u[$m+1-$i]; } $u[$m] = $s / ($m + 1); }
  # check X(a+u(zeta)) - X(a) + c2 zeta^2 = O(zeta^13) using the Taylor coefficients X_j
  { my $Lc = 12; my $t1 = $tau / ($tau + 1); my @X;
    for my $j (2..$Lc) { $X[$j] = R((-1) ** ($j - 1)) / $j * ($tau ** $j - ($tau + 1) * $t1 ** $j); }
    my @pw = (R(0), map { $u[$_] } 1..$Lc); my @acc = map { R(0) } 0..$Lc; my @cur = @pw;
    for my $j (2..$Lc) { @cur = @{ ps_mul(\@cur, \@pw, $Lc) }; $acc[$_] += $X[$j] * $cur[$_] for 0..$Lc; }
    $acc[2] += $c2; for my $i (0..$Lc) { die "u(zeta) check failed at order $i" unless $acc[$i]->is_zero; } }
  printf STDERR "  [setup] u done %.1f s\n", Time::HiRes::time() - $T0 if $VERB;
  # r = zeta/u ; xi_0 = -(1/tau) - (1/tau^2) u^(-1) = -(1/tau) - (1/tau^2) zeta^(-1) r
  my @uz = @u[1..$L+2]; my $r = ps_inv(\@uz, $L + 1);
  my @x0; for my $e (-1..$L) { my $cf = (R(-1) / ($tau * $tau)) * $r->[$e + 1]; $cf += R(-1) / $tau if $e == 0; push @x0, $cf; }
  @xi = (ls(-1, \@x0));
  for my $k (1..$kmax + 1) { my $P = $xi[$k-1]; my $v = $P->{v} - 2; my @c; my $f = R(-1) / (2 * $c2);
    for my $e ($v..ls_top($P) - 2) { push @c, $f * ($e + 2) * ls_coef($P, $e + 2); } $xi[$k] = ls($v, \@c); }
  # X_k = d xi_k / dzeta (Laurent, from -2k-2), and X_k(-zeta)
  @Xs = (); @Xm = ();
  for my $k (0..$kmax) { my $P = $xi[$k]; my @c; my $v = $P->{v} - 1;
    for my $e ($v..ls_top($P) - 1) { push @c, ($e + 1) * ls_coef($P, $e + 1); }
    $Xs[$k] = ls($v, \@c); $Xm[$k] = ls_negarg($Xs[$k]);
    for my $ee (-(2*$k+2)..-1) { next unless $ee % 2; die "X_$k has odd polar term" unless ls_coef($Xs[$k], $ee)->is_zero; } }
  printf STDERR "  [setup] xi_k done %.1f s\n", Time::HiRes::time() - $T0 if $VERB;
  # T_{kl} = [z^(-2l-2)] X_k / (2l+1), lower triangular; Tinv
  @T = (); @Tinv = ();
  for my $k (0..$kmax) { for my $l (0..$kmax) { $T[$k][$l] = $l <= $k ? ls_coef($Xs[$k], -2*$l - 2) / (2*$l + 1) : R(0); $Tinv[$k][$l] = R(0); } }
  for my $k (0..$kmax) { $Tinv[$k][$k] = R(1) / $T[$k][$k];
    for my $l (reverse 0..$k-1) { my $s = R(0); for my $m ($l+1..$k) { $s += $Tinv[$k][$m] * $T[$m][$l]; } $Tinv[$k][$l] = -$s / $T[$l][$l]; } }
  # y and invYh
  my @y = yseries($ycase);
  my @Yh = map { ($_ % 2) ? R(0) : 2 * $y[$_ + 1] } 0..$L;
  die "dy(a) = 0: irregular curve" if $Yh[0]->is_zero;
  @invYh = @{ ps_inv(\@Yh, $L) };
  printf STDERR "  [setup] y done %.1f s\n", Time::HiRes::time() - $T0 if $VERB;
  # F11 = B(q, sigma q)/dzeta^2 = -(1/4) z^-2 u'(z) u'(-z) / w(z)^2,  w = (u(z) - u(-z))/(2z);  only exponents <= 0 needed
  my $Lf = 4; my @w = map { (($_ + 1) % 2) ? $u[$_ + 1] : R(0) } 0..$Lf;
  my @up = map { $u[$_ + 1] * ($_ + 1) } 0..$Lf;
  my $f = ps_mul(ps_mul(\@up, ps_negarg(\@up), $Lf), ps_inv(ps_mul(\@w, \@w, $Lf), $Lf), $Lf);
  $F11 = ls(-2, [ map { -$_ / 4 } @$f ]);
  return 1; }

sub yseries { my $case = shift; my @y;
  my @den = map { $u[$_] / $a } 0..$L+1; $den[0] = R(1);
  my $inv = ps_inv(\@den, $L + 1); my @Y0 = map { -$_ / $a } @$inv;     # -1/(a+u)
  my @YB;
  if ($case ne 'base') {   # log(1 + u/(1+a)), constant dropped
    my @t = map { $u[$_] / (1 + $a) } 0..$L+1; my @t1 = @t; $t1[0] = R(1);
    my $inv1 = ps_inv(\@t1, $L + 1); my $dYB = ps_mul(ps_deriv(\@t), $inv1, $L);
    @YB = (R(0), map { $dYB->[$_ - 1] / $_ } 1..$L+1); }
  if    ($case eq 'base') { @y = @Y0; }
  elsif ($case eq 'bkmp') { @y = @YB; }
  elsif ($case eq 'zhou') { @y = map { -$_ } @YB; }
  elsif ($case eq 'mix')  { @y = map { $YB[$_] + $eparm * $Y0[$_] } 0..$L+1; }
  else { die "unknown y case $case" }
  return @y; }

# raw residue: [z^(-2k0)] ( X_x(z) X_y(-z) invYh ) / (2 c2 (2k0+1))   (e-coefficient k0 of the output)
sub rawAB { my ($k0, $x, $y) = @_; ($x, $y) = ($y, $x) if $x > $y; my $key = "$k0,$x,$y"; return $rAB{$key} if exists $rAB{$key};
  my ($A, $B) = ($Xs[$x], $Xm[$y]); my $s = R(0);
  # Sum over exponents i (of A), j (of B), m (of invYh) with i + j + m = -2k0
  for my $i ($A->{v}..ls_top($A)) { my $ac = ls_coef($A, $i); next if $ac->is_zero;
    for my $j ($B->{v}..ls_top($B)) { my $m = -2*$k0 - $i - $j; last if $m < 0; next if $m > $#invYh; my $bc = ls_coef($B, $j); next if $bc->is_zero;
      next if $invYh[$m]->is_zero; $s += $ac * $bc * $invYh[$m]; } }
  $s /= 2 * $c2 * (2*$k0 + 1); $rAB{$key} = $s; return $s; }
# rhoAB(k, x, y): xi-coefficient k of the output = Sum_{l >= k} rawAB(l,x,y) Tinv[l][k]
sub rhoAB { my ($k, $x, $y) = @_; ($x, $y) = ($y, $x) if $x > $y; my $key = "$k,$x,$y"; return $rhoAB{$key} if exists $rhoAB{$key};
  my $s = R(0); my $lmax = $x + $y + 2; $lmax = $kmax if $lmax > $kmax;
  for my $l ($k..$lmax) { my $t = $Tinv[$l][$k]; next if $t->is_zero; my $r = rawAB($l, $x, $y); next if $r->is_zero; $s += $r * $t; }
  $rhoAB{$key} = $s; return $s; }
# raw (C): [z^(-2k0)] ( z^(2l) (X_b(z) + X_b(-z)) invYh ) / (2 c2 (2k0+1)), nonzero only if k0 + l <= b + 1
sub rawC { my ($k0, $l, $b) = @_; my $key = "$k0,$l,$b"; return $rC{$key} if exists $rC{$key};
  my $s = R(0);
  for my $m (0..$#invYh) { my $e = -2*$k0 - 2*$l - $m; last if $e < $Xs[$b]{v}; next if $invYh[$m]->is_zero;
    my $pc = ls_coef($Xs[$b], $e) + ls_coef($Xm[$b], $e); next if $pc->is_zero; $s += $pc * $invYh[$m]; }
  $s /= 2 * $c2 * (2*$k0 + 1); $rC{$key} = $s; return $s; }
# rhoC(k, kp, b): output xi-coefficient k (slot p0) and k' (slot p_j)
sub rhoC { my ($k, $kp, $b) = @_; my $key = "$k,$kp,$b"; return $rhoC{$key} if exists $rhoC{$key};
  my $s = R(0);
  for my $k0 ($k..$b+1) { last if $k0 > $kmax; my $t0 = $Tinv[$k0][$k]; next if $t0->is_zero;
    for my $l ($kp..$b+1-$k0) { last if $l > $kmax; my $t1 = $Tinv[$l][$kp]; next if $t1->is_zero; my $r = rawC($k0, $l, $b); next if $r->is_zero; $s += $r * $t0 * $t1; } }
  $rhoC{$key} = $s; return $s; }

our %W;
sub add_to { my ($h, $k, $v) = @_; return if $v->is_zero; if (exists $h->{$k}) { $h->{$k} += $v; delete $h->{$k} if $h->{$k}->is_zero; } else { $h->{$k} = $v->copy; } }
sub stable { my ($g, $n) = @_; return 2*$g - 2 + $n > 0; }
sub initial { %W = ();
  # omega_{0,3}: e-coefficient [0,0,0] = invYh_0/c2, e_0 = d xi_0 / T_00
  my $t = $Tinv[0][0]; $W{"0,3"} = { "0,0,0" => $invYh[0] / $c2 * $t * $t * $t };
  my %e11; for my $k0 (0, 1) { my $s = R(0); for my $m (0..$#invYh) { my $e = -2*$k0 - $m; last if $e < -2; $s += ls_coef($F11, $e) * $invYh[$m]; }
    $e11{$k0} = -$s / (2 * $c2 * (2*$k0 + 1)); }
  my %w11; for my $k (0, 1) { my $s = R(0); for my $l ($k..1) { $s += $e11{$l} * $Tinv[$l][$k]; } add_to(\%w11, "$k", $s); }
  $W{"1,1"} = \%w11; }
sub compute { my ($g, $np1) = @_; my $n = $np1 - 1; my %res; my $K0 = 3*$g - 3 + $np1;
  if ($g >= 1 && stable($g - 1, $n + 2)) { my $Wa = $W{($g-1) . "," . ($n+2)} or die "missing";
    for my $key (keys %$Wa) { my @k = split /,/, $key; my $x = shift @k; my $y = shift @k; my $rest = join(",", @k); my $cf = $Wa->{$key};
      my $kk = $K0; for (@k) { $kk -= $_ } ;
      for my $k0 (0..$kk) { my $r = rhoAB($k0, $x, $y); next if $r->is_zero; add_to(\%res, ($rest eq '' ? "$k0" : "$k0,$rest"), $cf * $r); } } }
  for my $h (0..$g) { for my $mask (0..(2**$n) - 1) {
      my @I = grep { ($mask >> $_) & 1 } 0..$n-1; my @Ic = grep { !(($mask >> $_) & 1) } 0..$n-1;
      next unless stable($h, 1 + @I) && stable($g - $h, 1 + @Ic);
      my $W1 = $W{"$h," . (1 + @I)} or die "missing"; my $W2 = $W{($g-$h) . "," . (1 + @Ic)} or die "missing";
      for my $k1 (keys %$W1) { my @p = split /,/, $k1; my $x = shift @p; my $c1 = $W1->{$k1}; my $s1 = 0; $s1 += $_ for @p;
        for my $k2 (keys %$W2) { my @qq = split /,/, $k2; my $y = shift @qq; my $cc = $c1 * $W2->{$k2}; my $s2 = 0; $s2 += $_ for @qq;
          my @kJ; my ($ia, $ib) = (0, 0); for my $j (0..$n-1) { if (($mask >> $j) & 1) { push @kJ, $p[$ia++] } else { push @kJ, $qq[$ib++] } }
          my $rest = join(",", @kJ); my $kk = $K0 - $s1 - $s2;
          for my $k0 (0..$kk) { my $r = rhoAB($k0, $x, $y); next if $r->is_zero; add_to(\%res, ($rest eq '' ? "$k0" : "$k0,$rest"), $cc * $r); } } } } }
  if ($n >= 1 && stable($g, $n)) { my $Wc = $W{"$g,$n"} or die "missing";
    for my $jpos (0..$n-1) { for my $key (keys %$Wc) { my @k = split /,/, $key; my $b = shift @k; my $cf = $Wc->{$key};
        for my $k0 (0..$b+1) { for my $kp (0..$b+1-$k0) { my $r = rhoC($k0, $kp, $b); next if $r->is_zero;
            my @kJ = @k; splice(@kJ, $jpos, 0, $kp); add_to(\%res, join(",", $k0, @kJ), $cf * $r); } } } } }
  $W{"$g,$np1"} = \%res; }
sub run_all { my ($S, $gmax, $nmax_top) = @_; initial(); my @order;
  for my $chi (2..3*$S) { for my $g (0..$gmax) { my $n = $chi - 2*$g + 2; next if $n < 1; next if $g + $n > $S;
      next if ($g == 0 && $n <= 3) || ($g == 1 && $n == 1); push @order, [$g, $n]; } }
  for my $o (@order) { my $t = Time::HiRes::time(); compute(@$o);
    printf STDERR "  [rec] (%d,%d): %d terms, %.1f s\n", @$o, scalar(keys %{$W{"$o->[0],$o->[1]"}}), Time::HiRes::time() - $t if $VERB; }
  return @order; }
sub Pi_mu { my $mu = shift; my $r = R(1); for my $t (1..$mu-1) { $r *= ($mu * $tau + $t) / $t; } return $r; }
1;
