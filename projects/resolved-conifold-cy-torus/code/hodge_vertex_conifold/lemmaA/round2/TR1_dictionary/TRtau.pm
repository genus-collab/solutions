package TRtau;
# Exact Eynard-Orantin recursion on the framed-vertex curve C_tau (Math::BigRat, rational tau).
#   Curve: P^1 with coordinate w, X = log x = log w - (tau+1) log(1+w), dX = (1 - tau w) dw / (w (1+w)),
#   single branch point a = 1/tau, B = dw dw' / (w - w')^2.
#   y is one of: 'base'  y = -1/w                     (the base curve C_tau of P2)
#                'bkmp'  y = log(1+w)                 (BKMP / Zhou background, up to the sign y -> -y)
#                'mix'   y = log(1+w) + e * (-1/w)    (the family Y_e of the referee route; e rational)
#                'zhou'  y = -log(1+w)                (Zhou's convention ln y, y_Zhou = 1/(1+w))
#   Kernel convention (as in review_checks/referee1_trC1.pl and Zhou arXiv:0910.4320 (21)-(23)):
#     K(p0,q) = (1/2) Int_{sigma q}^{q} B(p0,.) / ((y(q) - y(sigma q)) dX(q)),
#     omega_{g,n+1}(p0,J) = Res_{q->a} K(p0,q) [ omega_{g-1,n+2}(q,sigma q,J) + Sum' omega omega ].
#   Local coordinate zeta at a:  X(a+u) - X(a) = -c2 zeta^2,  c2 = tau^3/(2(tau+1)),  zeta = u + O(u^2), sigma: zeta -> -zeta.
#   Local basis e_k(p) := Res_{q->a} B(p,q) zeta(q)^(-2k-1) = (2k+1) Sum_{j=1}^{2k+1} A_{k,j} dw/(w-a)^(j+1),  A_{k,j} = [zeta^(2k+1)] u^j
#     (in zeta: e_k = E_k(zeta) d zeta, E_k = (2k+1) zeta^(-2k-2) + holomorphic).
#   Laplace basis: xi_0 = w/(1 - tau w), xi_k = D^k xi_0, D = d/dX;  [x^mu] xi_k = mu^k Pi_mu(tau),  Pi_mu = C(mu(tau+1)-1, mu-1).
#   Correlators are stored as hashes  W{"g,n"}{"k1,...,kn"} = coefficient of prod e_{k_i}(p_i) (all orderings stored).
use strict; use warnings;
use Math::BigRat; use Time::HiRes;
our $ZERO = Math::BigRat->new(0);
our $VERB = $ENV{TRVERB} // 0; our $T0 = 0;
sub R { Math::BigRat->new(@_) }

# ---------------- power series: array refs, index = exponent ----------------
sub ps_mul { my ($a, $b, $L) = @_; my @c = map { R(0) } 0..$L;
  for my $i (0..$#$a) { last if $i > $L; my $ai = $a->[$i]; next if $ai->is_zero;
    for my $j (0..$#$b) { last if $i + $j > $L; my $bj = $b->[$j]; next if $bj->is_zero; $c[$i+$j] += $ai * $bj; } }
  return \@c; }
sub ps_inv { my ($a, $L) = @_; die "ps_inv: a0 = 0" if $a->[0]->is_zero; my $i0 = R(1) / $a->[0]; my @b = ($i0);
  for my $n (1..$L) { my $s = R(0); for my $k (1..$n) { last if $k > $#$a; next if $a->[$k]->is_zero; $s += $a->[$k] * $b[$n-$k]; } $b[$n] = -$s * $i0; }
  return \@b; }
sub ps_powr { my ($q, $alpha, $L) = @_;   # q^alpha for q0 = 1, alpha rational
  die "ps_powr: q0 != 1" unless $q->[0] == 1; my @f = (R(1));
  for my $m (1..$L) { my $s = R(0);
    for my $k (1..$m) { last if $k > $#$q; next if $q->[$k]->is_zero; $s += ($alpha * $k - ($m - $k)) * $q->[$k] * $f[$m-$k]; }
    $f[$m] = $s / $m; }
  return \@f; }
sub ps_deriv { my $a = shift; return [ map { $a->[$_] * $_ } 1..$#$a ]; }
sub ps_negarg { my $a = shift; return [ map { ($_ % 2) ? -$a->[$_] : $a->[$_]->copy } 0..$#$a ]; }

# ---------------- Laurent series: { v => lowest exponent, c => [coefficients] } ----------------
sub ls { my ($v, $c) = @_; return { v => $v, c => $c }; }
sub ls_top { my $A = shift; return $A->{v} + $#{$A->{c}}; }
sub ls_coef { my ($A, $e) = @_; my $i = $e - $A->{v}; return ($i < 0 || $i > $#{$A->{c}}) ? $ZERO : $A->{c}[$i]; }
sub ls_mul { my ($A, $B, $M) = @_; my $v = $A->{v} + $B->{v}; my $top = ls_top($A) + ls_top($B);
  # valid only up to min(topA + vB, topB + vA)
  my $valid = ls_top($A) + $B->{v}; my $v2 = ls_top($B) + $A->{v}; $valid = $v2 if $v2 < $valid; $M = $valid if $M > $valid;
  my $len = $M - $v; return ls($v, []) if $len < 0;
  my @c = map { R(0) } 0..$len; my ($a, $b) = ($A->{c}, $B->{c});
  for my $i (0..$#$a) { last if $i > $len; my $ai = $a->[$i]; next if $ai->is_zero;
    for my $j (0..$#$b) { last if $i + $j > $len; my $bj = $b->[$j]; next if $bj->is_zero; $c[$i+$j] += $ai * $bj; } }
  return ls($v, \@c); }
sub ls_negarg { my $A = shift; my $v = $A->{v};
  return ls($v, [ map { (($v + $_) % 2) ? -$A->{c}[$_] : $A->{c}[$_]->copy } 0..$#{$A->{c}} ]); }
sub ls_add { my ($A, $B) = @_; my $v = $A->{v} < $B->{v} ? $A->{v} : $B->{v}; my $top = ls_top($A) < ls_top($B) ? ls_top($A) : ls_top($B);
  return ls($v, [ map { ls_coef($A, $v + $_) + ls_coef($B, $v + $_) } 0..($top - $v) ]); }
sub ls_scal { my ($s, $A) = @_; return ls($A->{v}, [ map { $_ * $s } @{$A->{c}} ]); }
sub ls_from_ps { my ($p, $shift) = @_; return ls($shift, [ map { $_->copy } @$p ]); }

# ---------------- curve set-up ----------------
our ($tau, $a, $c2, $kmax, $ME, $L, $ycase, $eparm);
our (@u, @up, @A, @E, @Em, @invYh, $F11, @xi, @T, @Tinv, $sigma);
our (%rhoAB, %rhoC, %Pab);

sub setup { my %o = @_;
  $T0 = Time::HiRes::time();
  $tau = R($o{tau}); $kmax = $o{kmax}; $ycase = $o{y} // 'base'; $eparm = R($o{e} // 0);
  die "tau = 0 or -1 not allowed" if $tau->is_zero || $tau == -1;
  $sigma = $tau * ($tau + 1);
  $a = R(1) / $tau;
  $ME = 2 * $kmax + 6;            # Laurent series E_k kept up to this exponent
  $L  = 4 * $kmax + 12;           # power series length
  %rhoAB = (); %rhoC = (); %Pab = ();
  # Taylor coefficients X_j of X(a+u) - X(a): X_j = ((-1)^(j-1)/j) (tau^j - (tau+1) (tau/(tau+1))^j)
  my @X; my $t1 = $tau / ($tau + 1);
  for my $j (2..$L+3) { $X[$j] = R((-1) ** ($j - 1)) / $j * ($tau ** $j - ($tau + 1) * $t1 ** $j); }
  $c2 = -$X[2];
  die "c2 mismatch" unless $c2 == $tau ** 3 / (2 * ($tau + 1));
  # u(zeta) from X(a+u) - X(a) = -c2 zeta^2.  Differentiating, X'(a+u) = -tau u/((a+u)(1+a+u)) gives the ODE
  #   u u' = kappa zeta (a+u)(1+a+u),  kappa = 2 c2/tau = tau^2/(tau+1),  u = zeta + O(zeta^2);
  # coefficient of zeta^m (m >= 2):  (m+1) u_m + Sum_{i=2}^{m-1} (m+1-i) u_i u_{m+1-i} = kappa [ (1+2a) u_{m-1} + Sum_{i+j=m-1} u_i u_j ].
  my $kap = 2 * $c2 / $tau; die "kappa" unless $kap * $a * (1 + $a) == 1;
  @u = (R(0), R(1));
  for my $m (2..$L+1) { my $s = (1 + 2*$a) * $u[$m-1];
    for my $i (1..$m-2) { $s += $u[$i] * $u[$m-1-$i]; }
    $s *= $kap;
    for my $i (2..$m-1) { $s -= ($m + 1 - $i) * $u[$i] * $u[$m+1-$i]; }
    $u[$m] = $s / ($m + 1); }
  # check against the Taylor coefficients of X: X(a+u(z)) - X(a) + c2 z^2 = O(z^(L+2))
  { my $Lc = $L + 1 < 14 ? $L + 1 : 14;   # (checked to order 14 only, to save time)
    my @pw = (R(0), map { $u[$_] } 1..$Lc); my @acc = map { R(0) } 0..$Lc; my @cur = @pw;
    for my $j (2..$Lc) { @cur = @{ ps_mul(\@cur, \@pw, $Lc) }; $acc[$_] += $X[$j] * $cur[$_] for 0..$Lc; }
    $acc[2] += $c2; for my $i (0..$Lc) { die "u(zeta) check failed at order $i" unless $acc[$i]->is_zero; } }
  printf STDERR "  [setup] u(zeta) done %.1f s\n", Time::HiRes::time() - $T0 if $VERB;
  @up = @{ ps_deriv(\@u) };                         # u'(zeta), length L+1
  # r = zeta/u (power series, r0 = 1)
  my @uz = @u[1..$#u]; my $r = ps_inv(\@uz, $L);
  # powers r^m, m = 1 .. 2 kmax + 2 ;  u^(-m) = zeta^(-m) r^m
  my @rp; $rp[0] = [ R(1), map { R(0) } 1..$L ]; for my $m (1..2*$kmax+2) { $rp[$m] = ps_mul($rp[$m-1], $r, $L); }
  # positive powers u^j truncated at 2 kmax + 1: A_{k,j} = [zeta^(2k+1)] u^j
  my @Up; $Up[1] = [ map { $u[$_] } 0..2*$kmax+1 ];
  for my $j (2..2*$kmax+1) { $Up[$j] = ps_mul($Up[$j-1], $Up[1], 2*$kmax+1); }
  @A = ();
  for my $k (0..$kmax) { for my $j (1..2*$k+1) { $A[$k][$j] = $Up[$j][2*$k+1]; } }
  printf STDERR "  [setup] powers done %.1f s\n", Time::HiRes::time() - $T0 if $VERB;
  # E_k(zeta) = u'(zeta) (2k+1) Sum_j A_{k,j} zeta^(-j-1) r^(j+1)
  @E = (); @Em = ();
  my $upl = ls_from_ps(\@up, 0);
  for my $k (0..$kmax) { my $v = -(2*$k + 2); my @c;
    for my $e ($v..$ME) { my $s = R(0);
      for my $j (1..2*$k+1) { my $idx = $e + $j + 1; next if $idx < 0; die "rp too short" if $idx > $L; $s += $A[$k][$j] * $rp[$j+1][$idx]; }
      push @c, $s; }
    my $S = ls($v, \@c);
    $E[$k] = ls_scal(R(2*$k + 1), ls_mul($upl, $S, $ME));
    die "E_$k leading" unless ls_coef($E[$k], $v) == 2*$k + 1;
    $Em[$k] = ls_negarg($E[$k]); }
  printf STDERR "  [setup] E_k done %.1f s\n", Time::HiRes::time() - $T0 if $VERB;
  # y(zeta)
  my @y = yseries($ycase);
  printf STDERR "  [setup] y done %.1f s\n", Time::HiRes::time() - $T0 if $VERB;
  # Yh = (y(z) - y(-z))/z (even), invYh
  my @Yh = map { ($_ % 2) ? R(0) : 2 * $y[$_ + 1] } 0..$L-1;
  die "dy(a) = 0: irregular curve" if $Yh[0]->is_zero;
  @invYh = @{ ps_inv(\@Yh, $L - 1) };
  # F11 = B(q, sigma q)/dzeta^2 = -u'(z) u'(-z) / (u(z) - u(-z))^2 = -(1/4) z^-2 u'(z)u'(-z)/w(z)^2, w = (u(z)-u(-z))/(2z)
  my @w = map { $u[$_ + 1] * ((($_ + 1) % 2) ? 1 : 0) } 0..$L-1;
  my $w2i = ps_inv(ps_mul(\@w, \@w, $L - 1), $L - 1);
  my $uu = ps_mul(\@up, ps_negarg(\@up), $L - 1);
  my $f = ps_mul($uu, $w2i, $L - 1);
  $F11 = ls(-2, [ map { -$_ / 4 } @$f ]);
  # Laplace basis xi_k as Laurent series: xi_0 = -(1/tau) - (1/tau^2) u^-1 ;  xi_{k+1} = -(1/(2 c2 z)) d/dz xi_k
  my $top = 2 * $kmax + 4;
  my @x0; for my $e (-1..$top) { my $cf = (-1/($tau*$tau)) * $rp[1][$e + 1]; $cf += -1 / $tau if $e == 0; push @x0, $cf; }
  @xi = (ls(-1, \@x0));
  for my $k (1..$kmax) { my $P = $xi[$k-1]; my $v = $P->{v} - 2; my @c;
    for my $e ($v..ls_top($P) - 2) { push @c, (-1 / (2 * $c2)) * ($e + 2) * ls_coef($P, $e + 2); }
    $xi[$k] = ls($v, \@c); }
  # T_{k,l}: d xi_k = Sum_l T_{k,l} e_l,  T_{k,l} = [z^(-2l-2)] (xi_k'(z)) / (2l+1)  ( = -2 c2 [z^(-2l-3)] xi_{k+1} / (2l+1) )
  @T = ();
  for my $k (0..$kmax) { for my $l (0..$kmax) {
      # [z^(-2l-2)] xi_k' = e c_e with e = -2l-1 (c_e = [z^e] xi_k)
      my $e = -2*$l - 1; $T[$k][$l] = ($l <= $k) ? ls_coef($xi[$k], $e) * $e / (2*$l + 1) : R(0); }
    # odd polar parts check: xi_k must have only odd negative powers
    for my $ee (-(2*$k+1)..-1) { next if $ee % 2; die "xi_$k has even polar term z^$ee" unless ls_coef($xi[$k], $ee)->is_zero; } }
  # inverse of the lower-triangular T
  @Tinv = ();
  for my $k (0..$kmax) { for my $l (0..$kmax) { $Tinv[$k][$l] = R(0); } }
  for my $k (0..$kmax) { $Tinv[$k][$k] = R(1) / $T[$k][$k];
    for my $l (reverse 0..$k-1) { my $s = R(0); for my $m ($l+1..$k) { $s += $Tinv[$k][$m] * $T[$m][$l]; } $Tinv[$k][$l] = -$s / $T[$l][$l]; } }
  printf STDERR "  [setup] all done %.1f s\n", Time::HiRes::time() - $T0 if $VERB;
  return 1; }

sub yseries { my $case = shift; my @y;
  my @Y0; my @YB;
  # Y0 = -1/(a+u) = -(1/a) / (1 + u/a)
  my @den = map { $u[$_] / $a } 0..$L+1; $den[0] = R(1);
  my $inv = ps_inv(\@den, $L + 1); @Y0 = map { -$_ / $a } @$inv;
  # YB = log(1 + t), t = u/(1+a) (constant dropped):  YB' = t'/(1+t)
  my @t = map { $u[$_] / (1 + $a) } 0..$L+1; my @t1 = @t; $t1[0] = R(1);
  my $inv1 = ps_inv(\@t1, $L + 1); my $dYB = ps_mul(ps_deriv(\@t), $inv1, $L);
  @YB = (R(0), map { $dYB->[$_ - 1] / $_ } 1..$L+1);
  if    ($case eq 'base') { @y = @Y0; }
  elsif ($case eq 'bkmp') { @y = @YB; }
  elsif ($case eq 'zhou') { @y = map { -$_ } @YB; }
  elsif ($case eq 'mix')  { @y = map { $YB[$_] + $eparm * $Y0[$_] } 0..$L+1; }
  else { die "unknown y case $case" }
  return @y; }

# ---------------- residue tables ----------------
# rhoAB(k0,a,b) = [z^(-2k0)] ( E_a(z) E_b(-z) invYh ) / (2 c2 (2k0+1))
sub Pab { my ($x, $y) = @_; my $key = "$x,$y"; return $Pab{$key} if exists $Pab{$key};
  my $P = ls_mul($E[$x], $Em[$y], 0); $Pab{$key} = $P; return $P; }
sub rhoAB { my ($k0, $x, $y) = @_; ($x, $y) = ($y, $x) if $x > $y;   # symmetric in (a,b) at even powers
  my $key = "$k0,$x,$y"; return $rhoAB{$key} if exists $rhoAB{$key};
  my $P = Pab($x, $y); my $s = R(0);
  for my $m (0..$#invYh) { my $e = -2*$k0 - $m; last if $e < $P->{v}; next if $invYh[$m]->is_zero; my $pc = ls_coef($P, $e); next if $pc->is_zero; $s += $pc * $invYh[$m]; }
  $s /= 2 * $c2 * (2*$k0 + 1); $rhoAB{$key} = $s; return $s; }
# rhoC(k0,l,b) = [z^(-2k0)] ( z^(2l) (E_b(z)+E_b(-z)) invYh ) / (2 c2 (2k0+1));  nonzero only if k0 + l <= b + 1
sub rhoC { my ($k0, $l, $b) = @_; my $key = "$k0,$l,$b"; return $rhoC{$key} if exists $rhoC{$key};
  my $s = R(0);
  for my $m (0..$#invYh) { my $e = -2*$k0 - 2*$l - $m; last if $e < $E[$b]{v}; next if $invYh[$m]->is_zero;
    my $pc = ls_coef($E[$b], $e) + ls_coef($Em[$b], $e); next if $pc->is_zero; $s += $pc * $invYh[$m]; }
  $s /= 2 * $c2 * (2*$k0 + 1); $rhoC{$key} = $s; return $s; }

# ---------------- recursion ----------------
our %W;
sub add_to { my ($h, $k, $v) = @_; return if $v->is_zero; if (exists $h->{$k}) { $h->{$k} += $v; delete $h->{$k} if $h->{$k}->is_zero; } else { $h->{$k} = $v->copy; } }
sub stable { my ($g, $n) = @_; return 2*$g - 2 + $n > 0; }
sub initial {
  %W = ();
  $W{"0,3"} = { "0,0,0" => $invYh[0] / $c2 };
  my %w11; for my $k0 (0, 1) { my $s = R(0);
    for my $m (0..$#invYh) { my $e = -2*$k0 - $m; last if $e < -2; $s += ls_coef($F11, $e) * $invYh[$m]; }
    add_to(\%w11, "$k0", -$s / (2 * $c2 * (2*$k0 + 1))); }
  $W{"1,1"} = \%w11; }
sub compute { my ($g, $np1) = @_; my $n = $np1 - 1; my %res; my $K0 = 3*$g - 3 + $np1;
  # (A) omega_{g-1,n+2}(q, sigma q, J)
  if ($g >= 1 && stable($g - 1, $n + 2)) { my $Wa = $W{($g-1) . "," . ($n+2)} or die "missing " . ($g-1) . "," . ($n+2);
    for my $key (keys %$Wa) { my @k = split /,/, $key; my $x = shift @k; my $y = shift @k; my $rest = join(",", @k); my $cf = $Wa->{$key};
      for my $k0 (0..$K0) { my $r = rhoAB($k0, $x, $y); next if $r->is_zero; add_to(\%res, ($rest eq '' ? "$k0" : "$k0,$rest"), $cf * $r); } } }
  # (B) stable products
  for my $h (0..$g) { for my $mask (0..(2**$n) - 1) {
      my @I = grep { ($mask >> $_) & 1 } 0..$n-1; my @Ic = grep { !(($mask >> $_) & 1) } 0..$n-1;
      next unless stable($h, 1 + @I) && stable($g - $h, 1 + @Ic);
      my $W1 = $W{"$h," . (1 + @I)} or die "missing $h," . (1 + @I); my $W2 = $W{($g-$h) . "," . (1 + @Ic)} or die "missing";
      for my $k1 (keys %$W1) { my @p = split /,/, $k1; my $x = shift @p; my $c1 = $W1->{$k1};
        for my $k2 (keys %$W2) { my @qq = split /,/, $k2; my $y = shift @qq; my $cc = $c1 * $W2->{$k2};
          my @kJ; my ($ia, $ib) = (0, 0); for my $j (0..$n-1) { if (($mask >> $j) & 1) { push @kJ, $p[$ia++] } else { push @kJ, $qq[$ib++] } }
          my $rest = join(",", @kJ);
          for my $k0 (0..$K0) { my $r = rhoAB($k0, $x, $y); next if $r->is_zero; add_to(\%res, ($rest eq '' ? "$k0" : "$k0,$rest"), $cc * $r); } } } } }
  # (C) omega_{0,2}(q,p_j) omega_{g,n}(sigma q, J-j) + (q <-> sigma q)
  if ($n >= 1 && stable($g, $n)) { my $Wc = $W{"$g,$n"} or die "missing $g,$n";
    for my $jpos (0..$n-1) { for my $key (keys %$Wc) { my @k = split /,/, $key; my $b = shift @k; my $cf = $Wc->{$key};
        for my $l (0..$b+1) { for my $k0 (0..$b+1-$l) { my $r = rhoC($k0, $l, $b); next if $r->is_zero;
            my @kJ = @k; splice(@kJ, $jpos, 0, $l); add_to(\%res, join(",", $k0, @kJ), $cf * $r); } } } } }
  $W{"$g,$np1"} = \%res; }

# all (g,n) with g + n <= S (and 2g-2+n >= 1), ordered by 2g-2+n
sub run_all { my ($S, $gmax) = @_; initial();
  my @order;
  for my $chi (2..3*$S) { for my $g (0..$gmax) { my $n = $chi - 2*$g + 2; next if $n < 1; next if $g + $n > $S;
      next if ($g == 0 && $n <= 3) || ($g == 1 && $n == 1); push @order, [$g, $n]; } }
  for my $o (@order) { compute(@$o); }
  return @order; }

# ---------------- conversion to the Laplace basis d xi_k ----------------
# C^e[l_1..l_n] = Sum_k C^xi[k] prod T_{k_i l_i}  =>  C^xi = C^e (Tinv) slotwise: C^xi[k] = Sum_l C^e[l] prod Tinv_{l_i k_i}
sub to_xi { my $We = shift; my %cur = %$We; my $n;
  for my $key (keys %cur) { $n = scalar(split /,/, $key); last; }
  return {} unless defined $n;
  for my $slot (0..$n-1) { my %nxt;
    for my $key (keys %cur) { my @l = split /,/, $key; my $lv = $l[$slot];
      for my $k (0..$lv) { my $t = $Tinv[$lv][$k]; next if $t->is_zero; my @kk = @l; $kk[$slot] = $k; add_to(\%nxt, join(",", @kk), $cur{$key} * $t); } }
    %cur = %nxt; }
  return \%cur; }

# ---------------- x-expansion at w = 0 (independent of the basis change) ----------------
sub binom_r { my ($N, $i) = @_; my $r = R(1); for my $t (1..$i) { $r *= ($N - $t + 1) / $t; } return $r; }
sub binom_i { my ($n, $k) = @_; return R(0) if $k < 0 || $k > $n; my $r = R(1); for my $t (1..$k) { $r *= R($n - $k + $t) / $t; } return $r; }
my %xcache;
# [x^mu] of the primitive of e_l (expansion at w = 0 along x = w/(1+w)^(tau+1))
sub xcoef_e { my ($l, $mu) = @_; my $key = "$l,$mu"; return $xcache{$key} if exists $xcache{$key};
  my $N = $mu * ($tau + 1); my $s = R(0);
  for my $j (1..2*$l+1) { my $acc = R(0);
    for my $i (0..$mu-1) { my $m = $mu - 1 - $i; $acc += binom_r($N, $i) * binom_i($j + $m, $m) / ($a ** $m); }
    $s += $A[$l][$j] * $acc * (-$a) ** (-$j - 1); }
  $s *= R(2*$l + 1) / $mu; $xcache{$key} = $s; return $s; }
sub Pi_mu { my $mu = shift; my $r = R(1); for my $t (1..$mu-1) { $r *= ($mu * $tau + $t) / $t; } return $r; }
sub xcoef { my ($We, @mu) = @_; my $tot = R(0);
  for my $key (keys %$We) { my @l = split /,/, $key; my $p = $We->{$key}->copy; for my $i (0..$#mu) { $p *= xcoef_e($l[$i], $mu[$i]); } $tot += $p; }
  return $tot; }
1;
