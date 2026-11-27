package R2TR;
# Referee 2: an independent exact implementation of the Eynard-Orantin recursion on C_tau (Math::BigRat, rational tau).
# Written from scratch; it shares no code with the author's TRxi.pm / TRtau.pm.  Same TR conventions as the brief:
#   X = log w - (tau+1) log(1+w), branch point a = 1/tau, B = dw dw'/(w-w')^2,
#   K(p0,q) = (1/2) Int_{sigma q}^{q} B(p0,.) / ((y(q) - y(sigma q)) dX(q)),  omega_{g,n+1} = Res_{q->a} K E_{g,n+1}.
# Method (differs from the author's on purpose):
#  * zeta from the Taylor series of X:  X(a+u) - X(a) = -c2 zeta^2  =>  zeta = u (1 + q(u))^(1/2), u(zeta) by Lagrange
#    inversion (author: ODE u u' = kappa zeta (a+u)(1+a+u));
#  * basis E_l (l >= 0): the global differential with principal part zeta^(-2l-2) dzeta at a and no other pole,
#      E_l(w) = (1/(2l+1)) Sum_{j=0}^{2l} (j+1) [zeta^(2l)](u' u^j) dw/(w-a)^(j+2)       (Res_{q->a} B(w,q) zeta^(-2l-1)/(2l+1)),
#    whose local expansions are computed from this global form (author: xi-basis, xi_k = D^k xi_0);
#  * output of a recursion step = principal part of E/omega~ (Lemma 1(a)):  omega(p0) = Sum_l E_l(p0) [zeta^(-2l-2)](E/omega~);
#  * B(q,p) = Sum_m (m+1) zeta^m dzeta Ehat_m(p), Ehat_{2l} = E_l; the odd m drop (their contribution is O(zeta^2) because
#    every correlator stored in the E-basis has holomorphic odd part);
#  * optional dual numbers: y -> y + eps*eta (exact first-order variation, for the variational formula);
#  * x-expansion by Lagrange inversion from the global forms (no xi-basis), for an end-to-end test of P2.
use strict; use warnings;
use Math::BigRat; use Time::HiRes qw(time);
sub R { Math::BigRat->new(@_) }
our $Z = R(0);
our $VERB = $ENV{R2VERB} // 0;
our ($tau, $a, $c2, $LMAX, $NU, $DUAL, $sigma);
our (@u, @up, @e, @eev, @invYh, @invYhE, @F11, @T, @Ti, @Aglob);
our (%PAB, %PC);
our %W;

# ---------- power series (array refs, index = exponent >= 0) ----------
sub pmul { my ($p, $q, $N) = @_; my @c = map { R(0) } 0..$N;
  for my $i (0..$#$p) { last if $i > $N; next if $p->[$i]->is_zero; my $pi = $p->[$i];
    for my $j (0..$#$q) { last if $i + $j > $N; next if $q->[$j]->is_zero; $c[$i+$j] += $pi * $q->[$j]; } }
  return \@c; }
sub ppow { my ($q, $al, $N) = @_;   # (1 + q)^al, q = [0, q1, q2, ...]
  my @f = (R(1)); $al = R($al) unless ref $al;
  for my $m (1..$N) { my $s = R(0);
    for my $k (1..$m) { last if $k > $#$q; next if $q->[$k]->is_zero; $s += ($al * $k - ($m - $k)) * $q->[$k] * $f[$m-$k]; }
    push @f, $s / $m; }
  return \@f; }
sub pinv { my ($p, $N) = @_; die "pinv: p0 = 0" if $p->[0]->is_zero; my $i0 = R(1) / $p->[0]; my @b = ($i0);
  for my $n (1..$N) { my $s = R(0); for my $k (1..$n) { last if $k > $#$p; next if $p->[$k]->is_zero; $s += $p->[$k] * $b[$n-$k]; } push @b, -$s * $i0; }
  return \@b; }
sub pder { my $p = shift; return [ map { $p->[$_] * $_ } 1..$#$p ]; }
# ---------- Laurent series { v => lowest exponent, c => [..] } ----------
sub L { my ($v, $c) = @_; return { v => $v, c => $c }; }
sub Ltop { return $_[0]{v} + $#{$_[0]{c}}; }
sub Lc { my ($A, $e) = @_; my $i = $e - $A->{v}; return ($i < 0 || $i > $#{$A->{c}}) ? $Z : $A->{c}[$i]; }
sub Lneg { my $A = shift; return L($A->{v}, [ map { (($A->{v} + $_) % 2) ? -$A->{c}[$_] : $A->{c}[$_]->copy } 0..$#{$A->{c}} ]); }

sub setup { my %o = @_; my $t0 = time;
  $tau = R($o{tau}); $LMAX = $o{lmax}; $DUAL = $o{dual} // 0;
  die "tau = 0, -1 excluded" if $tau->is_zero || $tau == -1;
  $a = R(1) / $tau; $c2 = $tau ** 3 / (2 * (1 + $tau)); $sigma = $tau * (1 + $tau);
  %PAB = (); %PC = (); %W = ();
  $NU = 4 * $LMAX + 10;
  # Taylor coefficients of X(a+u) - X(a) = Sum_{j>=2} X_j u^j,  X_j = ((-1)^(j-1)/j) (a^-j - (tau+1)(1+a)^-j)
  my @X = (R(0), R(0)); for my $j (2..$NU+3) { push @X, R((-1) ** ($j - 1)) / $j * ($a ** (-$j) - (1 + $tau) * (1 + $a) ** (-$j)); }
  die "X_2 != -c2" unless $X[2] == -$c2;
  # zeta = u s(u), s = (1 + q)^(1/2), q_i = X_{i+2}/X_2
  my @q = (R(0)); push @q, $X[$_ + 2] / $X[2] for 1..$NU+1;
  # Lagrange: [zeta^n] u = (1/n) [u^(n-1)] s^(-n)
  @u = (R(0));
  for my $n (1..$NU+2) { my $f = ppow(\@q, R(-$n) / 2, $n - 1); push @u, $f->[$n-1] / $n; }
  # check: X(a+u(zeta)) - X(a) + c2 zeta^2 = O(zeta^(NC+1)), NC = min(NU+2, 30)
  { my $NC = $NU + 2; $NC = 30 if $NC > 30; my $acc = [ map { R(0) } 0..$NC ]; my $pw = [ @u[0..$NC] ]; my $cur = [ R(1) ];
    for my $j (1..$NC) { $cur = pmul($cur, $pw, $NC); next if $j < 2; $acc->[$_] += $X[$j] * $cur->[$_] for 0..$NC; }
    $acc->[2] += $c2; for my $i (0..$NC) { die "u(zeta) wrong at order $i" unless $acc->[$i]->is_zero; } }
  @up = @{ pder(\@u) };
  printf STDERR "  [R2TR] u(zeta) to order %d, %.1f s\n", $NU + 2, time - $t0 if $VERB;
  # r = u/zeta; powers r^m, m = -(2LMAX+2) .. 2LMAX
  my @r = map { $u[$_ + 1] } 0..$NU+1; my @rq = (R(0), @r[1..$#r]);
  my %RP; for my $m (-(2*$LMAX + 2) .. 2*$LMAX) { $RP{$m} = ppow(\@rq, R($m), $NU); }
  # global coefficients A_{l,j} = [zeta^(2l)] (u' u^j),  u^j = zeta^j r^j
  @Aglob = ();
  for my $l (0..$LMAX) { for my $j (0..2*$l) { my $s = R(0);
      for my $i (0..2*$l - $j) { my $rc = $RP{$j}[$i] // $Z; next if $rc->is_zero; $s += $up[2*$l - $j - $i] * $rc; }
      $Aglob[$l][$j] = $s; } }
  # local expansions e_l(zeta) = E_l/dzeta = (1/(2l+1)) Sum_j (j+1) A_{l,j} u'(zeta) zeta^(-j-2) r^(-j-2), up to exponent EM
  my $EM = 2 * $LMAX + 4; @e = (); @eev = (); my %UPR;
  for my $l (0..$LMAX) { my $lo = -2*$l - 2; my @c = map { R(0) } $lo..$EM;
    for my $j (0..2*$l) { my $g = (($j + 1) * $Aglob[$l][$j]) / (2*$l + 1); next if $g->is_zero;
      my $pr = $UPR{$j} //= pmul(\@up, $RP{-$j-2}, $EM + 2*$LMAX + 2);     # u' r^(-j-2) (cached), times zeta^(-j-2)
      for my $i (0..$#$pr) { my $ex = $i - $j - 2; next if $ex < $lo || $ex > $EM; $c[$ex - $lo] += $g * $pr->[$i]; } }
    my $E = L($lo, \@c);
    for my $ex ($lo..-1) { my $want = ($ex == $lo) ? R(1) : R(0); die "E_$l: wrong principal part at zeta^$ex" unless Lc($E, $ex) == $want; }
    push @e, $E;
    my $En = Lneg($E); push @eev, L($lo, [ map { (Lc($E, $lo + $_) + Lc($En, $lo + $_)) / 2 } 0..$EM-$lo ]); }
  printf STDERR "  [R2TR] basis E_l, l <= %d, %.1f s\n", $LMAX, time - $t0 if $VERB;
  # Laplace basis: xi_0 = -1/tau - (1/tau^2) zeta^(-1) r^(-1);  xi_{k+1} = -(1/(2 c2 zeta)) d xi_k/dzeta;  T_{kl} = [zeta^(-2l-2)] d xi_k/dzeta
  { my $rinv = $RP{-1}; my $lo = -1; my @c = map { R(0) } 0..$NU;   # exponents -1 .. NU-1
    for my $i (0..$NU) { $c[$i] = -$rinv->[$i] / ($tau * $tau); } $c[1] += -R(1) / $tau;   # exponent 0 is index 1
    my $xi = L(-1, \@c); @T = ();
    for my $k (0..$LMAX) {
      my $dxi = L($xi->{v} - 1, [ map { ($xi->{v} + $_) * $xi->{c}[$_] } 0..$#{$xi->{c}} ]);
      for my $ex ($dxi->{v}..-1) { next unless ($ex % 2); die "d xi_$k has odd principal part" unless Lc($dxi, $ex)->is_zero; }
      for my $l (0..$LMAX) { $T[$k][$l] = Lc($dxi, -2*$l - 2); }
      die "T not triangular" if grep { !$T[$k][$_]->is_zero } ($k+1..$LMAX);
      my $f = R(-1) / (2 * $c2);
      $xi = L($dxi->{v} - 1, [ map { $f * $_ } @{$dxi->{c}} ]); } }
  @Ti = (); for my $k (0..$LMAX) { for my $l (0..$LMAX) { $Ti[$k][$l] = R(0); } }
  for my $k (0..$LMAX) { $Ti[$k][$k] = R(1) / $T[$k][$k];
    for my $l (reverse 0..$k-1) { my $s = R(0); for my $m ($l+1..$k) { $s += $Ti[$k][$m] * $T[$m][$l]; } $Ti[$k][$l] = -$s / $T[$l][$l]; } }
  # F11 = B(q, sigma q)/dzeta^2 = -u'(z)u'(-z)/(u(z)-u(-z))^2  (Laurent from z^-2)
  { my $N = 2*$LMAX + 6; my @uo = map { ($_ % 2) ? $u[$_]->copy : R(0) } 0..$N+2;   # odd part of u
    my @w = map { $uo[$_ + 1] } 0..$N; # (u(z)-u(-z))/(2z) = Sum uo_{i+1} z^i
    my @upm = map { ($_ % 2) ? -$up[$_] : $up[$_]->copy } 0..$N;
    my $num = pmul([@up[0..$N]], \@upm, $N); my $den = pinv(pmul(\@w, \@w, $N), $N); my $f = pmul($num, $den, $N);
    @F11 = map { -$_ / 4 } @$f; }   # F11 = z^-2 * Sum F11[i] z^i
  return 1; }

# y(zeta) for y = base (-1/w), bkmp (log(1+w)), mix (log(1+w) + e(-1/w)); returns power series (constant term irrelevant)
sub yseries { my ($case, $e) = @_; my $N = 4 * $LMAX + 8;
  my @uu = @u[0..$N+1];
  # Y0 = -1/(a+u) = -tau/(1 + tau u)
  my $den = pinv([ R(1), map { $tau * $uu[$_] } 1..$N+1 ], $N + 1); my @Y0 = map { -$tau * $_ } @$den;
  # YB = log(1 + u/(1+a)) via derivative
  my @vv = (R(0), map { $uu[$_] / (1 + $a) } 1..$N+1);
  my $dv = pder(\@vv); my $iv = pinv([R(1), @vv[1..$#vv]], $N + 1); my $d = pmul($dv, $iv, $N);
  my @YB = (R(0), map { $d->[$_ - 1] / $_ } 1..$N+1);
  return @Y0 if $case eq 'base'; return @YB if $case eq 'bkmp';
  return map { -$_ } @YB if $case eq 'zhou';
  if ($case eq 'mix') { my $eR = R($e); return map { $YB[$_] + $eR * $Y0[$_] } 0..$N+1; }
  die "unknown y $case"; }
# install y (and, in dual mode, the direction eta: y + eps eta)
sub set_y { my ($y, $eta) = @_; my $N = 4 * $LMAX + 8;
  my @Yh = map { ($_ % 2) ? R(0) : 2 * $y->[$_ + 1] } 0..$N;   # (y(z)-y(-z))/z
  die "dy(a) = 0" if $Yh[0]->is_zero;
  @invYh = @{ pinv(\@Yh, $N) }; @invYhE = ();
  if ($DUAL) { my @Eh = map { ($_ % 2) ? R(0) : 2 * $eta->[$_ + 1] } 0..$N;   # d/deps invYh = -invYh^2 Eh
    my $sq = pmul(\@invYh, \@invYh, $N); my $pr = pmul($sq, \@Eh, $N); @invYhE = map { -$_ } @$pr; }
  %PAB = (); %PC = (); %W = (); }

# ---------- numbers: BigRat, or dual [re, eps] ----------
sub nmul { my ($x, $y) = @_; return $x * $y unless $DUAL;
  my $re = $x->[0] * $y->[0]; my $ep = $x->[0] * $y->[1]; $ep += $x->[1] * $y->[0]; return [$re, $ep]; }
sub nacc { my ($h, $k, $v) = @_;
  if (!$DUAL) { return if $v->is_zero; if (exists $h->{$k}) { $h->{$k} += $v; } else { $h->{$k} = $v->copy; } return; }
  if (exists $h->{$k}) { $h->{$k}[0] += $v->[0]; $h->{$k}[1] += $v->[1]; } else { $h->{$k} = [ $v->[0]->copy, $v->[1]->copy ]; } }
sub nzero { my $x = shift; return $DUAL ? ($x->[0]->is_zero && $x->[1]->is_zero) : $x->is_zero; }
sub coefpair { my ($A, $B, $inv, $ex) = @_;   # [zeta^ex] (A(z) * B(z) * inv(z)), inv a power series
  my $s = R(0);
  for my $i ($A->{v}..Ltop($A)) { my $ac = Lc($A, $i); next if $ac->is_zero;
    for my $j ($B->{v}..Ltop($B)) { my $m = $ex - $i - $j; last if $m < 0; next if $m > $#$inv; my $bc = Lc($B, $j); next if $bc->is_zero;
      next if $inv->[$m]->is_zero; $s += $ac * $bc * $inv->[$m]; } }
  return $s; }
# PAB(l0; l1, l2) = (1/(2 c2)) [zeta^(-2 l0)] ( e_l1(z) e_l2(-z) invYh(z) )
sub PAB { my ($l0, $l1, $l2) = @_; my $key = "$l0,$l1,$l2"; return $PAB{$key} if exists $PAB{$key};
  my $En = Lneg($e[$l2]);
  my $v = coefpair($e[$l1], $En, \@invYh, -2*$l0) / (2 * $c2);
  if ($DUAL) { my $w = coefpair($e[$l1], $En, \@invYhE, -2*$l0) / (2 * $c2); $v = [$v, $w]; }
  $PAB{$key} = $v; return $v; }
# PC(l0, l, l') = ((2l+1)/c2) [zeta^(-2 l0 - 2 l)] ( eev_l'(z) invYh(z) )
sub PC { my ($l0, $l, $lp) = @_; my $key = "$l0,$l,$lp"; return $PC{$key} if exists $PC{$key};
  my $ex = -2*$l0 - 2*$l; my ($s, $t) = (R(0), R(0));
  for my $i ($eev[$lp]{v}..$ex) { my $c = Lc($eev[$lp], $i); next if $c->is_zero; my $m = $ex - $i;
    $s += $c * $invYh[$m] if $m <= $#invYh; $t += $c * $invYhE[$m] if $DUAL && $m <= $#invYhE; }
  my $f = R(2*$l + 1) / $c2; my $v = $s * $f; $v = [$v, $t * $f] if $DUAL;
  $PC{$key} = $v; return $v; }
sub stable { return 2*$_[0] - 2 + $_[1] > 0; }
sub base_cases { %W = ();
  my $c03 = $invYh[0] / $c2; $c03 = [$c03, $invYhE[0] / $c2] if $DUAL; $W{"0,3"} = { "0,0,0" => $c03 };
  my %w11;
  for my $l0 (0, 1) { my ($s, $t) = (R(0), R(0)); my $ex = -2*$l0;   # [z^ex] F11 invYh, F11 = z^-2 Sum F11[i] z^i
    for my $i (0..$#F11) { my $m = $ex + 2 - $i; next if $m < 0; $s += $F11[$i] * $invYh[$m]; $t += $F11[$i] * $invYhE[$m] if $DUAL; }
    my $v = -$s / (2 * $c2); $v = [$v, -$t / (2 * $c2)] if $DUAL; nacc(\%w11, "$l0", $v); }
  $W{"1,1"} = \%w11; }
sub step { my ($g, $np1) = @_; my $n = $np1 - 1; my %res;
  # (T1) omega_{g-1,n+2}(q, sigma q, J)
  if ($g >= 1 && stable($g - 1, $n + 2)) { my $Wa = $W{($g-1) . "," . ($n+2)};
    for my $key (keys %$Wa) { my ($l1, $l2, @L) = split /,/, $key; my $cf = $Wa->{$key}; my $rest = join(",", @L);
      for my $l0 (0..$l1 + $l2 + 2) { my $p = PAB($l0, $l1, $l2); next if nzero($p);
        nacc(\%res, ($rest eq '' ? "$l0" : "$l0,$rest"), nmul($cf, $p)); } } }
  # (T2) stable products
  for my $h (0..$g) { for my $mask (0..(2**$n) - 1) {
      my @I = grep { ($mask >> $_) & 1 } 0..$n-1; my @Ic = grep { !(($mask >> $_) & 1) } 0..$n-1;
      next unless stable($h, 1 + @I) && stable($g - $h, 1 + @Ic);
      my $W1 = $W{"$h," . (1 + @I)}; my $W2 = $W{($g - $h) . "," . (1 + @Ic)};
      for my $k1 (keys %$W1) { my ($l1, @p) = split /,/, $k1; my $c1 = $W1->{$k1};
        for my $k2 (keys %$W2) { my ($l2, @qq) = split /,/, $k2; my $cc = nmul($c1, $W2->{$k2});
          my @LJ; my ($ia, $ib) = (0, 0); for my $j (0..$n-1) { push @LJ, (($mask >> $j) & 1) ? $p[$ia++] : $qq[$ib++]; }
          my $rest = join(",", @LJ);
          for my $l0 (0..$l1 + $l2 + 2) { my $pp = PAB($l0, $l1, $l2); next if nzero($pp);
            nacc(\%res, ($rest eq '' ? "$l0" : "$l0,$rest"), nmul($cc, $pp)); } } } } }
  # (T3) B(q,p_j) omega_{g,n}(sigma q, J-j) + omega_{g,n}(q, J-j) B(sigma q, p_j)
  if ($n >= 1 && stable($g, $n)) { my $Wc = $W{"$g,$n"};
    for my $jpos (0..$n-1) { for my $key (keys %$Wc) { my ($lp, @L) = split /,/, $key; my $cf = $Wc->{$key};
        for my $l0 (0..$lp + 1) { for my $l (0..$lp + 1 - $l0) { my $pc = PC($l0, $l, $lp); next if nzero($pc);
            my @LL = @L; splice(@LL, $jpos, 0, $l); nacc(\%res, join(",", $l0, @LL), nmul($cf, $pc)); } } } } }
  for my $k (keys %res) { delete $res{$k} if nzero($res{$k}); }
  $W{"$g,$np1"} = \%res; }
sub run { my ($S, $gmax) = @_; base_cases(); my @done;
  for my $chi (2..3*$S) { for my $g (0..$gmax) { my $n = $chi - 2*$g + 2; next if $n < 1 || $g + $n > $S;
      next if ($g == 0 && $n <= 3) || ($g == 1 && $n == 1); my $t = time; step($g, $n); push @done, [$g, $n];
      printf STDERR "  [R2TR] (%d,%d) %d terms %.1f s\n", $g, $n, scalar(keys %{$W{"$g,$n"}}), time - $t if $VERB; } }
  return @done; }
# conversion E-basis -> xi-basis: c^xi[K] = Sum_L c^E[L] prod Ti[L_i][K_i]   (real part, or a chosen part in dual mode)
sub to_xi { my ($Wx, $part) = @_; my %out;
  for my $key (keys %$Wx) { my @Lv = split /,/, $key; my $cf = $DUAL ? $Wx->{$key}[$part // 0] : $Wx->{$key}; next if $cf->is_zero;
    my @acc = ([[], $cf]);
    for my $i (0..$#Lv) { my @new; for my $o (@acc) { for my $k (0..$Lv[$i]) { my $t = $Ti[$Lv[$i]][$k]; next if $t->is_zero; push @new, [[@{$o->[0]}, $k], $o->[1] * $t]; } } @acc = @new; }
    for my $o (@acc) { my $k = join(",", @{$o->[0]}); if (exists $out{$k}) { $out{$k} += $o->[1]; } else { $out{$k} = $o->[1]->copy; } } }
  for (keys %out) { delete $out{$_} if $out{$_}->is_zero; }
  return \%out; }
# x-expansion: P_l(m) = [x^m] Prim(E_l) = (1/m) [w^(m-1)] (E_l/dw) (1+w)^(m(tau+1))   (Lagrange, x = w (1+w)^(-tau-1))
our %PL;
sub binR { my ($b, $i) = @_; my $r = R(1); for my $t (0..$i-1) { $r *= ($b - $t) / ($t + 1); } return $r; }
sub binI { my ($n, $k) = @_; my $r = R(1); for my $t (1..$k) { $r = $r * ($n - $k + $t) / $t; } return $r; }
sub PLm { my ($l, $m) = @_; my $key = "$l,$m"; return $PL{$key} if exists $PL{$key};
  my $beta = $m * (1 + $tau); my $s = R(0);
  for my $j (0..2*$l) { my $g = (($j + 1) * $Aglob[$l][$j]) / (2*$l + 1); next if $g->is_zero;
    for my $nn (0..$m-1) { my $i = $m - 1 - $nn;
      $s += $g * R((-1) ** $j) * binI($nn + $j + 1, $j + 1) * $tau ** ($j + 2 + $nn) * binR($beta, $i); } }
  $s /= $m; $PL{$key} = $s; return $s; }
sub Pi_mu { my $mu = shift; my $r = R(1); for my $t (1..$mu-1) { $r *= ($mu * $tau + $t) / $t; } return $r; }
1;
