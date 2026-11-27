package rk_pert;
# rk.pm -- referee 2 (route R1, round 2): independent modular library.
# Written from the statements of Route A Theorem D (all R) and leaf_structure_proof.md Theorem 11.1 only; shares no code
# with the author's kpmod.pm / kp_point.pl / r1lib.pm.
#
# Conventions (brief s.1, c = 1): sigma = tau(1+tau), Pi_j = prod_{a=1}^{j-1}(j tau + a)/(j-1)!,
# Y(T) := Pi_T Xi_G(T) (T a multiset of parts >= 2), y_L := Y(L+1) (one-leaf value at level L), unit parts by string
# Y(M + {1}) = J(M) Y(M), Y(1^k) = 0.
# Genus-0 kernel: exp(D(z)D(w)F_0) = (Y_P(z)-Y_P(w))/(z-w) = 1 + sum_{i,j>=1} z^-i w^-j d_1 d_{i+j-1} F_0, and
#   d^Q d_1 d_n F_0 (0) = sigma Pi_n Pi^Q Xi_0(Q + {1,n}),  Xi_0(T) = Gamma(J)/Gamma(J-r+3)  (unit parts included).
# Theorem D at [P^R z^-al w^-be] (T = R + {al,be}):
#   Y(T) = (J(R)+al+be-1) Y(R + {al+be-1}) - sum_{a<al,b<be} sum_{S sub R} w(S) Y(S+{a,b}) E0(al-a+be-b-1; R-S),
#   w(S) = prod_j binom(r_j, s_j)  (labelled Leibniz rule).
# The propagation is G-independent; genus G only fixes the range of validity (levels <= 2G-3).
use strict; use warnings;
use Exporter 'import';
our @EXPORT = qw($P md mm mi mfrac mpow setp gbin binom_int fct ifct
                 parts_level all_T keyof Jof Dof mults submultisets
                 kp_setup leaf_series ser_mul ser_inv ser_compose ser_log ser_exp T_series hyp_series kappa_list);
our $P = 1000000007;
sub setp { $P = shift; }
sub md { my $x = shift; $x %= $P; $x += $P if $x < 0; return $x; }
sub mm { return ($_[0] * $_[1]) % $P; }
sub mpow { my ($b, $e) = @_; $b = md($b); my $r = 1; while ($e > 0) { $r = mm($r, $b) if $e & 1; $b = mm($b, $b); $e >>= 1; } return $r; }
sub mi { my $a = md(shift); die "mi: zero\n" unless $a; return mpow($a, $P - 2); }
sub mfrac { my ($n, $d) = @_; return mm(md($n), mi($d)); }
my @F = (1); sub fct { my $n = shift; for my $i (scalar(@F)..$n) { $F[$i] = mm($F[$i-1], $i); } return $F[$n]; }
sub ifct { return mi(fct(shift)); }
sub binom_int { my ($n, $k) = @_; return 0 if $k < 0 || $k > $n; return mm(fct($n), mi(mm(fct($k), fct($n-$k)))); }
# generalised binomial binom(x, d) for x mod p, integer d >= 0
sub gbin { my ($x, $d) = @_; return 0 if $d < 0; my $v = 1; for my $i (0..$d-1) { $v = mm($v, md($x - $i)); } return mm($v, ifct($d)); }

# ---------- multisets ----------
sub keyof { return join(',', sort { $b <=> $a } @_); }
sub Jof { my $s = 0; $s += $_ for @_; return $s; }
sub Dof { my $s = 0; $s += $_ - 1 for @_; return $s; }
sub mults { my %m; $m{$_}++ for @_; return \%m; }
# partitions of d into parts <= mx (descending)
my %PL; sub parts_level { my ($d, $mx) = @_; return ([]) if $d == 0; my $k = "$d,$mx";
  return @{ $PL{$k} //= do { my @o; for my $a (reverse 1..($mx < $d ? $mx : $d)) { push @o, [$a, @$_] for parts_level($d - $a, $a); } \@o } }; }
# all multisets T of parts >= 2 with 1 <= level <= K (descending arrays)
sub all_T { my $K = shift; my @T; for my $d (1..$K) { for my $q (parts_level($d, $d)) { push @T, [ map { $_ + 1 } @$q ]; } } return \@T; }
# sub-multisets of a multiset (array): list of [S, R-S, weight]
my %SM; sub submultisets { my @R = @_; my $k = keyof(@R); return @{ $SM{$k} //= do {
  my $m = mults(@R); my @v = sort { $b <=> $a } keys %$m; my @out = ([[], [], 1]);
  for my $x (@v) { my $n = $m->{$x}; my @new; for my $o (@out) { for my $i (0..$n) { push @new, [[@{$o->[0]}, ($x) x $i], [@{$o->[1]}, ($x) x ($n - $i)], mm($o->[2], binom_int($n, $i))]; } } @out = @new; }
  \@out } }; }

# ---------- truncated power series (arrays, index = power) ----------
sub ser_mul { my ($a, $b, $N) = @_; my @c = (0) x ($N+1); for my $i (0..$N) { next if $i > $#$a || !$a->[$i]; for my $j (0..$N-$i) { last if $j > $#$b; next unless $b->[$j]; $c[$i+$j] = ($c[$i+$j] + $a->[$i]*$b->[$j]) % $P; } } return \@c; }
sub ser_inv { my ($a, $N) = @_; my $i0 = mi($a->[0]); my @b = ($i0); for my $n (1..$N) { my $s = 0; for my $k (1..$n) { last if $k > $#$a; $s = ($s + $a->[$k]*$b[$n-$k]) % $P; } $b[$n] = md(-mm($s, $i0)); } return \@b; }
sub ser_compose { my ($f, $g, $N) = @_; die "compose: g(0) != 0\n" if $g->[0]; my @r = (0) x ($N+1); my $pw = [1, (0) x $N];
  for my $k (0..$N) { last if $k > $#$f; if ($f->[$k]) { for my $i (0..$N) { $r[$i] = ($r[$i] + $f->[$k]*$pw->[$i]) % $P; } } $pw = ser_mul($pw, $g, $N); } return \@r; }
sub ser_log { my ($a, $N) = @_; die "log: a0 != 1\n" unless $a->[0] == 1; my $ia = ser_inv($a, $N); my @da = map { $_+1 <= $#$a ? mm($_+1, $a->[$_+1]) : 0 } 0..$N;
  my $q = ser_mul(\@da, $ia, $N); my @b = (0); for my $n (1..$N) { $b[$n] = mm($q->[$n-1], mi($n)); } return \@b; }
sub ser_exp { my ($a, $N) = @_; die "exp: a0 != 0\n" if $a->[0]; my @b = (1); for my $n (1..$N) { my $s = 0; for my $k (1..$n) { last if $k > $#$a; $s = ($s + mm($k, $a->[$k])*$b[$n-$k]) % $P; } $b[$n] = mm($s, mi($n)); } return \@b; }
# T(t) = 1 - e^-t
sub T_series { my $N = shift; my @T = (0); for my $n (1..$N) { $T[$n] = mm(($n % 2) ? 1 : $P - 1, ifct($n)); } return \@T; }
# 2F1(a,b;c;x) as a series in x (a, b, c mod p; c not a non-positive integer)
sub hyp_series { my ($a, $b, $c, $N) = @_; my @s = (1); my $t = 1; for my $n (1..$N) { $t = mm($t, mm(mm(md($a+$n-1), md($b+$n-1)), mi(mm(md($c+$n-1), $n)))); push @s, $t; } return \@s; }
# kappa_k = [x^k] (tau x + log((e^x-1)/x)), k = 0..N
sub kappa_list { my ($tau, $N) = @_; my @e = map { ifct($_ + 1) } 0..$N; my $l = ser_log(\@e, $N); $l->[1] = md($l->[1] + $tau); return $l; }

# ---------- leaf curve (Theorem 11.1 + string elimination of lam_1) ----------
# P_j(t) = (-1)^j (Pi_j/j) t^(j-1) H_j(T)/H_1(T)^j,  H_j = 2F1(-j tau, -j tau; j; T),  T = 1 - e^-t;  returns {j => series to t^N}
sub leaf_series { my ($tau, $Pi, $N, $jmax) = @_; my $T = T_series($N);
  my $H1 = ser_compose(hyp_series(md(-$tau), md(-$tau), 1, $N), $T, $N); $H1 = [1, (0) x $N] if ($ENV{PERT} // "") eq "noH1"; if (($ENV{PERT} // "") eq "Tt") { $T = [0, 1, (0) x ($N-1)]; $H1 = ser_compose(hyp_series(md(-$tau), md(-$tau), 1, $N), $T, $N); } my $iH1 = ser_inv($H1, $N);
  my %Pj; my $pw = [1, (0) x $N];
  for my $j (1..$jmax) { $pw = ser_mul($pw, $iH1, $N); next if $j < 2; my $Hj = ser_compose(hyp_series(md(-$j*$tau), md(-$j*$tau), $j, $N), $T, $N);
    my $r = ser_mul($Hj, $pw, $N); my $c = mm(($j % 2) ? $P - 1 : 1, mm($Pi->[$j], mi($j))); my @s = (0) x ($N+1);
    for my $o (0..$N-$j+1) { $s[$o+$j-1] = mm($c, $r->[$o]); } $Pj{$j} = \@s; }
  return \%Pj; }

# ---------- KP propagation ----------
# kp_setup(tau, K, yvals, strategy): yvals->[L] = arrayref (length m) for L = 1..K (one-leaf values y_L); returns {Y => {key => arrayref}, ...}
# strategy: 'big' (pair = two largest parts), 'small' (two smallest), 'mix' (largest and smallest)
sub kp_setup { my ($tau, $K, $yv, $strat) = @_; $strat //= 'big'; my $m = scalar @{$yv->[1]};
  my $sig = mm($tau, md($tau + 1));
  my @Pi = (1, 1); for my $j (2..2*$K+4) { my $v = 1; for my $a (1..$j-1) { $v = mm($v, md($j*$tau + $a)); } $Pi[$j] = mm($v, ifct($j-1)); }
  my %Y; my $zero = [ (0) x $m ];
  for my $L (1..$K) { $Y{$L+1} = $yv->[$L]; }
  my %E0C; my $E0 = sub { my ($n, $Q) = @_; my $k = "$n;" . join(',', @$Q); return $E0C{$k} //= do {
      my $v = mm($sig, $Pi[$n]); $v = mm($v, 2) if ($ENV{PERT} // "") eq "E0" && $n == 3 && @$Q == 1; $v = mm($v, $Pi[$_]) for @$Q;
      if (!@$Q) { $v = mm($v, mi($n + 1)); } else { my $Jp = Jof(@$Q) + $n + 1; for my $i (1..scalar(@$Q)-1) { $v = mm($v, md($Jp - $i)); } }
      $v } };
  my $Yg = sub { my @M = @_; my @Tp = sort { $b <=> $a } grep { $_ >= 2 } @M; my $u = @M - @Tp; return $zero unless @Tp;
    my $v = $Y{join(',', @Tp)}; die "Y(@Tp) not yet computed\n" unless defined $v; return $v if $u == 0;
    my $J = Jof(@Tp); my $f = 1; $f = mm($f, $J + $_) for 0..$u-1; return [ map { mm($_, $f) } @$v ]; };
  my $Ts = all_T($K); my @order = sort { Jof(@$a) <=> Jof(@$b) } grep { @$_ >= 2 } @$Ts;
  for my $T (@order) { my @t = @$T; my ($al, $be, @R);
    if ($strat eq 'big') { ($be, $al, @R) = @t; }
    elsif ($strat eq 'small') { $al = pop @t; $be = pop @t; @R = @t; }
    else { $be = shift @t; $al = pop @t; @R = @t; }
    my $JR = Jof(@R); my @acc = (0) x $m;
    my $add = sub { my ($c, $v) = @_; return unless $c; for my $i (0..$m-1) { $acc[$i] = ($acc[$i] + $c*$v->[$i]) % $P if $v->[$i]; } };
    $add->(md($JR + $al + $be - 1), $Yg->(@R, $al + $be - 1));
    for my $sb (submultisets(@R)) { my ($S, $Q, $w) = @$sb;
      for my $a (1..$al-1) { for my $b (1..$be-1) { my $n = ($al - $a) + ($be - $b) - 1; my $e = $E0->($n, $Q); next unless $e;
          $add->(md(-mm($w, $e)), $Yg->(@$S, $a, $b)); } } }
    $Y{keyof(@$T)} = \@acc; }
  return { tau => $tau, sig => $sig, Pi => \@Pi, Y => \%Y, K => $K, m => $m, Ts => $Ts }; }
1;
